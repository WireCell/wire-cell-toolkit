/**
   custard filters

   Provide boost::iostreams input and output filters to convert
   between custard and tar codec.
 */

#ifndef CUSTARD_BOOST_HPP
#define CUSTARD_BOOST_HPP

#include "custard.hpp"
#include "custard_stream.hpp"

#ifdef CUSTARD_BOOST_USE_MINIZ
#define MINIZ_NO_ZLIB_APIS
#include "miniz.h"
#endif

#include <boost/iostreams/concepts.hpp>    // multichar_output_filter
#include <boost/iostreams/operations.hpp> // write
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
// Optional: throws link errors on some systems
#ifdef CUSTARD_BOOST_USE_LZMA
#include <boost/iostreams/filter/lzma.hpp>
#endif
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include <boost/filesystem.hpp>
// Optional: boost-interal compiler warning
#ifdef CUSTARD_BOOST_USE_PROC
#include <boost/process.hpp>
#endif
// #define BOOST_SPIRIT_DEBUG
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/adapted/std_pair.hpp>

#include <vector>

#include <regex>
#include <sstream>

#include <iostream> // debug

namespace custard {

    namespace qi = boost::spirit::qi;

    using buffer_t = std::vector<char>;

    template <typename Iterator>
    struct parser : qi::grammar<Iterator, dictionary_t()>
    {
        using keyval_t = std::pair<std::string, std::string>;

        parser() : parser::base_type(header)
        {
            header = qi::no_skip[ ini >> *opt >> fin ];
            ini = qi::no_skip[qi::string("name") >> ' ' >> value >> '\n'];
            opt = qi::no_skip[keyword >> ' ' >> value >> '\n'];
            fin = qi::no_skip[qi::string("body") >> ' ' >> value >> '\n'];
            keyword = qi::no_skip[ qi::string("mode") | qi::string("mtime") | qi::string("uid") | qi::string("gid") | qi::string("uname") | qi::string("gname") ];
            value = qi::no_skip[+~qi::lit('\n')];

            BOOST_SPIRIT_DEBUG_NODE(header);
            BOOST_SPIRIT_DEBUG_NODE(ini);
            BOOST_SPIRIT_DEBUG_NODE(opt);
            BOOST_SPIRIT_DEBUG_NODE(fin);
            BOOST_SPIRIT_DEBUG_NODE(keyword);
            BOOST_SPIRIT_DEBUG_NODE(value);
        }

      private:
        qi::rule<Iterator, dictionary_t()> header;
        qi::rule<Iterator, keyval_t()> ini, opt, fin;
        qi::rule<Iterator, std::string()> keyword, value;
    };    

    template <typename Iterator>
    bool parse_vars(Iterator &first, Iterator last, dictionary_t & vars)
    {
        parser<Iterator> p;
        return boost::spirit::qi::parse(first, last, p, vars);
    }

    class stream_parser {
      public:

        using parsed_t = std::vector<std::pair<dictionary_t, buffer_t>>;

        // Call with new data
        parsed_t feed(const char* s, std::streamsize n)
        {
            // always slurp everything
            buf.insert(buf.end(), s, s+n);
            while (proc()) { /* empty */ }
            return std::move(collected);
        }

      private:
        // return true if more on the plate, false if starved.
        // Anything which is done gets collected.
        bool proc()
        {
            if (state == State::inhead) {
                buffer_t::const_iterator beg(buf.begin());
                buffer_t::const_iterator end(buf.end());
                params.clear();
                bool ok = custard::parse_vars(beg, end, params);
                if (!ok) {
                    return false;
                }
                body_size = lval(params, "body");
                buf = buffer_t(beg, end); // more for later
                state = State::inbody;
                return true;
            }

            if (state == State::inbody) {
                if (buf.size() < body_size) {
                    return false; // keep slurping body
                }
                auto mid = buf.begin()+body_size;
                buffer_t keep(mid, buf.end());
                buf.resize(mid-buf.begin());
                collect(std::move(buf));
                buf = std::move(keep);
                state = State::inhead;
                return true;
            }
            return false;
        }

        void collect(buffer_t&& give)
        {
            collected.push_back(std::make_pair(params, std::move(give)));
            body_size = 0;
            params.clear();
        }

      private:

        buffer_t buf;           // buffer fed data
        dictionary_t params;    // current member
        size_t body_size{0};    // current member

        enum class State { inhead, inbody };
        State state{State::inhead};

        parsed_t collected;     // temporary return
    };

    
#ifdef CUSTARD_BOOST_USE_MINIZ
    class miniz_sink {
      public:
        typedef char char_type;
        struct category :
            public boost::iostreams::sink_tag,
            public boost::iostreams::closable_tag
        { };
        
        miniz_sink(const std::string& outname)
            : zip(std::make_shared<mz_zip_archive>())
        {
            memset(zip.get(), 0, sizeof(mz_zip_archive));
            if (!mz_zip_writer_init_file(zip.get(), outname.c_str(), 0)) {
                //std::cerr << "miniz_sink: failed to initialize miniz for " << outname << "\n";
                throw std::runtime_error("failed to initialize miniz writing to " + outname);
            }
        }

        std::streamsize write(const char* s, std::streamsize n)
        {
            auto got = parser.feed(s, n);
            for (auto& [params, buf] : got) {
                auto name = sval(params, "name");

                auto ok = mz_zip_writer_add_mem(
                    zip.get(), name.c_str(), buf.data(), buf.size(),
                    MZ_BEST_SPEED);
                if (!ok) {
                    //std::cerr << "miniz_sink: failed to write " << name << "\n";
                    throw std::runtime_error("failed to write " + name);
                }
                // std::cerr << "miniz_sink: write: " << name << " ["<<buf.size()<<"]\n";
            }
            return n;
        }

        void close()
        {
            mz_zip_writer_finalize_archive(zip.get());
            mz_zip_writer_end(zip.get());
            zip = nullptr;
        }


      private:
        // this class must be copyable so put non-copyable into shared ptr
        stream_parser parser;
        std::shared_ptr<mz_zip_archive> zip;
    };

    class miniz_source {
      public:
        typedef char char_type;
        struct category :
            public boost::iostreams::source_tag,
            public boost::iostreams::closable_tag
        { };

        miniz_source(const std::string& inname)
            : zip(std::make_shared<mz_zip_archive>())
        {
            memset(zip.get(), 0, sizeof(mz_zip_archive));
            if (!mz_zip_reader_init_file(zip.get(), inname.c_str(), 0)) {
                //std::cerr << "miniz_sink: failed to initialize miniz for " << inname << "\n";
                auto err = mz_zip_get_last_error(zip.get());
                std::string errs = mz_zip_get_error_string(err);
                throw std::runtime_error("failed to initialize miniz reading from " + inname + ": " + errs);
            }
            memnum = mz_zip_reader_get_num_files(zip.get());
            memind = 0;
            memptr = member.end();
        }

        // Read up to n, return number read or -1 for EOF.  Will
        // return a "short read" when we reach end of a member.
        std::streamsize read(char* s, std::streamsize n)
        {
            if (memptr == member.end()) {
                next_member();
            }
            if (memptr == member.end()) {
                // still then at EOF
                return -1;
            }

            std::streamsize left = member.end() - memptr;
            std::streamsize take = std::min(left, n);
            std::copy(memptr, memptr+take, s);
            memptr += take;
            return take;            
        }

        void close()
        {
            member.clear();
            memptr = member.end();
            mz_zip_reader_end(zip.get());
            zip = nullptr;
        }

      private:
        std::shared_ptr<mz_zip_archive> zip;
        size_t memnum{0}, memind{0};
        std::vector<char> member;
        std::vector<char>::iterator memptr{member.end()};

        // Read in the next member, return its index.
        size_t next_member() {
            if (memind >= memnum) {
                member.clear();
                memptr = member.end();
                return memnum;
            }

            const size_t index = memind;
            ++memind;

            mz_zip_archive_file_stat fs;
            if (!mz_zip_reader_file_stat(zip.get(), index, &fs)) {
                mz_zip_reader_end(zip.get());
                throw std::runtime_error("failed to stat miniz member");
            }
            bool isdir = mz_zip_reader_is_file_a_directory(zip.get(), index);
            if (isdir) {
                return next_member();
            }
            const size_t bsize = fs.m_uncomp_size;

            std::stringstream ss;
            ss << "name " << fs.m_filename << "\nbody " << bsize << "\n";
            const std::string head = ss.str();
            const size_t hsize = head.size();

            member.resize(hsize + bsize, 0);
            memptr = member.begin();

            std::copy(head.begin(), head.end(), member.begin());
            auto ok = mz_zip_reader_extract_to_mem(
                zip.get(), index, member.data()+hsize, bsize, 0);
            if (!ok) {
                mz_zip_reader_end(zip.get());
                throw std::runtime_error("failed to read miniz member");
            }
            return index;
        }
    };

#endif  // miniz

#ifdef CUSTARD_BOOST_USE_PROC
    class proc_sink {
      public:
        typedef char char_type;
        struct category :
            public boost::iostreams::sink_tag,
            public boost::iostreams::closable_tag
        { };
        proc_sink(const std::string& cmd, const std::vector<std::string>& args)
        {
            p = std::make_shared<impl>(cmd, args);
        }
        void close()
        {
            p->in.close();
            p->in.pipe().close();
            p->child.wait();
        }
        std::streamsize write(const char* s, std::streamsize n)
        {
            if (!p->child.running()) {
                return 0;
            }
            if (!p->in) {
                return 0;
            }
            p->in.write(s, n);
            if (!p->in) {
                return 0;
            }
            return n;
        }
      private:
        struct impl {
            impl(const std::string& cmd, const std::vector<std::string>& args)
                : in{}, child(boost::process::search_path(cmd),
                              boost::process::args(args),
                              boost::process::std_out > stdout,
                              boost::process::std_err > stderr,
                              boost::process::std_in < in)
            {
            }
            boost::process::opstream in;
            boost::process::child child;
        };
        std::shared_ptr<impl> p;
    };
#endif

    inline
    bool assuredir(const std::string& pathname)
    {
        boost::filesystem::path p(pathname);
        if ( ! p.extension().empty() ) {
            p = p.parent_path();
        }
        if (p.empty()) {
            return false;
        }
        return boost::filesystem::create_directories(p);
    }

    // This filter inputs a custard stream and outputs a tar stream.
    class tar_writer : public boost::iostreams::multichar_output_filter {
      public:

        template<typename Sink>
        std::streamsize write(Sink& dest, const char* buf, std::streamsize bufsiz)
        {
            auto got = parser.feed(buf, bufsiz);
            for (auto& [params, buf] : got) {
                auto th = make_header(params);
                // std::cerr << "tar_writer: " << th.name()
                //           << " [" << th.size() << "] = [" << buf.size() << "]\n";
                boost::iostreams::write(dest, (char*)th.as_bytes(), 512);

                boost::iostreams::write(dest, buf.data(), buf.size());
                const size_t extra = buf.size() % 512;
                if (extra) {
                    const size_t pad = 512 - extra;
                    for (size_t ind=0; ind<pad; ++ind) {
                        boost::iostreams::put(dest, '\0');
                    };
                }
            }
            return bufsiz;
        }

        template<typename Sink>
        void close(Sink& dest) {
            buffer_t buf(1024, 0);
            boost::iostreams::write(dest, buf.data(), buf.size());
        }

      private:
        stream_parser parser;

    };


    // This filter inputs a tar stream and produces custard stream.
    class tar_reader : public boost::iostreams::multichar_input_filter {
      public:

        template<typename Source>
        std::streamsize slurp_header(Source& src)
        {
            // std::cerr << "tar_reader: slurp_header\n";
            std::streamsize got{0};

            // Either we have a next non-zero header or we have zero
            // headers until we hit EOF.
            while (true) {
                got = boost::iostreams::read(src, th.as_bytes(), 512);
                if (got < 0) {
                    // std::cerr << "tar_reader: EOF in header\n";
                    return got;
                }
                if (th.size() == 0) {
                    // std::cerr << "tar_reader: empty header\n";
                    continue;
                }
                break;
            }
            assert (got == 512);

            // std::cerr << "tar_reader: |" << th.name() << "| " << th.size() <<"|\n";

            bodyleft = th.size();
            std::ostringstream ss;
            custard::write(ss, th);
            headstring = ss.str();

            // std::cerr << "\theadstring:\n|" << headstring << "|\n";

            return 512;
        }

        template<typename Source>
        std::streamsize read_one(Source& src, char* buf, std::streamsize bufsiz)
        {
            // std::cerr << "tar_reader: read_one at state "
            //           << (int)state << " bufsize "<<bufsiz<<"\n";
            if (state == State::gethead) {
                auto got = slurp_header(src);
                if (got < 0) {
                    // std::cerr << "tar_reader: forwarding EOF in read_one\n";
                    return got;
                }
                state = State::sendhead;
                return read_one(src, buf, bufsiz);
            }

            // The rest may take several calls if bufsiz is too small.

            if (state == State::sendhead) { 
                std::streamsize tosend = std::min<size_t>(bufsiz, headstring.size());
                std::memcpy(buf, headstring.data(), tosend);
                headstring = headstring.substr(tosend);
                if (headstring.empty()) {
                    state = State::body;
                }
                return tosend;
            }

            if (state == State::body) {
                std::streamsize want = std::min<size_t>(bufsiz, bodyleft);
                auto put = boost::iostreams::read(src, buf, want);

                // std::cerr << "bodyleft: " << bodyleft << " " << bufsiz
                //           << " " << want << " " << put << "\n";
                if (put < 0) {
                    return put;
                }

                bodyleft -= put;
                if (bodyleft == 0) {
                    state = State::filedone;
                    // slurp padding; a member whose size is an exact multiple
                    // of 512 has NO padding (jump must be 0, not 512, else the
                    // next member's tar header block is swallowed).
                    const size_t extra = th.size()%512;
                    const size_t jump = extra ? 512 - extra : 0;
                    if (jump) {
                        std::string pad(jump, 0);
                        auto got = boost::iostreams::read(src, &pad[0], jump);
                        // std::cerr << "padleft: " << jump << " " << got << "\n";
                        if (got < 0) {
                            return got;
                        }
                    }
                }
                return put;
            }
            return -1;
        }


        template<typename Source>
        std::streamsize read(Source& src, char* buf, std::streamsize bufsiz)
        {
            std::streamsize filled = 0; // number filled in buf
            char* ptr = buf;
            while (filled < bufsiz) {
                std::streamsize left = bufsiz - filled;
                auto got = read_one(src, ptr, left);
                if (got < 0) {
                    return got;
                }
                if (!got) {
                    return filled;
                }
                filled += got;
                ptr += got;
                if (state == State::filedone) {
                    state = State::gethead;
                    return filled;
                }
            }
            return filled;
        }

      private:
        custard::Header th;
        // temp buffer for custard codec of header
        std::string headstring{""};
        // track how much of the body is left to read
        size_t bodyleft{0};
        enum class State { gethead, sendhead, body, getpadd, filedone };
        State state{State::gethead};
        
    };

    // Build input filtering stream based on parsing the filename.
    inline
    void input_filters(boost::iostreams::filtering_istream& in,
                       std::string inname)
    {
        auto has = [&](const std::string& things) {
            return std::regex_search(inname, std::regex("[_.]("+things+")\\b"));
        };

        if (has("tar|tar.gz|tgz|tar.bz2|tbz|tbz2|tar.xz|txz|tar.pixz|tix|tpxz")) {
            in.push(custard::tar_reader());
        }

        if (has("gz|tgz")) {
            in.push(boost::iostreams::gzip_decompressor());
        }
        else if (has("bz2|tbz|tbz2")) {
            in.push(boost::iostreams::bzip2_decompressor());
        }
#ifdef CUSTARD_BOOST_USE_LZMA
        else if (has("xz|txz|pixz|tix|tpxz")) {
            // pixz makes xz-compatible files
            in.push(boost::iostreams::lzma_decompressor());
        }
#endif
        else if (has("zip|npz")) {
#ifdef CUSTARD_BOOST_USE_MINIZ
            in.push(custard::miniz_source(inname));
            return;
#else
            throw std::runtime_error("zip/npz file support requires compilation with CUSTARD_BOOST_USE_MINIZ");
#endif
        }

        in.push(boost::iostreams::file_source(inname));
    }

    /// Parse outname and based complete the filter ostream.  If
    /// parsing fails, nothing is added to "out".  The compression
    /// level is interpreted as being from 0-9, higher is more,
    /// slower.
    inline
    void output_filters(boost::iostreams::filtering_ostream& out,
                        std::string outname,
                        int level = 1)
    {
        auto has = [&](const std::string& things) {
            return std::regex_search(outname, std::regex("[_.]("+things+")\\b"));
        };

        assuredir(outname);

        // Add tar writer if we see tar at the end.
        if (has("tar|tar.gz|tgz|tar.bz2|tbz|tbz2|tar.xz|txz|tar.pixz|tix|tpxz")) {
            out.push(custard::tar_writer());
        }

        // In addition, if compression is wanted, add the appropriate
        // filter next.
        if (has("gz|tgz")) {
            out.push(boost::iostreams::gzip_compressor(level));
        }
        else if (has("bz2|tbz|tbz2")) {
            out.push(boost::iostreams::bzip2_compressor(level));
        }
#ifdef CUSTARD_BOOST_USE_LZMA
        else if (has("xz|txz")) {
            out.push(boost::iostreams::lzma_compressor(level));
        }
#endif
#ifdef CUSTARD_BOOST_USE_PROC
        else if (has("pixz|tix|tpxz")) {
            // std::stringstream ss;
            // ss << "-p 1 -" << level << " -o " << outname;
            // std::cerr << ss.str() << std::endl;
            out.push(custard::proc_sink("pixz", {"-p", "1", "-1", "-o", outname}));
            return;
        }
#endif
        else if (has("zip|npz")) {
#ifdef CUSTARD_BOOST_USE_MINIZ
            out.push(custard::miniz_sink(outname));
            return;
#else
            throw std::runtime_error("zip/npz file support requires compilation with CUSTARD_BOOST_USE_MINIZ");
#endif
        }


        // finally, we save to the actual file
        out.push(boost::iostreams::file_sink(outname));
    }
}
#endif
