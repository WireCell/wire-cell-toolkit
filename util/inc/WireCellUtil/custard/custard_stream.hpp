/**
   custard streams

   Provide functions to read/write custard header from/to std::stream
   or derived.

   stream = *member
   member = name *option body
   name = "name" SP namestring LF
   body = "body" SP numstring LF data
   data = *OCTET
   option = ("mode" / "mtime" / "uid" / "gid") SP numstring LF
   option /= ("uname" / "gname" ) SP identstring LF

 */
#ifndef CUSTARD_STREAM_HPP
#define CUSTARD_STREAM_HPP

#include "custard.hpp"

namespace custard {

    // inline
    // std::istream& get_line(std::istream& si, std::string& s, char delim='\n')
    // {
    //     while (true) {
    //         char c{0};
    //         si.get(c);
    //         if (!si or c == delim) {
    //             break;
    //         }
    //         s.push_back(c);
    //         std::cerr << "get_line: |" << s << "|\n";
    //     }
    //     return si;
    // }

    inline
    std::string get_string(std::istream& si, char delim=' ')
    {
        std::string s;
        std::getline(si, s, delim);
        if (!si) {
            // std::cerr << "get_string broke stream (delim=\""<<delim<<"\")\n";
            return "";
        }
        return s;
    }

    inline
    size_t get_number(std::istream& si)
    {
        std::string nstr = get_string(si, '\n');
        if (!si) {
            // std::cerr << "get_number broke stream\n";
            return 0;
        }

        return std::stol(nstr);
    }

    // Read header from stream.  Extract characters from istream
    // filling head until reaching the start of the member body data.
    // Caller MUST continue to read exactly the number of bytes
    // indicated by head.size() in order to position stream at the
    // start of the next archive member.
    inline
    std::istream& read(std::istream& si, Header& head)
    {
        head.clear();
        head.init();
        while (true) {
            std::string key = get_string(si, ' ');
            if (!si) {
                // std::cerr << "custard::read(header) error \""
                //           << strerror(errno) 
                //           << "\" in reading key " << key << "\n";
                return si;
            }
            // std::cerr << "custard::read(header) input stream at key |"<<key<<"|\n";
            if (key == "name") {
                head.set_name(get_string(si, '\n'));
                // std::cerr << "\t|" << head.name() << "|\n";
            }
            else if (key == "mode") {
                head.set_mode(get_number(si));
            }
            else if (key == "mtime") {
                head.set_mtime(get_number(si));
            }
            else if (key == "uid") {
                head.set_uid(get_number(si));
            }
            else if (key == "gid") {
                head.set_gid(get_number(si));
            }
            else if (key == "uname") {
                head.set_uname(get_string(si, '\n'));
            }
            else if (key == "gname") {
                head.set_gname(get_string(si, '\n'));
            }
            else if (key == "body") {
                head.set_size(get_number(si));
                head.gensum();
                // std::cerr << "custard::read(header) got body size "
                //           << head.size() << " " << head.chksum() << " " << head.checksum() << std::endl;
                return si;
            }
            if (!si) {
                // std::cerr << "custard::read(header) input stream is bad: " << strerror(errno) << std::endl;
                return si;
            }
        }
        return si;        
    }
    // The minimum header read.
    inline
    std::istream& read(std::istream& si,
                       std::string& filename, size_t& filesize)
    {
        if (!si) {
            // std::cerr << "custard::read given bad stream\n";
            return si;
        }
        custard::Header head;
        read(si, head);
        if (!si) {
            // std::cerr << "custard::read head broke stream\n";
            return si;
        }
        filename = head.name();
        filesize = head.size();
        return si;
    }



    // Stream header to ostream, leaving stream ready to write member
    // body data.  The caller MUST write additional number of bytes as
    // indicated by head.size() in order to avoid a corrupted stream.
    inline
    std::ostream& write(std::ostream& so, const Header& head)
    {
        so << "name " << head.name() << "\n"
           << "mode " << head.mode() << "\n"
           << "mtime " << head.mtime() << "\n"
           << "uid " << head.uid() << "\n"
           << "gid " << head.gid() << "\n"
           << "uname " << head.uname() << "\n"
           << "gname " << head.gname() << "\n"
           << "body " << head.size() << "\n";

        return so;
    }

    // The minimum header write:
    inline
    std::ostream& write(std::ostream& so,
                        const std::string& filename, size_t filesize)
    {
        so << "name " << filename << "\n"
           << "body " << filesize << "\n";
        return so;
    }


}

#endif
