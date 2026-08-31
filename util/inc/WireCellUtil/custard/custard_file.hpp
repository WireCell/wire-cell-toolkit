/**
   custard file

   Provide support for directly reading/writing tar files.

 */

#ifndef CUSTARD_FILE_HPP
#define CUSTARD_FILE_HPP

#include "custard.hpp"

namespace custard {
    class File {
        size_t loc{0}; // current read/write file location
        Header head;
      public:

        const Header& header() const {
            return head;
        }
        Header& header() {
            return head;
        }

        void clear() {
            loc = 0;
            head.clear();
        }

        // Start a file of given file name and target file size.
        void write_start(std::ostream& so, std::string filename,
                         size_t datasize)
        {
            loc = 0;
            head.init(filename, datasize);
            so.write((char*)head.as_bytes(), 512);
        }

        // Write some data to current file, not overflowing file size.
        // Return how data consumed.
        size_t write_data(std::ostream& so, const char* data,
                        size_t datasize)
        {
            size_t remain = head.size() - loc;
            size_t take = std::min(datasize, remain);
            so.write(data, take);
            loc += take;
            if (loc == head.size()) {
                write_pad(so);
                clear();
            }

            return take;
        }

        // Start a file and write entire or intial file data.  Return
        // how much of data consumed.
        size_t write_file(std::ostream& so, std::string filename,
                        const char* data, size_t datasize)
        {
            write_start(so, filename, datasize);
            size_t ret = write_data(so, data, datasize);
            return ret;
        }

        // Write cpad to next 512 boundary
        void write_pad(std::ostream& so)
        {
            // already aligned -> no padding (512-loc%512 would wrongly give 512)
            size_t pad = loc%512 ? 512 - loc%512 : 0;
            if (pad == 0) {
                return;
            }
            for (size_t ind=0; ind<pad; ++ind) {
                so.put('\0');
            };
        }
        
        // GNU manual on tar file format: at the end of the archive
        // file there are two 512-byte blocks filled with binary zeros
        // as an end-of-file marker.
        //
        // GNU tar, actual: lol, nah, I'mma gonna write 20! lololulz
        //
        // With only 2, GNU tar reads the archive fine and will
        // happily write 2 with "tar -b2".  With -b0 tar complains
        // about invalid blocking factor but with -b1 it quietly still
        // uses 2.  Or, may be it wants to have an even size mod 1024.
        void write_finish(std::ostream& so, size_t nblocks=2)
        {
            nblocks = std::max(2UL, nblocks);
            for (size_t ind=0; ind<nblocks*512; ++ind) {
                so.put('\0');
            };
        }

        // Begin a file read, return file size or zero on error.
        size_t read_start(std::istream& si)
        {
            clear();
            si.read(head.as_bytes(), 512);
            if (! si) {
                clear();
                return 0;
            }
            return head.size();
        }

        size_t read_data(std::istream& si, char* buf, size_t bufsize)
        {
            size_t remain = head.size() - loc;
            size_t take = std::min(bufsize, remain);
            si.read(buf, take);
            if (!si) {
                take = si.gcount();
            }
            loc += take;
            if (loc == head.size()) {
                read_pad(si);
                clear();
            }
            return take;
        }

        size_t read_file(std::istream& si, char* buf, size_t bufsize)
        {
            size_t got = read_start(si);
            if (got == 0) {
                return 0;
            }
            return read_data(si, buf, bufsize);
        }
        
        void read_pad(std::istream& si) {
            // already aligned -> no padding (512-loc%512 would wrongly give 512)
            const size_t jump = loc%512 ? 512 - loc%512 : 0;
            if (jump) {
                si.seekg(jump, si.cur);
            }
        }

    };

}
#endif
