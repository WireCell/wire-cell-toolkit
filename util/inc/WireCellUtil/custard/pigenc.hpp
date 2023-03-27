/**
   pigenc is an array serialization codec for Numpy format.
 */

#ifndef PIGENC_HPP
#define PIGENC_HPP

#include <iostream>             // testing
#include <istream>
#include <ostream>
#include <string>
#include <regex>
#include <complex>

#include "nlohmann/json.hpp"

namespace pigenc {


    using json = nlohmann::json;

    using shape_t = std::vector<size_t>;

    inline
    json parse_header(std::string header) {

        std::regex q("(')");
        std::regex op("(\\()");
        std::regex cp("(\\))");
        std::regex truify("(True)");
        std::regex falsify("(False)");
        std::regex trailcomma("(, ?)([}\\]])");

        header = std::regex_replace(header, q, "\"");
        header = std::regex_replace(header, op, "[");
        header = std::regex_replace(header, cp, "]");
        header = std::regex_replace(header, truify, "true");
        header = std::regex_replace(header, falsify, "false");
        header = std::regex_replace(header, trailcomma, "$2");

        return json::parse(header);
    }

    inline
    json read_header(std::istream& si) {
        std::string mver(8, 0);
        si.read(&mver[0], 8);
        if (!si) return json();
        uint32_t hsize{0};
        if (mver[6] == '\x01') {
            uint16_t s{0};
            si.read((char*)&s, 2);
            if (!si) return json();
            hsize = s;
        }
        else {
            si.read((char*)&hsize, 4);
            if (!si) return json();
        }

        std::string header(hsize, 0);
        si.read(&header[0], hsize);
        if (!si) return json();

        auto ret = parse_header(header);
        ret["version"] = {(int)mver[6], (int)mver[7]};
        return ret;
    }


    inline 
    std::string make_header(const std::string& desc, 
                            const shape_t& shape,
                            bool fortran_order=false)
    {
        std::string tf = fortran_order ? "True" : "False";
        std::stringstream ss;
        ss << "{'descr': '" << desc << "', "
           << "'fortran_order': " << tf << ", "
           << "'shape': (";
        for (auto s : shape) {
            ss << s << ", ";
        }
        ss << "), }";
        std::string dict = ss.str();

        int siz = 10 + dict.size();
        int pad = 64 - siz%64;
        uint16_t hsize = pad + siz - 10;

        std::string ret = "\x93NUMPY";
        ret.push_back( '\x01' );
        ret.push_back( '\x00' );
        ret.push_back( *(((char*)&hsize)+0) );
        ret.push_back( *(((char*)&hsize)+1) );
        ret += dict;
        for (int ind=0; ind<pad-1; ++ind) {
            ret.push_back(' ');
        }
        ret.push_back('\n');
        return ret;
    }

    template<typename T> std::string dtype()        {
        std::string tn = typeid(T).name();
        throw std::domain_error("pigenc: unsupported type: " + tn);
    }
    template<> inline std::string dtype<char>()     { return "c"; }
    template<> inline std::string dtype<int8_t>()   { return "<i1"; }
    template<> inline std::string dtype<uint8_t>()  { return "<u1"; }
    template<> inline std::string dtype<int16_t>()  { return "<i2"; }
    template<> inline std::string dtype<uint16_t>() { return "<u2"; }
    template<> inline std::string dtype<int32_t>()  { return "<i4"; }
    template<> inline std::string dtype<uint32_t>() { return "<u4"; }
    template<> inline std::string dtype<int64_t>()  { return "<i8"; }
    template<> inline std::string dtype<uint64_t>() { return "<u8"; }
    template<> inline std::string dtype<float>()    { return "<f4"; }
    template<> inline std::string dtype<std::complex<float>>()    { return "<c8"; }
    template<> inline std::string dtype<double>()   { return "<f8"; }

    // This assumes dtype string like above eg "...NN" where NN is number of bytes.
    inline size_t dtype_size(std::string dt)
    {
        if (dt.empty()) {
            throw std::domain_error("pigenc::dtype_size: empty dtype string");
        }
        if (dt.size() == 1 and dt[0] == 'c') {
            return 1;           // char has no number.
        }
        while (dt.size() and (dt[0]-'0' < 0 or dt[0]-'0' > 9)) {
            dt.erase(dt.begin());
        }
        if (dt.empty()) {
            throw std::domain_error("pigenc::dtype_size: unknown dtype");
        }
        return std::stol(dt);
    }

    /// Write bytes with all specified.
    inline
    std::ostream& write_data(std::ostream& os,
                             char* bytes,
                             const shape_t& shape,
                             const std::string& dtype_name,
                             size_t dtype_size,
                             bool fortran_order=false)
    {
        auto head = make_header(dtype_name, shape, fortran_order);
        size_t nbytes = dtype_size;
        for (auto s: shape) {
            nbytes *= s;
        }
        os.write(head.data(), head.size());
        os.write(bytes, nbytes);
        return os;
    }

    /// Write simple, raw array of elements type Type
    template<typename Type>
    std::ostream& write_data(std::ostream& os,
                             Type* data,
                             const shape_t& shape,
                             bool fortran_order=false)
    {
        return write_data(os, (char*)data, shape, dtype<Type>(),
                          sizeof(Type), fortran_order);

        // auto head = make_header(dtype<Type>(), shape, fortran_order);
        // size_t size = 1;
        // for (auto s: shape) {
        //     size *= s;
        // }
        // os.write(head.data(), head.size());
        // os.write((char*)data, size*sizeof(Type));
        // return os;
    }
        

    /// Represent the npy header for a simple array (not structured).
    class Header {
        // String rep of Python dictionary describing the array
        std::string descr_{""};

        // Array dimension sizes
        shape_t shape_{};

        // storage order
        bool fortran_order_{false};

        // derived values:

        std::string header_string_{""};

      public:

        // reset to empty/initial value
        void clear()
        {
            descr_ = "";
            shape_.clear();
            fortran_order_ = false;
            header_string_ = "";
        }

        // Set the header values given C++ element type
        template<typename Type>
        void set(const shape_t& shape,
                 bool fortran_order=false)
        {
            shape_ = shape;
            fortran_order_ = fortran_order;
            descr_ = pigenc::dtype<Type>();
            update();
        }

        // Set the header given descr as string
        void set(const shape_t& shape,
                 const std::string& descr,
                 bool fortran_order=false)
        {
            shape_ = shape;
            fortran_order_ = fortran_order;
            descr_ = descr;
            update();
        }

        size_t type_size() const {
            // note, only simple arrays have descr as dtype
            return dtype_size(descr_);
        }

        // Number of elements in the array
        size_t array_size() const {
            size_t size = 1;
            for (auto s: shape_) {
                size *= s;
            }
            return size;
        }
        
        // The size of the array part in bytes.  
        size_t data_size() const {
            return type_size() * array_size();
        }

        size_t header_size() const {
            return str().size();
        }

        // Size of the entire file
        size_t file_size() const {
            return header_size() + data_size();
        }

        bool fortran_order() const {
            return fortran_order_;
        }
        const shape_t& shape() const {
            return shape_;
        }
        shape_t shape() {
            return shape_;
        }

        // Update derived values
        void update() {
            std::string tf = fortran_order_ ? "True" : "False";
            std::stringstream ss;
            ss << "{'descr': '" << descr_ << "', "
               << "'fortran_order': " << tf << ", "
               << "'shape': (";
            if (shape_.size() == 1) {
                ss << shape_[0] << ",), }";
            }
            else {
                std::string comma = "";
                for (auto s : shape_) {
                    ss << comma << s;
                    comma = ", ";
                }
                ss << "), }";
            }
            std::string dict = ss.str();

            int siz = 10 + dict.size();
            int pad = 64 - siz%64;
            uint16_t hlen = pad + siz - 10;

            std::string ret = "\x93NUMPY";
            ret.push_back( '\x01' );
            ret.push_back( '\x00' );
            ret.push_back( *(((char*)&hlen)+0) );
            ret.push_back( *(((char*)&hlen)+1) );
            ret += dict;
            for (int ind=0; ind<pad-1; ++ind) {
                ret.push_back(' ');
            }
            ret.push_back('\n');
            header_string_ = ret;
            //assert(hlen+10 == str().size());
        }

        std::string str() const {
            return header_string_;
        }

        std::string dtype() const {
            return descr_;
        }

        // Read self from istream
        std::istream& read(std::istream& si)
        {
            auto dat = read_header(si);
            if (dat.is_null() or !si) {
                si.setstate(std::ios::failbit);
                return si;
            }
            descr_ = dat["descr"];
            fortran_order_ = dat["fortran_order"];
            shape_.clear();
            for (auto s : dat["shape"]) {
                size_t siz = s;
                shape_.push_back(siz);
            }
            update();
            return si;
        }

        std::ostream& write(std::ostream& so)
        {
            so.write(header_string_.data(), header_string_.size());
            return so;
        }

    };

    // Embody one pigenc/npy file.  It may be streamed in/out or have
    // data set/get in app.  Note that data size MUST be made
    // consistent with dtype and shape in the header.
    class File {
    public:
        using data_t = std::vector<char>;

        File() {}

        /// Create with consistent pairing of head and data.  
        File(const Header& head, const data_t& data) : head_(head), data_(data) {}

        template<typename Scalar>
        void set(const data_t& data, const shape_t& shape, bool fortran_order=false)
        {
            head_.set<Scalar>(shape, fortran_order);
            data_ = data;
        }

        void clear()
        {
            head_.clear();
            data_.clear();
        }

        Header& header() { return head_; }
        const Header& header() const { return head_; }

        data_t& data() { return data_; }
        data_t data() const { return data_; }

        template<typename Type>
        const Type* as_type() const
        {
            if (pigenc::dtype<Type>() == head_.dtype()) {
                return reinterpret_cast<const Type*>(data_.data());
            }
            return nullptr;
        }

        template<typename Type>
        std::vector<Type> as_vec() const
        {
            const Type* dat = as_type<Type>();
            if (!dat) {
                return std::vector<Type>();
            }
            std::vector<Type> ret(dat, dat + head_.array_size());
            return ret;
        }

        std::istream& read(std::istream& si)
        {
            head_.read(si);
            if (!si) { return si; }
            data_.resize(head_.data_size(), 0);
            si.read(&data_[0], data_.size());
            return si;
        }

        std::ostream& write(std::ostream& so)
        {
            head_.write(so);
            if (!so) { return so; }
            // std::cerr << "pigenc::File::write " << data_.size() << "\n";
            so.write(&data_[0], data_.size());
            return so;
        }
        
    private:
        Header head_;
        data_t data_;

    };

}

#endif

