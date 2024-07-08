#include <stdio.h>
#include <iostream>
#include <vector>
#include <chrono>

class Vec_vec {
    public:
    Vec_vec (const double& x, const double& y, const double& z) : vec(3,0) {
        vec.resize(3, 0);
        vec[0] = x;
        vec[1] = y;
        vec[2] = z;
    }
    std::vector<double> vec{3};

    friend std::ostream& operator<<(std::ostream& os, const Vec_vec& v) {
        os << " capacity " << v.vec.capacity();
        os << " elems ";
        for (const auto& elem : v.vec) {
            os << elem << " ";
        }
        return os;
    }
};

class Vec_arr {
    public:
    Vec_arr (const double& x, const double& y, const double& z) {
        vec[0] = x;
        vec[1] = y;
        vec[2] = z;
    }
    double vec[3];

    friend std::ostream& operator<<(std::ostream& os, const Vec_arr& v) {
        os << " elems ";
        for (size_t i = 0; i < 3; ++i) {
            os << v.vec[i] << " ";
        }
        return os;
    }
};

int main()
{
    std::vector<double> v(3, 128);
    std::cout << v.at(0) << " " << v.at(1) << " " << v.at(2) << std::endl;

    // vector ctor/dtor performance
    auto start = std::chrono::high_resolution_clock::now();
    // Vec_vec vv(42,128,256);
    for (int i = 0; i < 1E7; ++i) {
        Vec_vec vv(42,128,256);
        // vv.vec.resize(3, 0);
        // vv.vec[0] = 42;
        // vv.vec[1] = 42;
        // vv.vec[2] = 42;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Vec_vec ctor/dtor " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    // vector assign performance
    start = std::chrono::high_resolution_clock::now();
    Vec_vec vv(0,0,0);
    for (int i = 0; i < 1E7; ++i) {
        vv.vec.resize(3, 0);
        vv.vec[0] = 42;
        vv.vec[1] = 128;
        vv.vec[2] = 256;
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Vec_vec assign " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    // array ctor/dtor performance
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1E7; ++i) {
        Vec_arr va(42,128,256);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Vec_arr ctor/dtor " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    

    return 0;
}
