#include "WireCellUtil/Testing.h"
#include "WireCellUtil/ExecMon.h"

#include "WireCellUtil/Eigen.h"
#include <iostream>

using namespace std;
using namespace Eigen;
using WireCell::ExecMon;

void test_matrix_by_array(ExecMon& em)
{
    em("mba: start");
    Matrix3f mat1, mat2, mat3, mat4, mat5;
    mat1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    mat2 << 1, 1, 1, 1, 1, 1, 1, 1, 1;
    mat3 = mat1.array() + mat2.array();
    mat4 = mat1 * mat2;
    // mat5 = mat1 * mat2.array(); compilation error, can't mix!

    cerr << "mat1:\n" << mat1 << endl;
    cerr << "mat2:\n" << mat2 << endl;
    cerr << "mat3:\n" << mat3 << endl;
    cerr << "mat4:\n" << mat4 << endl;
}

int main()
{
    ExecMon em;
    test_matrix_by_array(em);

    cerr << em.summary() << endl;

    return 0;
}
