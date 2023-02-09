#include <iostream>
#include <stdio.h>
#include <limits>
#include <float.h>

#include "debug-helper.cpp"

#include "WireCellCuda/cuFftDFT.h"
#include "debug-helper.cpp"

using namespace WireCell;

using complex_t = IDFT::complex_t;

// Read comment for `is_equal()` below for details
float EQUAL_EPSILON = 1000*FLT_EPSILON;
float ARR_EQUAL_THRESH = 0.005;

// Tests whether a and b are about equal ( difference is a result of rounding errors )
bool is_equal(float a, float b)
{ return fabs(a-b) <= max(fabs(a), fabs(b))*EQUAL_EPSILON; }

bool is_equal(complex_t a, complex_t b)
{ return is_equal(std::real(a), std::real(b)) && is_equal(std::imag(a), std::imag(b)); }

// Returns the number of corresponding values of `A` and `B` 
// that are not equal ( according to `is_equal` )
bool is_arr_equal(complex_t *A, complex_t *B, int n)
{
    int m = 0;
    for( int i = 0; i < n; ++i )
        if( !is_equal(A[i], B[i]) )
            ++m;
            
    return (float)m/n < ARR_EQUAL_THRESH;
}

// Returns the number of corresponding values of `A` that  
// are not equal to `1+0i` ( according to `is_equal` )
bool is_arr_all_1(complex_t *A, int n)
{
    int m = 0;
    for( int i = 0; i < n; ++i )
        if( !is_equal(A[i], complex_t(1, 0)) )
            ++m;

    return (float)m/n < ARR_EQUAL_THRESH;
}


// Generates an array of random floats/complex numbers between 0 and 1
void alloc_rand_arr( complex_t **arr, int n )
{
    *arr = (complex_t *) malloc( sizeof(complex_t) * n );
    for( int i = 0; i < n; ++i ) 
        (*arr)[i] = complex_t(
            (float)rand()/float(RAND_MAX),
            (float)rand()/float(RAND_MAX)
        );
}

// Divides every element of `arr` by `m`
void div_arr(complex_t *arr, float m, int n)
{
    for( int i = 0; i < n; ++i )
        arr[i] /= m;
}

void test_1d_unit()
{
    complex_t *in, *out;
    Cuda::cuFftDFT fftObj;
    int n = 24; // Size of the array to test with

    in = (complex_t *)calloc(n, sizeof(complex_t));
    out = (complex_t *)malloc(n*sizeof(complex_t));
    in[0] = 1;

    fftObj.fwd1d(in, out, n);

    if(!is_arr_all_1(out, n))
        throw runtime_error("1d transform unit test failed");

    free(in);
    free(out);
}

void test_1d_rand()
{
    complex_t *in, *out;
    Cuda::cuFftDFT fftObj;
    int n = 24; // Size of the array to test with

    // Test 1: [1, 0, 0, ...]
    alloc_rand_arr(&in, n);
    out = (complex_t *)malloc(n*sizeof(complex_t));

    fftObj.fwd1d(in, out, n);
    fftObj.inv1d(out, out, n);

    div_arr(out, n, n); // Normalize transform output

    if(!is_arr_equal(in, out, n))
        throw runtime_error("1d transform rand test failed");

    free(in);
    free(out);
}

void test_1b_unit()
{
    Cuda::cuFftDFT fftObj;
    complex_t *in, *out;

    // Length of each signal and batch size
    int n = 8;
    int m = 3; 

    in = (complex_t *)calloc(n*m, sizeof(complex_t));
    out = (complex_t *)malloc( sizeof(complex_t)*n*m );

    // Axis 0
    for( int i = 0; i < m; ++i )
        in[i*n] = 1;

    fftObj.fwd1b(in, out, m, n, 0);

    if(!is_arr_all_1(out, n*m))
        throw runtime_error("1b transform axis 0 unit test failed");

    for( int i = 0; i < n*m; ++i )        
        in[i] = 0;

    // Axis 1
    for( int i = 0; i < n; ++i )
        in[i] = 1;

    fftObj.fwd1b(in, out, n, m, 1);

    if(!is_arr_all_1(out, n*m))
        throw runtime_error("1b transform axis 1 unit test failed");
    
    free(in);
    free(out);
}

void test_1b_rand()
{
    Cuda::cuFftDFT fftObj;
    complex_t *in, *out;

    // Length of each signal and batch size
    int n = 8;
    int m = 3; 

    alloc_rand_arr(&in, n*m);
    out = (complex_t *)malloc( sizeof(complex_t)*n*m );

    fftObj.fwd1b(in, out, m, n, 0);
    fftObj.inv1b(out, out, m, n, 0);
    div_arr(out, n, n*m);

    if(!is_arr_equal(in, out, n*m))
        throw runtime_error("1b transform axis 0 rand test failed");

    // Axis 1
    fftObj.fwd1b(in, out, n, m, 1);
    fftObj.inv1b(out, out, n, m, 1);
    div_arr(out, m, n*m);

    if(!is_arr_equal(in, out, n*m))
        throw runtime_error("1b transform axis 1 rand test failed");
    
    free(in);
    free(out);
}

void test_2d_unit()
{
    Cuda::cuFftDFT fftObj;
    complex_t *in, *out;

    // Length of each signal and batch size
    int n = 8;
    int m = 3; 

    in = (complex_t *)calloc(n*m, sizeof(complex_t));
    out = (complex_t *)malloc(sizeof(complex_t)*n*m);

    in[0] = 1;

    fftObj.fwd2d(in, out, n, m);

    if(!is_arr_all_1(out, n*m))
        throw runtime_error("2d transform unit test failed");
    
    free(in);
    free(out);

}

void test_2d_rand()
{
    Cuda::cuFftDFT fftObj;
    complex_t *in, *out;

    // Length of each signal and batch size
    int n = 8;
    int m = 3; 

    alloc_rand_arr(&in, n*m);
    out = (complex_t *)malloc(sizeof(complex_t)*n*m);

    fftObj.fwd2d(in, out, n, m);
    fftObj.inv2d(out, out, n, m);
    div_arr(out, n*m, n*m);

    if(!is_arr_equal(in, out, n*m))
        throw runtime_error("2d transform rand test failed");
    
    free(in);
    free(out);
}

int main()
{
    srand (time(NULL));

    // Test 1d transform
    test_1d_unit();
    test_1d_rand();

    // Test 1b transform
    test_1b_unit();
    test_1b_rand();

    // Test 2d transform
    test_2d_unit();
    test_2d_rand();

    return 0;
}
