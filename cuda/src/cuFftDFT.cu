#include "WireCellCuda/cuFftDFT.h"

#include <cuda.h>
#include <cufft.h>
#include <assert.h>
#include <iostream>

using namespace WireCell;
using complex_t = IDFT::complex_t;

/* Because otherwise a more explicit casting will be needed which
 * will reduce performance. */
static_assert(sizeof(cufftComplex) == sizeof(complex_t));

Cuda::memArgs<cufftComplex>
Cuda::cuFftDFT::preTransformTasks(
    Cuda::memArgs<complex_t> hmem, int size
) const
{
    // Allocate memory for the input and output arrays in the GPU   
    static Cuda::memArgs<cufftComplex> dmem;

    assert(cudaMalloc( &dmem.in, size*sizeof(cufftComplex) ) == cudaSuccess);

    if( dmem.in == dmem.out ) dmem.out = const_cast<cufftComplex*>(dmem.in);
    else assert(cudaMalloc( &dmem.out, size*sizeof(cufftComplex) ) == cudaSuccess);

    // Copy input array in CPU to the GPU
    assert(cudaMemcpy(
        (void *)dmem.in, (void *)hmem.in, 
        size*sizeof(cufftComplex), 
        cudaMemcpyHostToDevice
    ) == cudaSuccess);

    return dmem;
}

void Cuda::cuFftDFT::postTransformTasks(
    memArgs<complex_t> hmem, memArgs<cufftComplex> dmem, int size
) const
{
    // Copy result back to the CPU
    assert(cudaMemcpy(
        (void *)hmem.out, (void *)dmem.out, 
        size*sizeof(complex_t), 
        cudaMemcpyDeviceToHost
    ) == cudaSuccess);

    // Free memory on GPU
    assert(cudaFree((void *)dmem.in) == cudaSuccess);
    if( dmem.in != dmem.out )
        assert(cudaFree((void *)dmem.out) == cudaSuccess);
}

void Cuda::cuFftDFT::gen1d(const complex_t* in, complex_t* out, int size, int dir) const
{
    static cufftHandle plan;
    Cuda::memArgs<complex_t> hmem{in, out};
    Cuda::memArgs<cufftComplex> dmem = preTransformTasks(hmem, size);

    // Create the plan and perform the transform
    assert(cufftPlan1d(&plan, size, CUFFT_C2C, 1) == CUFFT_SUCCESS);
    assert(cufftExecC2C(plan, const_cast<cufftComplex*>(dmem.in), dmem.out, dir) == CUFFT_SUCCESS);
    
    postTransformTasks(hmem, dmem, size);
}

void Cuda::cuFftDFT::gen1b(const complex_t* in, complex_t* out, int nrows, int ncols, int axis, int dir) const
{
    static cufftHandle plan;
    Cuda::memArgs<complex_t> hmem{in, out};
    Cuda::memArgs<cufftComplex> dmem = preTransformTasks(hmem, nrows*ncols);

    // Create the plan and perform the transform
    assert(axis == 0 || axis == 1);
    if( axis == 0 )
    {
        assert(cufftPlanMany(
            &plan, 1, new int[2] {ncols, nrows}, 
            &nrows, 1, ncols, 
            &nrows, 1, ncols, 
            CUFFT_C2C, nrows 
        ) == CUFFT_SUCCESS);
    } else
    {
        assert(cufftPlanMany(
            &plan, 1, new int[2] {ncols, nrows}, 
            &ncols, nrows, 1, 
            &ncols, nrows, 1, 
            CUFFT_C2C, nrows 
        ) == CUFFT_SUCCESS);
    }

    assert(cufftExecC2C(plan, const_cast<cufftComplex*>(dmem.in), dmem.out, dir) == CUFFT_SUCCESS);
    
    postTransformTasks(hmem, dmem, nrows*ncols);
}

void Cuda::cuFftDFT::gen2d(const complex_t* in, complex_t* out, int nrows, int ncols, int dir) const
{
    static cufftHandle plan;
    Cuda::memArgs<complex_t> hmem{in, out};
    Cuda::memArgs<cufftComplex> dmem = preTransformTasks(hmem, nrows*ncols);

    assert(cufftPlan2d(&plan, nrows, ncols, CUFFT_C2C) == CUFFT_SUCCESS);
    assert(cufftExecC2C(plan, const_cast<cufftComplex*>(dmem.in), dmem.out, dir) == CUFFT_SUCCESS);
    
    postTransformTasks(hmem, dmem, nrows*ncols);
}

void Cuda::cuFftDFT::fwd1d(const complex_t* in, complex_t* out, int size) const
{ gen1d(in, out, size, CUFFT_FORWARD); }

void Cuda::cuFftDFT::inv1d(const complex_t* in, complex_t* out, int size) const 
{ gen1d(in, out, size, CUFFT_INVERSE); }

void Cuda::cuFftDFT::fwd1b(const complex_t* in, complex_t* out, int nrows, int ncols, int axis) const 
{ gen1b(in, out, nrows, ncols, axis, CUFFT_FORWARD); }

void Cuda::cuFftDFT::inv1b(const complex_t* in, complex_t* out, int nrows, int ncols, int axis) const 
{ gen1b(in, out, nrows, ncols, axis, CUFFT_INVERSE); }

void Cuda::cuFftDFT::fwd2d(const complex_t* in, complex_t* out, int nrows, int ncols) const 
{ gen2d(in, out, nrows, ncols, CUFFT_FORWARD); }

void Cuda::cuFftDFT::inv2d(const complex_t* in, complex_t* out, int nrows, int ncols) const 
{ gen2d(in, out, nrows, ncols, CUFFT_INVERSE); }

