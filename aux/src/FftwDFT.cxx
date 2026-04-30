#include "WireCellAux/FftwDFT.h"
#include "WireCellUtil/NamedFactory.h"

#include <fftw3.h>
#include <mutex>
#include <shared_mutex>
#include <unordered_map>

WIRECELL_FACTORY(FftwDFT, WireCell::Aux::FftwDFT, WireCell::IDFT)


using namespace WireCell;

using plan_key_t = int64_t;
using plan_type = fftwf_plan;
using plan_map_t = std::unordered_map<plan_key_t, plan_type>;
using plan_val_t = fftwf_complex;

// Make a key by which a plan is known.  dir should be FFTW_FORWARD or
// FFTW_BACKWARD and "axis" is -1 for all or in {0,1} for one of 2D.
// For 1D, use the default axis=-1.
//
// Imp note: The key is slightly over-specified as we keep one
// independent cache for each of the six methods.  The "dir" is
// thus redundant.
static
plan_key_t make_key(const void * src, void * dst, int nrows, int ncols, int dir, int axis=-1)
{
    ++axis;                     // need three positive values, default is both axis
    bool inverse = dir == FFTW_BACKWARD;
    bool inplace = (dst==src);
    // Drop 'aligned' from the key: plans are created with FFTW_UNALIGNED so a
    // single plan serves any buffer alignment. Removing this address-derived bit
    // from the key is what makes the cache reproducible across runs (heap
    // addresses vary between processes under ASLR and even without it).
    int64_t key = ( ( (((int64_t)nrows) << 32)| (ncols<<5 ) | (axis<<3) | (inverse<<2) | (inplace<<1) ) << 1 ) + 1;
    return key;
}

// Look up a plan by key or return NULL
static
plan_type get_plan(std::shared_mutex& mutex, plan_map_t& plans, plan_key_t key)
{
    std::shared_lock lock(mutex);
    auto it = plans.find(key);
    if (it == plans.end()) {
        return NULL;
    }
    return it->second;
}


// #include <iostream>             // debugging

using planner_function = std::function<plan_type()>;

// This wraps plan lookup, possible plan creation and subsequent plan
// execution so that we get thread-safe plan caching.
template<typename ValueType>
void doit(std::shared_mutex& mutex, plan_map_t& plans, plan_key_t key,
          ValueType* src, ValueType* dst,
          planner_function make_plan,
          std::function<void(const plan_type, ValueType*, ValueType*)> exec_plan)
{
    auto plan = get_plan(mutex, plans, key);
    if (!plan) {
        std::unique_lock lock(mutex);
        // Check again in case another thread snakes us.
        auto it = plans.find(key);
        if (it == plans.end()) {
            //std::cerr << "make plan for " << key << std::endl;
            plan = make_plan();
            plans[key] = plan;
        }
        else {
            plan = it->second;
        }
    }
    //fftwf_execute_dft(plan, src, dst);
    exec_plan(plan, src, dst);
}


static
plan_val_t* pval_cast( const IDFT::complex_t * p)
{ 
    return const_cast<plan_val_t*>( reinterpret_cast<const plan_val_t*>(p) ); 
}


void Aux::FftwDFT::fwd1d(const complex_t* in, complex_t* out, int ncols) const
{
    static std::shared_mutex mutex;
    static plan_map_t plans;
    static const int dir = FFTW_FORWARD;
    auto src = pval_cast(in);
    auto dst = pval_cast(out);
    auto key = make_key(src, dst, 1, ncols, dir);
    doit<plan_val_t>(mutex, plans, key, src, dst, [&]( ) {
        return fftwf_plan_dft_1d(ncols, src, dst, dir, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
    }, fftwf_execute_dft);
}
void Aux::FftwDFT::inv1d(const complex_t* in, complex_t* out, int ncols) const
{
    static std::shared_mutex mutex;
    static plan_map_t plans;
    static const int dir = FFTW_BACKWARD;
    auto src = pval_cast(in);
    auto dst = pval_cast(out);
    auto key = make_key(src, dst, 1, ncols, dir);

    doit<plan_val_t>(mutex, plans, key, src, dst, [&]( ) {
        return fftwf_plan_dft_1d(ncols, src, dst, dir, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
    }, fftwf_execute_dft);

    // Apply 1/n normalization
    for (int ind=0; ind<ncols; ++ind) {
        out[ind] /= ncols;
    }
}


fftwf_plan plan_1b(fftwf_complex *in, fftwf_complex *out,
                   int nrows, int ncols, int sign, int axis)
{
    // (r,c) element at in + r*stride + c*dist

    const int rank = 1;         // dimension of transform
    int n = ncols;              // along rows
    int howmany = nrows;
    int stride = 1;
    int dist = ncols;
    if (axis == 0) {            // along columns
        n = nrows;
        howmany = ncols;
        stride = ncols;
        dist = 1;
    }
    int *inembed=&n, *onembed=&n;

    unsigned int flags =  FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED;

    return fftwf_plan_many_dft(rank, &n, howmany,
                               in, inembed,
                               stride, dist,
                               out, onembed,
                               stride, dist,
                               sign, flags);
}


void Aux::FftwDFT::fwd1b(const complex_t* in, complex_t* out, int nrows, int ncols, int axis) const
{
    static std::shared_mutex mutex;
    static plan_map_t plans;
    static const int dir = FFTW_FORWARD;
    auto src = pval_cast(in);
    auto dst = pval_cast(out);
    auto key = make_key(src, dst, nrows, ncols, dir, axis);

    doit<plan_val_t>(mutex, plans, key, src, dst, [&]( ) {
        return plan_1b(src, dst, nrows, ncols, dir, axis);
    }, fftwf_execute_dft);
}


void Aux::FftwDFT::inv1b(const complex_t* in, complex_t* out, int nrows, int ncols, int axis) const
{
    static std::shared_mutex mutex;
    static plan_map_t plans;
    static const int dir = FFTW_BACKWARD;
    auto src = pval_cast(in);
    auto dst = pval_cast(out);
    auto key = make_key(src, dst, nrows, ncols, dir, axis);

    doit<plan_val_t>(mutex, plans, key, src, dst, [&]( ) {
        return plan_1b(src, dst, nrows, ncols, dir, axis);
    }, fftwf_execute_dft);

    // 1/n normalization
    const int norm = axis ? ncols : nrows;
    const int ntot = ncols*nrows;
    for (int ind=0; ind<ntot; ++ind) {
        out[ind] /= norm;
    }
}


void Aux::FftwDFT::fwd2d(const complex_t* in, complex_t* out, int nrows, int ncols) const
{
    static std::shared_mutex mutex;
    static plan_map_t plans;
    static const int dir = FFTW_FORWARD;
    auto src = pval_cast(in);
    auto dst = pval_cast(out);
    auto key = make_key(src, dst, nrows, ncols, dir);
    doit<plan_val_t>(mutex, plans, key, src, dst, [&]( ) {
        return fftwf_plan_dft_2d(ncols, nrows, src, dst, dir, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
    }, fftwf_execute_dft);
}


void Aux::FftwDFT::inv2d(const complex_t* in, complex_t* out, int nrows, int ncols) const
{
    static std::shared_mutex mutex;
    static plan_map_t plans;
    static const int dir = FFTW_BACKWARD;
    auto src = pval_cast(in);
    auto dst = pval_cast(out);
    auto key = make_key(src, dst, nrows, ncols, dir);
    doit<plan_val_t>(mutex, plans, key, src, dst, [&]( ) {
        return fftwf_plan_dft_2d(ncols, nrows, src, dst, dir, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
    }, fftwf_execute_dft);

    // reverse normalization
    const int ntot = ncols*nrows;
    for (int ind=0; ind<ntot; ++ind) {
        out[ind] /= ntot;
    }
}


// based on example from fftw3 faq
static
plan_type transpose_plan_complex(plan_val_t *in, plan_val_t *out, int rows, int cols)
{
    const unsigned flags = FFTW_ESTIMATE|FFTW_UNALIGNED; /* other flags are possible */
    fftw_iodim howmany_dims[2];

    howmany_dims[0].n  = rows;
    howmany_dims[0].is = cols;
    howmany_dims[0].os = 1;

    howmany_dims[1].n  = cols;
    howmany_dims[1].is = 1;
    howmany_dims[1].os = rows;

    return fftwf_plan_guru_dft(/*rank=*/ 0, /*dims=*/ NULL,
                               /*howmany_rank=*/ 2, howmany_dims,
                               in, out, /*sign=*/ 0, flags);
}
void Aux::FftwDFT::transpose(const complex_t* in, complex_t* out,
                             int nrows, int ncols) const
{
    static std::shared_mutex mutex;
    static plan_map_t plans;
    static const int dir = 0;
    auto src = pval_cast(in);
    auto dst = pval_cast(out);
    auto key = make_key(src, dst, nrows, ncols, dir);
    doit<plan_val_t>(mutex, plans, key, src, dst, [&]( ) {
        return transpose_plan_complex(src, dst, nrows, ncols);
    }, fftwf_execute_dft);
}

static
plan_type transpose_plan_real(float *in, float *out, int rows, int cols)
{
    const unsigned flags = FFTW_ESTIMATE|FFTW_UNALIGNED; /* other flags are possible */
    fftw_iodim howmany_dims[2];

    howmany_dims[0].n  = rows;
    howmany_dims[0].is = cols;
    howmany_dims[0].os = 1;

    howmany_dims[1].n  = cols;
    howmany_dims[1].is = 1;
    howmany_dims[1].os = rows;

    return fftwf_plan_guru_r2r(/*rank=*/ 0, /*dims=*/ NULL,
                               /*howmany_rank=*/ 2, howmany_dims,
                               in, out, /*kind=*/ NULL, flags);
}
void Aux::FftwDFT::transpose(const scalar_t* in, scalar_t* out,
                             int nrows, int ncols) const
{
    static std::shared_mutex mutex;
    static plan_map_t plans;
    static const int dir = 0;
    auto src = const_cast<scalar_t*>(in);
    auto dst = out;
    auto key = make_key(src, dst, nrows, ncols, dir);
    doit<float>(mutex, plans, key, src, dst, [&]( ) {
        return transpose_plan_real(src, dst, nrows, ncols);
    }, fftwf_execute_r2r);
}

Aux::FftwDFT::FftwDFT()
{
}
Aux::FftwDFT::~FftwDFT()
{
}

