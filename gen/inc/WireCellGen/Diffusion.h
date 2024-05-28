#ifndef WIRECELL_DIFFUSION
#define WIRECELL_DIFFUSION

#include "WireCellIface/IDiffusion.h"
#include "WireCellIface/IDepo.h"
#include "WireCellUtil/MultiArray.h"

namespace WireCell {

    class Diffusion : public IDiffusion {
        IDepo::pointer m_depo;
        boost::multi_array<double, 2> array;
        double lmin, tmin, lmax, tmax, lbin, tbin;

       public:
        Diffusion(IDepo::pointer depo, int nlong, int ntrans, double lmin, double tmin, double lmax, double tmax);

        Diffusion(const Diffusion& other);
        Diffusion& operator=(const Diffusion& other);

        virtual ~Diffusion();

        virtual IDepo::pointer depo() const;

        // Size of diffusion patch in both directions.
        virtual int lsize() const;
        virtual int tsize() const;

        virtual double get(int lind, int tind) const;
        virtual double set(int lind, int tind, double value);

        // Longitudinal position at index with extra offset 0.5 is bin center.
        virtual double lpos(int ind, double offset = 0.0) const;
        // Transverse position at index with extra offset 0.5 is bin center.
        virtual double tpos(int ind, double offset = 0.0) const;
    };
}  // namespace WireCell
#endif
