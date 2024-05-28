#ifndef WIRECELLSIGPROC_RESPONSE
#define WIRECELLSIGPROC_RESPONSE

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Binning.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Array.h"

#include "WireCellUtil/Point.h"

namespace WireCell {

    namespace Response {

        //// Units notice: all quantities are expressed in the WCT
        //// system of unis.  In particular, time is not seconds.

        // These objects correspond to those defined in the Wire Cell
        // field response transfer file format schema.
        namespace Schema {

            // FIXME: this schema is very specific to Garfield 2D
            // results.  The namespace should reflect that and a more
            // generic interface should hide it.

            /// Hold information about the induced current response
            /// due to passage of a charge along one drift path.
            struct PathResponse {
                /// An array holding the induced current for the path on the wire-of-interest.
                WireCell::Waveform::realseq_t current;

                /// The position in the pitch direction to the starting point of the path.
                double pitchpos;

                /// The position along the wire direction to the starting point of the path.
                double wirepos;

                PathResponse()
                  : pitchpos(0.0)
                  , wirepos(-99999.0)
                {
                }
                PathResponse(const WireCell::Waveform::realseq_t& c, double p, double w)
                  : current(c)
                  , pitchpos(p)
                  , wirepos(w)
                {
                }

                ~PathResponse();
            };

            /// Hold information about the collection of induced
            /// current responses on one wire plane.
            struct PlaneResponse {
                /// List of PathResponse objects.
                std::vector<PathResponse> paths;

                /// A numerical identifier for the plane.
                int planeid;

                /// location, in direction of drift, of this plane (in
                /// same coordinate system as used by
                /// FieldResponse::origin).
                double location;

                /// The pitch distance between neighboring wires.
                double pitch;

                PlaneResponse()
                  : planeid(-1)
                  , location(0.0)
                  , pitch(0.0)
                {
                }
                PlaneResponse(const std::vector<PathResponse>& paths, int pid, double l, double p)
                  : paths(paths)
                  , planeid(pid)
                  , location(l)
                  , pitch(p)
                {
                }

                ~PlaneResponse();
            };

            /// Hold info about multiple plane responses in the detector.
            struct FieldResponse {
                /// List of PlaneResponse objects.
                std::vector<PlaneResponse> planes;

                /// A normalized 3-vector giving direction of axis
                /// (anti)parallel to nominal drift direction.
                WireCell::Vector axis;

                /// The location on the X-axis where drift paths
                /// begin.  See PlaneResponse::location.
                double origin;

                /// Time at which drift paths begin.
                double tstart;

                /// The sampling period of the response function.
                double period;

                /// The nominal drift speed.
                double speed;

                PlaneResponse* plane(int ident)
                {
                    for (auto& pr : planes) {
                        if (pr.planeid == ident) {
                            return &pr;
                        }
                    }
                    return nullptr;
                }
                const PlaneResponse* plane(int ident) const
                {
                    for (auto& pr : planes) {
                        if (pr.planeid == ident) {
                            return &pr;
                        }
                    }
                    return nullptr;
                }

                FieldResponse()
                  : origin(-999.0)
                  , tstart(-999.0)
                  , period(0.0)
                  , speed(0.0)
                {
                }
                FieldResponse(const std::vector<PlaneResponse>& planes, const WireCell::Vector& adir, double o,
                              double t, double p, double s)
                  : planes(planes)
                  , axis(adir)
                  , origin(o)
                  , tstart(t)
                  , period(p)
                  , speed(s)
                {
                }
                ~FieldResponse();
            };

            FieldResponse load(const char* filename);
            void dump(const char* filename, const FieldResponse& fr);

        }  // namespace Schema

        /// Return a reduced FieldResponse structure where the
        /// Path::Response::current arrays are reduced by averaging
        /// over each wire region.
        Schema::FieldResponse wire_region_average(const Schema::FieldResponse& fr);

        Schema::FieldResponse average_1D(const Schema::FieldResponse& fr);

        /// Return the plane's response as a 2D array.  This is a
        /// straight copy of the plane's current vectors into rows of
        /// the returned array.  The first "path" will be in row 0.
        /// Each column thus contains the same sample (aka tick)
        /// across all current waveforms.  Column 0 is first sample.
        /// The plane response input data is taken at face value an
        /// not attempt to resolve any implicit symmetries is made.
        Array::array_xxf as_array(const Schema::PlaneResponse& pr);
        Array::array_xxf as_array(const Schema::PlaneResponse& pr, int set_nrows, int set_ncols);

        /// The cold electronics response function.
        double coldelec(double time, double gain = 7.8, double shaping = 1.0 * units::us);

	/// The parameterized cold electronics response function.
	double paramscoldelec(double time, double gain = 7.8, double shaping = 1.0 * units::us, double k3=0.1, double k4=0.1, double k5=0.0, double k6=0.0);

        /// The warm electronics response function.
        double warmelec(double time, double gain = 30, double shaping = 1.3 * units::us);

        // HF filter format
        double hf_filter(double freq, double sigma = 1, double power = 2, bool zero_freq_removal = true);

        // LF filter format
        double lf_filter(double freq, double tau = 0.02);

        class Generator {
           public:
            virtual ~Generator();
            virtual double operator()(double time) const = 0;

            /// FIXME: eradicate Domain in favor of Binning.
            WireCell::Waveform::realseq_t generate(const WireCell::Waveform::Domain& domain, int nsamples);
            /// Lay down the function into a binned waveform.
            WireCell::Waveform::realseq_t generate(const WireCell::Binning& tbins);
        };

        /// A functional object caching gain and shape.
        class ColdElec : public Generator {
            const double _g, _s;

           public:
            // Create cold electronics response function.  Gain is an
            // arbitrary scale, typically in [voltage/charge], and
            // shaping time in WCT system of units.
            ColdElec(double gain = 14 * units::mV / units::fC, double shaping = 1.0 * units::us);
            virtual ~ColdElec();

            // Return the response at given time.  Time is in WCT
            // system of units.
            virtual double operator()(double time) const;
        };

        /// A functional object caching gain, shape and other params
        class ParamsColdElec : public Generator {
            const double _g, _s, _k3, _k4, _k5, _k6;

           public:
            // Create parameterized cold electronics response function.  Gain is an
            // arbitrary scale, typically in [voltage/charge], and
            // shaping time in WCT system of units.
            ParamsColdElec(double gain = 14 * units::mV / units::fC, double shaping = 2.0 * units::us, double k3=0.1, double k4=0.1, double k5=0.0, double k6=0.0);
            virtual ~ParamsColdElec();

            // Return the response at given time.  Time is in WCT
            // system of units.
            virtual double operator()(double time) const;
        };


        /// A functional object caching gain and shape.
        /// ICARUS warm electronics
        class WarmElec : public Generator {
            const double _g, _s;

           public:
            // Create warm electronics response function.  Gain is an
            // arbitrary scale, typically in [voltage/charge], and
            // shaping time in WCT system of units.
            WarmElec(double gain = 30 * units::mV / units::fC, double shaping = 1.3 * units::us);
            virtual ~WarmElec();

            // Return the response at given time.  Time is in WCT
            // system of units.
            virtual double operator()(double time) const;
        };

        /// A functional object giving the response as a function of
        /// time to a simple RC circuit.
        class SimpleRC : public Generator {
            const double _width, _offset, _tick;

           public:
            // Create (current) response function for a simple RC
            // circuit where a unit of charge is placed on the cap at
            // time offset and circuit has RC time constant of given
            // width.  Time is in WCT system of units.
            SimpleRC(double width = 1.0 * units::ms, double tick = 0.5 * units::us, double offset = 0.0);
            virtual ~SimpleRC();

            // Return the response at a given time.  Time in WCT
            // system of units.  Warning: to get the delta function,
            // one must call *exactly* at the offset time.
            virtual double operator()(double time) const;
        };

        class SysResp : public Generator {
            const double _tick, _mag, _smear, _offset;

           public:
            SysResp(double tick = 0.5 * units::us, double magnitude = 1.0, double smear = 0.0 * units::us,
                    double offset = 0.0 * units::us);
            virtual ~SysResp();
            virtual double operator()(double time) const;
        };

        class LfFilter : public Generator {
            const double _tau;

           public:
            LfFilter(double tau);
            virtual ~LfFilter();
            virtual double operator()(double freq) const;
        };

        class HfFilter : public Generator {
            const double _sigma, _power;
            const bool _flag;

           public:
            HfFilter(double sigma, double power, bool flag);
            virtual ~HfFilter();
            virtual double operator()(double freq) const;
        };

    }  // namespace Response
}  // namespace WireCell

#endif
