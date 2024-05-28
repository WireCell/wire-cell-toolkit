#include "WireCellUtil/Persist.h"
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/Response.h"
#include "WireCellUtil/Logging.h"
#include <cmath>
#include <set>

using spdlog::error;

using namespace WireCell;

Response::Schema::FieldResponse::~FieldResponse() {}

Response::Schema::PlaneResponse::~PlaneResponse() {}

Response::Schema::PathResponse::~PathResponse() {}

/**
 frdict['FieldResponse'].keys()
 ['origin', 'axis', 'period', 'tstart', 'planes']

 frdict['FieldResponse']['planes'][0]['PlaneResponse'].keys()
 ['paths', 'planeid', 'pitch']

 frdict['FieldResponse']['planes'][0]['PlaneResponse']['paths'][0]['PathResponse'].keys()
 ['current', 'wirepos', 'pitchpos']

 frdict['FieldResponse']['planes'][0]['PlaneResponse']['paths'][0]['PathResponse']['current']['array'].keys()
 ['shape', 'elements']

 */
WireCell::Response::Schema::FieldResponse WireCell::Response::Schema::load(const char* filename)
{
    if (!filename) {
        error("Response::Schema::load(): empty field response file name");
        return FieldResponse();
    }
    Json::Value top = WireCell::Persist::load(filename);
    if (top.isNull()) {
        error("Response::Schema::load(): failed to load {}", filename);
        return FieldResponse();
    }
    Json::Value fr = top["FieldResponse"];

    using namespace WireCell::Response::Schema;

    std::vector<PlaneResponse> planes;
    for (auto plane : fr["planes"]) {
        auto plr = plane["PlaneResponse"];

        std::vector<PathResponse> paths;
        for (auto path : plr["paths"]) {
            auto par = path["PathResponse"];
            WireCell::Waveform::realseq_t current;
            for (auto c : par["current"]["array"]["elements"]) {
                current.push_back(c.asDouble());
            }
            auto wcpath = PathResponse(current, par["pitchpos"].asDouble(), par["wirepos"].asDouble());
            paths.push_back(wcpath);
        }

        PlaneResponse wcplr(paths, plr["planeid"].asInt(), plr["location"].asDouble(), plr["pitch"].asDouble());

        planes.push_back(wcplr);
    }

    auto adir = fr["axis"];
    auto axis = WireCell::Vector(adir[0].asDouble(), adir[1].asDouble(), adir[2].asDouble());
    auto ret = FieldResponse(planes, axis, fr["origin"].asDouble(), fr["tstart"].asDouble(), fr["period"].asDouble(),
                             fr["speed"].asDouble());
    return ret;
}

void Response::Schema::dump(const char* filename, const Response::Schema::FieldResponse& fr) {}

/// Warning!  this function is NOT GENERAL.  It is actually specific
/// to Garfield 1D line of paths with half the impact positions
/// represented!
Response::Schema::FieldResponse Response::wire_region_average(const Response::Schema::FieldResponse& fr)
{
    using namespace WireCell::Waveform;
    using namespace WireCell::Response::Schema;

    std::vector<PlaneResponse> newplanes;
    for (auto plane : fr.planes) {
        std::vector<PathResponse> newpaths;

        double pitch = plane.pitch;

        std::map<int, realseq_t> avgs;
        // std::map<int, int> nums;

        std::map<int, realseq_t> fresp_map;
        std::map<int, std::pair<double, double>> pitch_pos_range_map;

        // figure out the range of each response ...

        int nsamples = 0;

        for (auto path : plane.paths) {
            int eff_num = path.pitchpos / (0.01 * pitch);
            if (fresp_map.find(eff_num) == fresp_map.end()) {
                fresp_map[eff_num] = path.current;
            }
            else {
                nsamples = path.current.size();
                for (int k = 0; k != nsamples; k++) {
                    fresp_map[eff_num].at(k) = (fresp_map[eff_num].at(k) + path.current.at(k)) / 2.;
                }
            }
            if (fresp_map.find(-eff_num) == fresp_map.end()) {
                fresp_map[-eff_num] = path.current;
            }
            else {
                int nsamples = path.current.size();
                for (int k = 0; k != nsamples; k++) {
                    fresp_map[-eff_num].at(k) = (fresp_map[-eff_num].at(k) + path.current.at(k)) / 2.;
                }
            }
        }

        std::vector<int> pitch_pos;
        for (auto it = fresp_map.begin(); it != fresp_map.end(); it++) {
            pitch_pos.push_back((*it).first);
        }

        double min = -1e9;
        double max = 1e9;
        for (size_t i = 0; i != pitch_pos.size(); i++) {
            int eff_num = pitch_pos.at(i);
            if (i == 0) {
                pitch_pos_range_map[eff_num] =
                    std::make_pair(min, (pitch_pos.at(i) + pitch_pos.at(i + 1)) / 2. * 0.01 * pitch);
            }
            else if (i == pitch_pos.size() - 1) {
                pitch_pos_range_map[eff_num] =
                    std::make_pair((pitch_pos.at(i) + pitch_pos.at(i - 1)) / 2. * 0.01 * pitch, max);
            }
            else {
                pitch_pos_range_map[eff_num] =
                    std::make_pair((pitch_pos.at(i) + pitch_pos.at(i - 1)) / 2. * 0.01 * pitch,
                                   (pitch_pos.at(i) + pitch_pos.at(i + 1)) / 2. * 0.01 * pitch);
            }
        }

        // figure out how many wires ...
        std::set<int> wire_regions;
        for (size_t i = 0; i != pitch_pos.size(); i++) {
            if (pitch_pos.at(i) > 0) {
                wire_regions.insert(round((pitch_pos.at(i) * 0.01 * pitch - 0.001 * pitch) / pitch));
            }
            else {
                wire_regions.insert(round((pitch_pos.at(i) * 0.01 * pitch + 0.001 * pitch) / pitch));
            }
        }

        // do the average ...
        for (auto it = wire_regions.begin(); it != wire_regions.end(); it++) {
            int wire_no = *it;
            if (avgs.find(wire_no) == avgs.end()) {
                avgs[wire_no] = realseq_t(nsamples);
            }
            for (auto it1 = fresp_map.begin(); it1 != fresp_map.end(); it1++) {
                int resp_num = (*it1).first;
                realseq_t& response = (*it1).second;
                double low_limit = pitch_pos_range_map[resp_num].first;
                double high_limit = pitch_pos_range_map[resp_num].second;
                if (low_limit < (wire_no - 0.5) * pitch) {
                    low_limit = (wire_no - 0.5) * pitch;
                }
                if (high_limit > (wire_no + 0.5) * pitch) {
                    high_limit = (wire_no + 0.5) * pitch;
                }

                //

                if (high_limit > low_limit) {
                    for (int k = 0; k != nsamples; k++) {
                        avgs[wire_no].at(k) += response.at(k) * (high_limit - low_limit) / pitch;
                    }
                }
            }
        }

        // do average.
        for (auto it : avgs) {
            int region = it.first;
            realseq_t& response = it.second;

            // pack up everything for return.
            newpaths.push_back(PathResponse(response, region * pitch, 0.0));
        }
        newplanes.push_back(PlaneResponse(newpaths, plane.planeid, plane.location, plane.pitch));
    }
    return FieldResponse(newplanes, fr.axis, fr.origin, fr.tstart, fr.period, fr.speed);
}

Response::Schema::FieldResponse Response::average_1D(const Response::Schema::FieldResponse& fr)
{
    using namespace WireCell::Waveform;
    using namespace WireCell::Response::Schema;

    FieldResponse fr_wire_avg = Response::wire_region_average(fr);

    std::vector<PlaneResponse> newplanes;
    for (auto plane : fr_wire_avg.planes) {
        std::vector<PathResponse> newpaths;

        int nsamples = Response::as_array(plane).cols();

        realseq_t ave_response(nsamples, 0);

        for (auto path : plane.paths) {
            for (int k = 0; k != nsamples; k++) {
                ave_response.at(k) += path.current.at(k);
            }
        }

        newpaths.push_back(PathResponse(ave_response, 0.0, 0.0));

        newplanes.push_back(PlaneResponse(newpaths, plane.planeid, plane.location, plane.pitch));
    }
    return FieldResponse(newplanes, fr.axis, fr.origin, fr.tstart, fr.period, fr.speed);
}

Array::array_xxf Response::as_array(const Schema::PlaneResponse& pr, int set_nrows, int set_ncols)
{
    int nrows = pr.paths.size();
    int ncols = pr.paths[0].current.size();
    Array::array_xxf ret = Array::array_xxf::Zero(set_nrows, set_ncols);  // warning, uninitialized

    if (set_nrows < nrows || set_ncols < ncols) {
        error("Response: array dimension not correct! ");
        return ret;
    }

    for (int irow = 0; irow < nrows; ++irow) {
        if (irow < set_nrows) {
            auto& path = pr.paths[irow];
            for (int icol = 0; icol < ncols; ++icol) {
                if (icol < set_ncols)
                    ret(irow, icol) = path.current[icol];  // maybe there is a fast way to do this copy?
            }
        }
    }
    return ret;
}

Array::array_xxf Response::as_array(const Schema::PlaneResponse& pr)
{
    int nrows = pr.paths.size();
    int ncols = pr.paths[0].current.size();
    Array::array_xxf ret(nrows, ncols);  // warning, uninitialized

    for (int irow = 0; irow < nrows; ++irow) {
        auto& path = pr.paths[irow];
        for (int icol = 0; icol < ncols; ++icol) {
            ret(irow, icol) = path.current[icol];  // maybe there is a fast way to do this copy?
        }
    }
    return ret;
}

Response::Generator::~Generator() {}

// FIXME: eradicate Domain in favor of Binning
WireCell::Waveform::realseq_t Response::Generator::generate(const WireCell::Waveform::Domain& domain, int nsamples)
{
    WireCell::Waveform::realseq_t ret(nsamples);
    const double tick = (domain.second - domain.first) / nsamples;
    for (int ind = 0; ind < nsamples; ++ind) {
        double t = domain.first + ind * tick;
        ret[ind] = (*this)(t);
    }
    return ret;
}
WireCell::Waveform::realseq_t Response::Generator::generate(const WireCell::Binning& tbins)
{
    const int nsamples = tbins.nbins();
    WireCell::Waveform::realseq_t ret(nsamples, 0.0);
    for (int ind = 0; ind < nsamples; ++ind) {
        const double time = tbins.center(ind);
        ret[ind] = (*this)(time);
    }
    return ret;
}

/*
  Cold Electronics response function.

  How was this function derived?

    1. Hucheng provided a transfer function of our electronics, which
       is the Laplace transformation of the shaping function.

    2. We did a anti-Laplace inverse of the shaping function

    3. The one you saw in the code is basically the result of the inversion

  - time is time in system of units

  - gain in units of volts/charge gives peak value of the response to
    a delta function of current integrating to unit charge.

  - shaping is the shaping time in system of units

  - the hard-coded numbers are the result of the inverting the
    Lapalace transformation in Mathematica.

 */
double Response::coldelec(double time, double gain, double shaping)
{
    if (time <= 0 || time >= 10 * units::microsecond) {  // range of validity
        return 0.0;
    }

    const double reltime = time / shaping;

    // a scaling is needed to make the anti-Lapalace peak match the expected gain
    // fixme: this scaling is slightly dependent on shaping time.  See response.py
    gain *= 10 * 1.012;

    return 4.31054 * exp(-2.94809 * reltime) * gain - 2.6202 * exp(-2.82833 * reltime) * cos(1.19361 * reltime) * gain -
           2.6202 * exp(-2.82833 * reltime) * cos(1.19361 * reltime) * cos(2.38722 * reltime) * gain +
           0.464924 * exp(-2.40318 * reltime) * cos(2.5928 * reltime) * gain +
           0.464924 * exp(-2.40318 * reltime) * cos(2.5928 * reltime) * cos(5.18561 * reltime) * gain +
           0.762456 * exp(-2.82833 * reltime) * sin(1.19361 * reltime) * gain -
           0.762456 * exp(-2.82833 * reltime) * cos(2.38722 * reltime) * sin(1.19361 * reltime) * gain +
           0.762456 * exp(-2.82833 * reltime) * cos(1.19361 * reltime) * sin(2.38722 * reltime) * gain -
           2.620200 * exp(-2.82833 * reltime) * sin(1.19361 * reltime) * sin(2.38722 * reltime) * gain -
           0.327684 * exp(-2.40318 * reltime) * sin(2.5928 * reltime) * gain +
           +0.327684 * exp(-2.40318 * reltime) * cos(5.18561 * reltime) * sin(2.5928 * reltime) * gain -
           0.327684 * exp(-2.40318 * reltime) * cos(2.5928 * reltime) * sin(5.18561 * reltime) * gain +
           0.464924 * exp(-2.40318 * reltime) * sin(2.5928 * reltime) * sin(5.18561 * reltime) * gain;
}

/*
 * A parameterized cold electronics response function (adapted from Xin Qian's formula)
 * Parameters (k3,k4,k5,k6) control undershoot/overshoot behaviors
 */
double Response::paramscoldelec(double time, double gain, double shaping, double k3, double k4, double k5, double k6)
{
    if (time <= 0 || time >= 50 * units::microsecond) {  // range of validity
        return 0.0;                                      // do we want to ignore long tail?
    }

    // const double reltime = time / shaping;
    double t = time / units::microsecond;
    double A0 = gain *10.0;
    double tp = shaping / units::microsecond;
    double CT = 1./1.996;
    double A = A0 * 2.7433/pow(tp*CT,4);
    double p0 = 1.477/tp/CT;
    double pr1 = 1.417/tp/CT;
    double pr2 = 1.204/tp/CT;
    double pi1 = 0.598/tp/CT;
    double pi2 = 1.299/tp/CT;
  
    return A*((-(k3*k4) + pow(k4,2) + k3*k5 - k4*k5)/(exp(k4*t)*(k4 - k6)*(k4 - p0)*(pow(k4,2) + pow(pi1,2) - 2*k4*pr1 + pow(pr1,2))*(pow(k4,2) + pow(pi2,2) - 2*k4*pr2 + pow(pr2,2))) +
             (-(k3*k5) + k3*k6 + k5*k6 - pow(k6,2))/(exp(k6*t)*(k4 - k6)*(k6 - p0)*(pow(k6,2) + pow(pi1,2) - 2*k6*pr1 + pow(pr1,2))*(pow(k6,2) + pow(pi2,2) - 2*k6*pr2 + pow(pr2,2))) +
             (-(k3*k5) + k3*p0 + k5*p0 - pow(p0,2))/(exp(p0*t)*(k4 - p0)*(-k6 + p0)*(pow(p0,2) + pow(pi1,2) - 2*p0*pr1 + pow(pr1,2))*(pow(p0,2) + pow(pi2,2) - 2*p0*pr2 + pow(pr2,2))) +
     (pi1*((pow(pi1,2) + pow(pr1,2))*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)) +
              k5*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) +
                 k6*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)))) + k4*((pow(pi1,2) + pow(pr1,2))*(2*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2))) +
              k5*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) - k6*p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2))) +
              k6*(-pow(pi1,4) + pow(pi1,2)*(pow(pi2,2) - 2*pow(pr1,2) + 2*p0*pr2 + pow(pr2,2)) + pr1*(pr1*(pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) - 2*p0*(pow(pi2,2) - pr1*pr2 + pow(pr2,2))))) +
           k3*(-((pow(pi1,2) + pow(pr1,2))*(4*pow(pi1,2)*pr1 - 2*pow(pi2,2)*pr1 - 4*pow(pr1,3) - 2*pow(pi1,2)*pr2 + 6*pow(pr1,2)*pr2 - 2*pr1*pow(pr2,2) + p0*(-pow(pi1,2) + pow(pi2,2) + 3*pow(pr1,2) - 4*pr1*pr2 + pow(pr2,2)) +
                   k6*(-pow(pi1,2) + pow(pi2,2) - (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)))) + k5*(-pow(pi1,4) + pow(pi1,2)*(pow(pi2,2) - 4*p0*pr1 + 10*pow(pr1,2) + 2*p0*pr2 - 8*pr1*pr2 + pow(pr2,2)) +
                 k6*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2))) +
                 pr1*(2*p0*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) - pr1*(3*pow(pi2,2) + 5*pow(pr1,2) - 8*pr1*pr2 + 3*pow(pr2,2)))) +
              k4*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)) +
                 k5*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) +
                    k6*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2))))))*cos(pi1*t) -
        ((pow(pi1,2) + pow(pr1,2))*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
              k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
              k5*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)))) +
           k4*((pow(pi1,2) + pow(pr1,2))*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
              k5*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
                 k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2)))) +
              k6*((pow(pi1,2) + pow(pr1,2))*(pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pr1 - 2*pr2)) -
                 p0*(pow(pi1,4) + pow(pr1,2)*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - 2*pow(pr1,2) + 2*pr1*pr2 + pow(pr2,2))))) +
           k3*((pow(pi1,2) + pow(pr1,2))*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2))) +
              k5*(5*pow(pi1,4)*pr1 - 3*pow(pi1,2)*pow(pi2,2)*pr1 - 10*pow(pi1,2)*pow(pr1,3) + pow(pi2,2)*pow(pr1,3) + pow(pr1,5) - 2*pow(pi1,4)*pr2 + 12*pow(pi1,2)*pow(pr1,2)*pr2 - 2*pow(pr1,4)*pr2 - 3*pow(pi1,2)*pr1*pow(pr2,2) +
                 pow(pr1,3)*pow(pr2,2) - p0*(pow(pi1,4) + pow(pr1,2)*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) + 6*pow(pr1,2) - 6*pr1*pr2 + pow(pr2,2))) +
                 k6*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)))) +
              k4*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
                 k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
                 k5*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                    k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2))))))*sin(pi1*t))/
             (exp(pr1*t)*pi1*(pow(k4,2) + pow(pi1,2) - 2*k4*pr1 + pow(pr1,2))*(pow(k6,2) + pow(pi1,2) - 2*k6*pr1 + pow(pr1,2))*(pow(p0,2) + pow(pi1,2) - 2*p0*pr1 + pow(pr1,2))*
        (pow(pi1,4) - 2*pow(pi1,2)*(pow(pi2,2) - pow(pr1 - pr2,2)) + pow(pow(pi2,2) + pow(pr1 - pr2,2),2))) +
     (-(pi2*(k4*(-((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - pow(pr2,2)) - 2*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)))) +
                k5*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2))) +
                k6*(pow(pi2,4) - pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - 2*pow(pr2,2)) - pow(pi1,2)*(pow(pi2,2) + pr2*(-2*p0 + pr2)) - (pr1 - pr2)*pr2*(-2*p0*pr1 + pr2*(pr1 + pr2)))) +
             (pow(pi2,2) + pow(pr2,2))*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) +
                k5*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2))))) +
             k3*((pow(pi2,2) + pow(pr2,2))*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))) +
                k5*(pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) + 4*p0*pow(pi2,2)*pr2 + 8*pow(pi2,2)*pr1*pr2 - 2*p0*pow(pr1,2)*pr2 - 10*pow(pi2,2)*pow(pr2,2) + 6*p0*pr1*pow(pr2,2) + 3*pow(pr1,2)*pow(pr2,2) -
                   4*p0*pow(pr2,3) - 8*pr1*pow(pr2,3) + 5*pow(pr2,4) - pow(pi1,2)*(pow(pi2,2) + (2*p0 - 3*pr2)*pr2) + k6*p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*k6*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))) +
                k4*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) +
                   k5*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                      2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))))))*cos(pi2*t)) +
        ((pow(pi2,2) + pow(pr2,2))*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
              k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
              k5*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 + pow(pi1,2)*pow(pr2,2) -
                 6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)))) +
           k4*(-(k6*(pow(pi2,2) + pow(pr2,2))*((pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + pr2))) +
              k6*p0*(pow(pi2,4) + pow(pr1 - pr2,2)*pow(pr2,2) - pow(pi2,2)*(pow(pr1,2) + 2*pr1*pr2 - 2*pow(pr2,2)) + pow(pi1,2)*(-pow(pi2,2) + pow(pr2,2))) +
              (pow(pi2,2) + pow(pr2,2))*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
              k5*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
                 k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))))) +
           k3*((pow(pi2,2) + pow(pr2,2))*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 +
                 pow(pi1,2)*pow(pr2,2) - 6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2))) -
              k5*(p0*pow(pi1,2)*pow(pi2,2) - p0*pow(pi2,4) - 2*pow(pi2,4)*pr1 + p0*pow(pi2,2)*pow(pr1,2) - 3*pow(pi1,2)*pow(pi2,2)*pr2 + 5*pow(pi2,4)*pr2 - 6*p0*pow(pi2,2)*pr1*pr2 - 3*pow(pi2,2)*pow(pr1,2)*pr2 -
                 p0*pow(pi1,2)*pow(pr2,2) + 6*p0*pow(pi2,2)*pow(pr2,2) + 12*pow(pi2,2)*pr1*pow(pr2,2) - p0*pow(pr1,2)*pow(pr2,2) + pow(pi1,2)*pow(pr2,3) - 10*pow(pi2,2)*pow(pr2,3) + 2*p0*pr1*pow(pr2,3) + pow(pr1,2)*pow(pr2,3) -
                 p0*pow(pr2,4) - 2*pr1*pow(pr2,4) + pow(pr2,5) + k6*(-pow(pi2,4) + (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi1,2)*(pow(pi2,2) + (p0 - pr2)*pr2) + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - 3*p0*pr2 - 6*pr1*pr2 + 6*pow(pr2,2)))) +
              k4*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
                 k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
                 k5*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 + pow(pi1,2)*pow(pr2,2) -
                    6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                    k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2))))))*sin(pi2*t))/
             (exp(pr2*t)*pi2*(pow(pi1,4) - 2*pow(pi1,2)*(pow(pi2,2) - pow(pr1 - pr2,2)) + pow(pow(pi2,2) + pow(pr1 - pr2,2),2))*(pow(k4,2) + pow(pi2,2) - 2*k4*pr2 + pow(pr2,2))*(pow(k6,2) + pow(pi2,2) - 2*k6*pr2 + pow(pr2,2))*
       (pow(p0,2) + pow(pi2,2) - 2*p0*pr2 + pow(pr2,2))));

}


double Response::hf_filter(double freq, double sigma, double power, bool flag)
{
    if (flag) {
        if (freq == 0) return 0;
    }

    return exp(-0.5 * pow(freq / sigma, power));
}

double Response::lf_filter(double freq, double tau) { return 1 - exp(-pow(freq / tau, 2)); }

Response::ColdElec::ColdElec(double gain, double shaping)
  : _g(gain)
  , _s(shaping)
{
}
Response::ColdElec::~ColdElec() {}

double Response::ColdElec::operator()(double time) const { return coldelec(time, _g, _s); }

Response::ParamsColdElec::ParamsColdElec(double gain, double shaping, double k3, double k4, double k5, double k6)
  : _g(gain)
  , _s(shaping)
  , _k3(k3)
  , _k4(k4)
  , _k5(k5)
  , _k6(k6)
{
}
Response::ParamsColdElec::~ParamsColdElec() {}

double Response::ParamsColdElec::operator()(double time) const { return paramscoldelec(time, _g, _s, _k3, _k4, _k5, _k6); }


/*
  Warm Electronics response function.

  How was this function derived?

  - from the article: arXiv:1805.03931 [physics.ins-det]
  - copy from ElectronicsResponseBesselApprox_tool.cc in icaruscode v08_50_00
  - SPICE simulation (ADC dynamic range 0-3.3V)
  qLSB (i.e. 1ADCxtick): 65e ~ 31mVxus/fC ~ 0.026fC/(ADCxus)
  Qpeak: 181ADC / 100fC ~ 9.517mV/fC
  - CERN test
  qLSB: 69.2e ~ 29.11mVxus/fC ~ 0.0277fC/(ADCxus)
  Qpeak: (4732/7.7)ADC/52.5fC = 9.433mV/fC

  (11.16.2020 W.Gu) The two results are roughly consistent, the measured result
  is recommended in the icaruscode (FCperADCMicroS: 2.7e-2), while we recommend
  the setting as follows.
        gain : 19.2855 mV/fC
        shaping : 1.2 microsecond
  As this setting provides a better agreement for both Qpeak and qLSB. The
  channel-to-channel variation can be calibrated in the future.
  CAVEAT: Here, the "gain" does not mean the absolute amplitude as in the
  ColdElectronics Reponse.
  
 */
double Response::warmelec(double time, double gain, double shaping)
{
    if (time <= 0 || time >= 10 * units::microsecond) {  // range of validity
        return 0.0;
    }

    const double reltime = time / shaping;

    double leftArg = (0.9 * reltime) * (0.9 * reltime);
    double rightArg = (0.5 * reltime) * (0.5 * reltime);

    return gain * ((1. - exp(-0.5 * leftArg)) * exp(-0.5 * rightArg));
}

Response::WarmElec::WarmElec(double gain, double shaping)
  : _g(gain)
  , _s(shaping)
{
}
Response::WarmElec::~WarmElec() {}

double Response::WarmElec::operator()(double time) const { return warmelec(time, _g, _s); }

Response::SimpleRC::SimpleRC(double width, double tick, double offset)
  : _width(width)
  , _offset(offset)
  , _tick(tick)

{
}
Response::SimpleRC::~SimpleRC() {}
double Response::SimpleRC::operator()(double time) const
{
    double ret = 0;
    if (_width > 0) {
        ret += -_tick / _width *
               exp(-(time - _offset) / _width);  // _tick here is to make this RC response integrated in each bin
    }
    if (time < _offset + _tick) {  // just the first bin
        ret += 1.0;                // delta function
    }
    return ret;
}

// Vary field response to study systematics
// Currently a Gaussian function
Response::SysResp::SysResp(double tick, double magnitude, double smear, double offset)
  : _tick(tick)
  , _mag(magnitude)
  , _smear(smear)
  , _offset(offset)
{
}
Response::SysResp::~SysResp() {}
double Response::SysResp::operator()(double time) const
{
    double ret = 0;
    if (_smear > 0) {
        ret = _tick * exp(-0.5 * pow((time - _offset) / _smear, 2)) / _smear * 0.3989422804;
    }
    else if (time < _tick + _offset && time >= _offset) {
        ret += 1.0;
    }
    else {
        ret = 0;
    }
    return ret * _mag;
}

Response::LfFilter::LfFilter(double tau)
  : _tau(tau)
{
}

Response::LfFilter::~LfFilter() {}

double Response::LfFilter::operator()(double freq) const { return lf_filter(freq, _tau); }

Response::HfFilter::HfFilter(double sigma, double power, bool flag)
  : _sigma(sigma)
  , _power(power)
  , _flag(flag)
{
}

Response::HfFilter::~HfFilter() {}

double Response::HfFilter::operator()(double freq) const { return hf_filter(freq, _sigma, _power, _flag); }
