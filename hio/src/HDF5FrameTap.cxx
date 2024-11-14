#include "WireCellHio/HDF5FrameTap.h"
#include "WireCellIface/ITrace.h"

#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"

#include <string>
#include <vector>


/// macro to register name - concrete pair in the NamedFactory
/// @param NAME - used to configure node in JSON/Jsonnet
/// @parame CONCRETE - C++ concrete type
/// @parame ... - interfaces
WIRECELL_FACTORY(HDF5FrameTap,
                 WireCell::Hio::HDF5FrameTap,
                 WireCell::INamed,
                 WireCell::IFrameFilter,
                 WireCell::ITerminal,
                 WireCell::IConfigurable)

using namespace WireCell;

Hio::HDF5FrameTap::HDF5FrameTap()
  : Aux::Logger("HDF5FrameTap","hio")
{
}

Hio::HDF5FrameTap::~HDF5FrameTap() {}


void Hio::HDF5FrameTap::configure(const WireCell::Configuration &cfg)
{
    std::vector<std::string> ignored = {"anode","high_throughput"};
    for (const auto& ign : ignored) {
        if (cfg[ign.c_str()].isNull()) {
            continue;
        }
        log->warn("The parameter \"{}\" is ignored, consider fixing your config", ign);
    }

    m_digitize = get<bool>(cfg, "digitize", m_digitize);
    m_baseline = get(cfg, "baseline", m_baseline);
    m_scale = get(cfg, "scale", m_scale);
    m_offset = get(cfg, "offset", m_offset);
    m_gzip = get(cfg, "gzip", m_gzip);

    if (m_gzip) {
        if (!H5Zfilter_avail(H5Z_FILTER_DEFLATE)) {
            raise<ValueError>("HDF5 gzip filter not available");
        }
        unsigned int filter_info;
        /*herr_t status = */H5Zget_filter_info (H5Z_FILTER_DEFLATE, &filter_info);
        if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
             !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED) ) {
            raise<ValueError>("HDF5 gzip filter not available for encoding and decoding");
        }

        auto jchunk = cfg["chunk"];
        if (jchunk.isInt()) {
            int chunk = jchunk.asInt();
            m_chunk[0] = chunk;
            m_chunk[1] = chunk;
        }
        else if (jchunk.isArray()) {
            m_chunk[0] = jchunk[0].asInt();
            m_chunk[1] = jchunk[1].asInt();
        }
    }

    log->debug("digitize={} baseline={} scale={} offset={} gzip={} chunking:[{},{}]",
               m_digitize, m_baseline, m_scale, m_offset, m_gzip, m_chunk[0], m_chunk[1]);

    std::string fn = cfg["filename"].asString();
    if (fn.empty()) {
        THROW(ValueError() << errmsg{"Must provide output filename to HDF5FrameTap"});
    }

    m_hfile = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (m_hfile == H5I_INVALID_HID) {
        raise<IOError>("Failed to create file %s", fn);
    }
    log->debug("created file {}", fn);

    m_cfg = cfg;


}
void Hio::HDF5FrameTap::finalize()
{
    H5Fclose(m_hfile);
}

WireCell::Configuration Hio::HDF5FrameTap::default_configuration() const
{
    Configuration cfg;

    // If digitize is true, then samples as 16 bit ints.  Otherwise
    // save as 32 bit floats.
    cfg["digitize"] = false;

    // This number is set to the waveform sample array before any
    // charge is added.
    cfg["baseline"] = m_baseline;

    // This number will be multiplied to each waveform sample before
    // casting to dtype.
    cfg["scale"] = m_scale;

    // This number will be added to each scaled waveform sample before
    // casting to dtype.
    cfg["offset"] = m_offset;

    // Set these only if you want to override what the frame actually provides.
    // cfg["nticks"] = 6000;
    // cfg["tick0"] = 0;

    // The frame tags to consider for saving.  If null or empty then all traces are used.
    cfg["trace_tags"] = Json::arrayValue;
    // The summary tags to consider for saving
    // cfg["summary_tags"] = Json::arrayValue;
    // The channel mask maps to consider for saving
    // cfg["chanmaskmaps"] = Json::arrayValue;

    // The output file name to write.  Only compressed (zipped) Numpy
    // files are supported.  Writing is always in "append" mode.  It's
    // up to the user to delete a previous instance of the file if
    // it's old contents are not wanted.
    cfg["filename"] = "wct-frame.hdf5";

    cfg["gzip"] = m_gzip;

    // ignored, ...for now?
    // cfg["chunk"] = Json::arrayValue;
    // cfg["high_throughput"] = true;

    return cfg;
}

bool Hio::HDF5FrameTap::operator()(const IFrame::pointer &inframe, IFrame::pointer &outframe)
{
    if (!inframe) {
        log->debug("EOS at call={}", m_calls++);
        outframe = nullptr;
        return true;
    }

    // Fixme: all this config machination should be in configure(), not here.

    // const int tick0 = m_cfg["tick0"].asInt();
    // const int nticks = m_cfg["nticks"].asInt();
    // const int tbeg = tick0;
    // const int tend = tick0 + nticks - 1;
    // auto channels = m_anode->channels();
    // const int cbeg = channels.front();
    // const int cend = channels.back();
    // log->debug("{}: t: {} - {}; c: {} - {}", m_cfg["anode"].asString(), tbeg, tend, cbeg, cend);

    outframe = inframe;  // pass through actual frame

    const std::string mode = "a";
    const std::string fname = m_cfg["filename"].asString();

    // Eigen3 array is indexed as (irow, icol) or (ichan, itick)
    // one row is one channel, one column is a tick.

    if (m_cfg["trace_tags"].isNull() or m_cfg["trace_tags"].empty()) {
        m_cfg["trace_tags"][0] = "";
    }

    // fixme: chunking?  gzipping?

    // fixme: replace this logging with call to taginfo().
    std::stringstream ss;
    ss << "HDF5FrameTap: see frame #" << inframe->ident() << " with " << inframe->traces()->size()
       << " traces with frame tags:";
    for (auto t : inframe->frame_tags()) {
        ss << " \"" << t << "\"";
    }
    ss << " and trace tags:";
    for (auto t : inframe->trace_tags()) {
        ss << " \"" << t << "\"";
    }
    ss << " looking for tags:";
    for (auto jt : m_cfg["trace_tags"]) {
        ss << " \"" << jt.asString() << "\"";
    }
    log->debug(ss.str());

    for (auto jtag : m_cfg["trace_tags"]) {
        const std::string tag = jtag.asString();
        auto traces = Aux::tagged_traces(inframe, tag);
        if (traces.empty()) {
            log->warn("no traces for tag: \"{}\"", tag);
            continue;
        }

        auto channels = Aux::channels(traces);
        size_t nrows = channels.size();

        auto tbinmm = Aux::tbin_range(traces);
        const int tbin = get(m_cfg, "tick0", tbinmm.first);
        size_t ncols = get(m_cfg, "nticks", tbinmm.second - tbin);


        Array::array_xxf arr = Array::array_xxf::Zero(nrows, ncols) + m_baseline;
        Aux::fill(arr, traces, channels.begin(), channels.end(), tbin);
        arr = arr * m_scale + m_offset;

        // log->debug("saving {}: ncols={} nrows={} with tbin={}, qtot={}",
        //            tag, ncols, nrows, tbin, arr.sum());

        int sequence = inframe->ident();
        {  // the 2D frame array, long hand....
            const std::string aname = String::format("/%d/frame_%s", sequence, tag.c_str());

            // link create property list
            hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
            if (lcpl == H5I_INVALID_HID) {
                raise<IOError>("failed to create link creation property list");
            }
            if (H5Pset_char_encoding(lcpl, H5T_CSET_UTF8) < 0) {
                raise<IOError>("failed to set UTF encoding");
            }
            if (H5Pset_create_intermediate_group(lcpl, 1) < 0) {
                raise<IOError>("failed to set creation of intermediate groups");
            }

            // data set create property list
            hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
            if (dcpl == H5I_INVALID_HID) {
                raise<IOError>("failed to create dataset creation property list");
            }
            if (m_gzip) {
                if (H5Pset_deflate (dcpl, m_gzip) < 0) {
                    raise<IOError>("failed set gzip compression of {}", m_gzip);
                }
                if (H5Pset_chunk (dcpl, 2, m_chunk.data()) < 0) {
                    raise<IOError>("failed set chunk of [{},{}]", m_chunk[0], m_chunk[1]);
                }
            }

            int ndims = 2;
            // BIG WARNING: HDF5 is sanely C-order (row-major, column iteration
            // is fastest) but eigen defaults to Fortran-order (column-major,
            // row iteration is fastest).  We thus do a little transpose dance
            // below, just before writing.
            hsize_t dims[2] = {nrows, ncols};
            hid_t dspace = H5Screate_simple(ndims, dims, NULL);
            if (dspace == H5I_INVALID_HID) {
                raise<IOError>("failed to create data space of %d x %d", dims[0], dims[1]);
            }

            hid_t dtype = 0;
            if (m_digitize) {
                dtype = H5Tcopy(H5T_NATIVE_SHORT);
            }
            else {
                dtype = H5Tcopy(H5T_NATIVE_FLOAT);
            }
            if (dtype == H5I_INVALID_HID) {
                raise<IOError>("failed to make dtype for {}", m_digitize ? "short int" : "float");
            }

            // The dataset
            hid_t dset = H5Dcreate2(m_hfile, aname.c_str(), dtype, dspace, lcpl, dcpl, H5P_DEFAULT);
            if (dset == H5I_INVALID_HID) {
                raise<IOError>("failed to create dataset %s", aname);
            }
            
            herr_t status;
            if (m_digitize) {
                using ROWI = Eigen::Array<int16_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
                ROWI sarr = arr.cast<short>();
                status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, sarr.data());
            }
            else {
                using ROWF = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
                ROWF farr = arr;
                status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, farr.data());
            }
            if (status == H5I_INVALID_HID) {
                raise<IOError>("HDF5FrameTap: failed to write samples");
            }

            H5Dclose(dset);
            H5Tclose(dtype);
            H5Sclose(dspace);
            H5Pclose(dcpl);
            H5Pclose(lcpl);

            log->debug("saved {} with {} channels {} ticks @t={} ms qtot={}", aname, nrows, ncols,
                     inframe->time() / units::ms, arr.sum());
        }

        {  // the channel array
            const std::string aname = String::format("/%d/channels_%s", sequence, tag.c_str());

            int ndims = 1;
            hsize_t dims[1] = {nrows};
            hid_t dspace = H5Screate_simple(ndims, dims, NULL);
            hid_t dtype = H5Tcopy(H5T_NATIVE_INT);
            hid_t dset = H5Dcreate1(m_hfile, aname.c_str(), dtype, dspace, H5P_DEFAULT);
            herr_t status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, channels.data());
            if (status == H5I_INVALID_HID) {
                raise<IOError>("HDF5FrameTap: failed to write channels");
            }
            H5Dclose(dset);
            H5Tclose(dtype);
            H5Sclose(dspace);
        }

        {  // the tick array
            const std::string aname = String::format("/%d/tickinfo_%s", sequence, tag.c_str());
            const std::vector<double> tickinfo{inframe->time(), inframe->tick(), (double) tbin};

            // fixme: this small array would perhaps be better as attributes.
            int ndims = 1;
            hsize_t dims[1] = {3};
            hid_t dspace = H5Screate_simple(ndims, dims, NULL);
            hid_t dtype = H5Tcopy(H5T_NATIVE_DOUBLE);
            hid_t dset = H5Dcreate1(m_hfile, aname.c_str(), dtype, dspace, H5P_DEFAULT);
            herr_t status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tickinfo.data());
            if (status == H5I_INVALID_HID) {
                raise<IOError>("HDF5FrameTap: failed to write tickinfo");
            }
            H5Dclose(dset);
            H5Tclose(dtype);
            H5Sclose(dspace);
        }
    }

    ++m_calls;
    return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
