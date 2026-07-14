#include "WireCellMatch/PhotonLibraryModel.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/cnpy.h"

#include <algorithm>
#include <cmath>

using namespace WireCell;
using WireCell::Match::PhotonLibraryModel;

PhotonLibraryModel::PhotonLibraryModel(const std::string& meta_file)
{
    const std::string meta_path = Persist::resolve(meta_file);
    if (meta_path.empty()) {
        raise<ValueError>("PhotonLibraryModel: cannot resolve meta file '%s'", meta_file);
    }
    const auto meta = Persist::load(meta_path);
    if (!meta.isObject() || !meta["origin_cm"].isArray() || !meta["step_cm"].isArray()
        || !meta["n"].isArray() || !meta["vis_npy"].isString()) {
        raise<ValueError>("PhotonLibraryModel: malformed meta file '%s'", meta_path);
    }
    for (int i = 0; i < 3; ++i) {
        m_origin[i] = meta["origin_cm"][i].asDouble();
        m_step[i] = meta["step_cm"][i].asDouble();
        m_n[i] = meta["n"][i].asInt();
        if (m_n[i] < 2 || m_step[i] <= 0.0) {
            raise<ValueError>("PhotonLibraryModel: bad grid axis %d in '%s'", i, meta_path);
        }
    }
    m_nch = meta["nchan"].asUInt();

    std::string npy = meta["vis_npy"].asString();
    if (!npy.empty() && npy[0] != '/') {
        const auto slash = meta_path.find_last_of('/');
        if (slash != std::string::npos) {
            npy = meta_path.substr(0, slash + 1) + npy;
        }
    }
    cnpy::NpyArray arr = cnpy::npy_load(npy);
    if (arr.word_size != sizeof(float) || arr.fortran_order
        || arr.shape.size() != 4
        || arr.shape[0] != (size_t) m_n[0] || arr.shape[1] != (size_t) m_n[1]
        || arr.shape[2] != (size_t) m_n[2] || arr.shape[3] != m_nch) {
        raise<ValueError>("PhotonLibraryModel: vis npy '%s' shape/dtype mismatch with meta", npy);
    }
    const float* data = arr.data<float>();
    m_vis.assign(data, data + arr.num_vals);
}

void PhotonLibraryModel::visibilities(std::vector<double>& vis, const WireCell::Point& p_cm,
                                      const std::vector<unsigned int>* mask) const
{
    vis.assign(m_nch, 0.0);

    // grid index space, clamped to the boundary cell
    double f[3];
    int i0[3];
    const double p[3] = {p_cm.x(), p_cm.y(), p_cm.z()};
    for (int a = 0; a < 3; ++a) {
        double t = (p[a] - m_origin[a]) / m_step[a];
        t = std::clamp(t, 0.0, (double) (m_n[a] - 1));
        int i = std::min((int) t, m_n[a] - 2);
        i0[a] = i;
        f[a] = t - i;
    }

    const size_t sz = m_nch;                     // stride ch
    const size_t sy = sz * m_n[2];               // stride z
    const size_t sx = sy * m_n[1];               // stride y
    const size_t base = i0[0] * sx + i0[1] * sy + i0[2] * sz;

    const double w[8] = {
        (1 - f[0]) * (1 - f[1]) * (1 - f[2]),
        (1 - f[0]) * (1 - f[1]) * f[2],
        (1 - f[0]) * f[1] * (1 - f[2]),
        (1 - f[0]) * f[1] * f[2],
        f[0] * (1 - f[1]) * (1 - f[2]),
        f[0] * (1 - f[1]) * f[2],
        f[0] * f[1] * (1 - f[2]),
        f[0] * f[1] * f[2],
    };
    const size_t off[8] = {
        0, sz, sy, sy + sz, sx, sx + sz, sx + sy, sx + sy + sz,
    };

    for (size_t ch = 0; ch < m_nch; ++ch) {
        if (mask && ch < mask->size() && (*mask)[ch] == 0) continue;
        double v = 0.0;
        for (int c = 0; c < 8; ++c) {
            v += w[c] * m_vis[base + off[c] + ch];
        }
        vis[ch] = v;
    }
}
