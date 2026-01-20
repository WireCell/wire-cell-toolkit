#include "WireCellSpng/CellViews.h"

namespace WireCell::SPNG {
    CellViews::CellViews()
        : TensorSetFilter("CellViews", "spng") {}

    void CellViews::configure(const WireCell::Configuration& config)
    {
        // Propagate and parse config
        this->TensorSetFilter::configure(config);
        from_json(m_cfg, config);

        auto ianode = Factory::find_tn<IAnodePlane>(m_cfg.anode);
        m_wire_basis.clear();
        IChannel::vector chans;
        for (auto face_ident in config.face_idents) {
            auto iface = ianode->face(face_ident);
            m_wire_basis.push_back(CellBasis::cell_basis(iface));
        }

    }

    WireCell::Configuration CellViews::default_configuration() const
    {
        auto cfg = this->TensorSetFilter::default_configuration();
        auto cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }

    ITorchTensor::pointer CellViews::filter_tensor(const ITorchTensor::pointer& in)
    {
    }
    

}
