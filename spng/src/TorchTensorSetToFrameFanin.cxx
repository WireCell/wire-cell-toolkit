#include "WireCellSpng/TorchTensorSetToFrameFanin.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellIface/INamed.h"

WIRECELL_FACTORY(TorchTensorSetToFrameFanin, WireCell::SPNG::TorchTensorSetToFrameFanin,
                 WireCell::INamed,
                 WireCell::SPNG::ITorchTensorSetToFrameFanin)


using namespace WireCell;

WireCell::Configuration SPNG::TorchTensorSetToFrameFanin::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    return cfg;
}

void SPNG::TorchTensorSetToFrameFanin::configure(const WireCell::Configuration& config)
{
    //Get the anode to make a channel map for output
    m_anode_tn = get(config, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);

    //Get output groups (map WirePlaneId --> output index)
    if (config.isMember("input_groups")) {
        auto groups = config["input_groups"];
        int i = 0;
        m_multiplicity = groups.size();
        for (auto group : groups) {
            for (auto wpid : group) {
                WirePlaneId the_wpid(wpid.asInt());
                m_input_groups[the_wpid] = i;
                m_per_group_channel_map.emplace_back();

                log->debug("WPID: {}", WirePlaneId(wpid.asInt()));
            }
            ++i;
        }
    }


    //Make a map to go from the channel ID to the output group
    //First Loop over the faces in the anode we're working with
    for (const auto & face : m_anode->faces()) {
        if (!face) {   // A null face means one sided AnodePlane.
            continue;  // Can be "back" or "front" face.
                       //Throw instead?
        }

        //Loop over each plane in the face
        for (const auto & plane : face->planes()) {

            //Each channel in the plane will be associated
            //with the corresponding WirePlaneId's output group
            for (const auto & channel : plane->channels()) {

                auto out_group = m_input_groups[plane->planeid()];
                //Within a given plane, the traces will be in order as they were
                //seen here. We have a map to determine what the output
                //size (nchannels) is when we make the tensors later.
                //We also save the reverse order map per each group for metadata
                auto output_channel = m_output_nchannels[out_group];
                m_per_group_channel_map[out_group][output_channel] = channel->ident();
                ++m_output_nchannels[out_group];
                // std::cout << "[hyu1]chmap: " << channel->ident() << " " << plane->ident() << " " << m_channel_map[channel->ident()] << std::endl;
            }
        }
    }
}

std::vector<std::string> SPNG::TorchTensorSetToFrameFanin::input_types()
{
    const std::string tname = std::string(typeid(ITorchTensorSet).name());
    // log->debug("Got {}", m_multiplicity);
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}


SPNG::TorchTensorSetToFrameFanin::TorchTensorSetToFrameFanin()
    : Aux::Logger("TorchTensorSetToFrameFanin", "spng") {}

bool SPNG::TorchTensorSetToFrameFanin::operator()(const input_vector& inv, output_pointer& out)
{
    out = nullptr;

    size_t neos = std::count(inv.begin(), inv.end(), nullptr);
    if (neos) {
        log->debug("EOS in {} of {} at call={}", neos, m_multiplicity, m_count);
        ++m_count;
        return true;
    }

    log->debug("multiplicity {} at call={}", m_multiplicity, m_count);
    auto itraces = std::make_shared<ITrace::vector>();
    //Loop over the TorchTensorSets in the TorchTensorSet vector
    for (size_t input_group = 0; input_group < inv.size(); ++input_group) {
        const auto & channel_map = m_per_group_channel_map[input_group];
        
        //Get the first tensor from the TorchTensorSet -- improve this later
        const auto & in = inv.at(input_group);
        auto tensor_clone = in->tensors()->at(0)->tensor().clone().to(torch::kCPU);

        auto accessor = tensor_clone.accessor<double, 3>();

        auto nrows = tensor_clone.size(1);
        auto ncols = tensor_clone.size(2);
        // reuse this temporary vector to hold charge for a channel.
        ITrace::ChargeSequence charge(ncols, 0.0);
        for (int i = 0; i < nrows; ++i) {
            const auto & channel = channel_map.at(i);
            for (int j = 0; j < ncols; ++j) {
                charge[j] = accessor[0][i][j];
            }
            auto trace = std::make_shared<Aux::SimpleTrace>(channel, 0/*tbin*/, charge);
            itraces->push_back(trace);
        }
    }

    log->debug("Make ident-less frame with {} traces at call={}", itraces->size(), m_count);
    out = std::make_shared<Aux::SimpleFrame>(0, 0., itraces, 0.);

    ++m_count;
    return true;
}
