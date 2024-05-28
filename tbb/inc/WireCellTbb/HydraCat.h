/**
   The hydra category adapts a WCT IHydraNodeBase to TBB.  

   A TBB subgraph supports the WCT hydra node.  This subgraph can be
   considered split into input and output parts.

                 [input] - [output]
   [indexer(Nin) - func] - [mfunc - sequencers(Nout)]

   The internal "func" is used to break the Nin x Nout template
   parameter complexity to Nin + Nout.  The mfunc calls the WCT hydra
   node operator() and marshals any results out its Nout ports.
   Sequencers assure output ordering.

 */

#ifndef WIRECELLTBB_HYDRACAT
#define WIRECELLTBB_HYDRACAT

#include "WireCellTbb/NodeWrapper.h"
#include "WireCellIface/IHydraNode.h"
#include "WireCellUtil/Exceptions.h"

namespace WireCellTbb {

    /// The "any" interface to IHydraNode acccepts this type for both
    /// input and output.
    using any_queue_vector = WireCell::IHydraNodeBase::any_queue_vector;

    /// The TBB hydra implementation consists of an input side and an
    /// output side:

    /// Hydra input


    template<typename T, std::size_t N, typename... REST>
    struct generate_indexer_node
    {
        typedef typename generate_indexer_node<T, N-1, T, REST...>::type type;
    };
    template<typename T, typename... REST>
    struct generate_indexer_node<T, 0, REST...> {
        typedef tbb::flow::indexer_node<T, REST...> type;
    };

    template <typename IndexerNode, std::size_t... Is>
    receiver_port_vector indexer_ports(IndexerNode& dxn, std::index_sequence<Is...>)
    {
        return {dynamic_cast<receiver_type*>(&tbb::flow::input_port<Is>(dxn))...};
    }

    template<std::size_t Nin>
    using indexer_node = typename generate_indexer_node<msg_t, Nin>::type;

    template<std::size_t Nin>
    using indexer_msg_t = typename indexer_node<Nin>::output_type;

    // Erase the Nin.
    using tagged_msg_t = std::pair<size_t, msg_t>;

    template<std::size_t Nin>
    struct HydraInputBody {

        tagged_msg_t operator()(const indexer_msg_t<Nin>& in) const
        {
            const std::size_t index = in.tag();
            const msg_t& msg = in.template cast_to<msg_t>();
            return std::make_pair(index, msg);
        }

    };

    // Construct the hydra input of port-cardinality Nin.  This is a
    // subgraph of an indexer node followed by a function node to
    // erase the Nin template argument to break Nin X Nout complexity.
    // The subgraph is added to nodes in order of input to output.
    template <std::size_t Nin>
    tbb::flow::sender<tagged_msg_t>* // the output of the input
    build_hydra_input(tbb::flow::graph& graph,
                      receiver_port_vector& receivers,
                      std::vector<tbb::flow::graph_node*>& nodes)
    {
        auto* dxn = new indexer_node<Nin>(graph);
        nodes.push_back(dxn);
        receivers = indexer_ports(*dxn, std::make_index_sequence<Nin>{});

        auto* fn = new tbb::flow::function_node<indexer_msg_t<Nin>, tagged_msg_t>(graph, 1, HydraInputBody<Nin>{});
        nodes.push_back(fn);
        make_edge(*dxn, *fn);

        return fn;
    }


    /// Hydra output


    template<typename InT, typename T, std::size_t N, typename... REST>
    struct generate_outdexer_node
    {
        typedef typename generate_outdexer_node<InT, T, N-1, T, REST...>::type type;
    };
    template<typename InT, typename T, typename... REST>
    struct generate_outdexer_node<InT, T, 0, REST...> {
        typedef tbb::flow::multifunction_node<InT, std::tuple<T, REST...>> type;
    };


    template<typename OutdexerNode, std::size_t... Is>
    sender_port_vector outdexer_ports(OutdexerNode& mfn, std::index_sequence<Is...>)
    {
        return {dynamic_cast<sender_type*>(&tbb::flow::output_port<Is>(mfn))...};
    }

    template<std::size_t Nout>
    using mfunc_node_type = typename generate_outdexer_node<tagged_msg_t, msg_t, Nout>::type;

    template<std::size_t Nout>
    using mfunc_ports_type = typename mfunc_node_type<Nout>::output_ports_type;

    template<std::size_t Nout, std::size_t N=0>
    void try_transfer_n(mfunc_ports_type<Nout>& out, 
                        const any_queue_vector& oqv,
                        std::vector<size_t>& counts)
    {
        const auto& oq = oqv[N];
        for (const auto& o : oq) {
            bool accepted = std::get<N>(out).try_put(msg_t(counts[N]++, o));
            if (!accepted) {
                std::cerr << "TbbFlow: hydra node input return false when try put on " << N << "\n";
            }
        }
        if constexpr (N + 1 < Nout) {
            try_transfer_n<Nout, N + 1>(out, oqv, counts);
        }
    }

    template <std::size_t Nout>
    class HydraOutputBody {
        WireCell::IHydraNodeBase::pointer m_wcnode;
        mutable std::vector<seqno_t> m_seqnos;
        const size_t nin;
        NodeInfo& m_info;
      public:
        HydraOutputBody(WireCell::INode::pointer wcnode, NodeInfo& info)
            : m_seqnos(Nout, 0)
            , nin(wcnode->input_types().size()) // number of input ports of the hydra node
            , m_info(info)
        {
            m_wcnode = std::dynamic_pointer_cast<WireCell::IHydraNodeBase>(wcnode);
        }
        void operator()(const tagged_msg_t& in, mfunc_ports_type<Nout>& out)
        {
            // in: {iport, {seqno, any}}
            // iqv: vector (nports) of queues (?) of any
            any_queue_vector iqv(nin), oqv(Nout);

            size_t index = in.first;
            iqv[index].push_back(in.second.second);

            m_info.start();
            bool ok = (*m_wcnode)(iqv, oqv);
            m_info.stop();
            if (!ok) {
                std::cerr << "TbbFlow: hydra body return false ignored\n";
            }
            try_transfer_n<Nout>(out, oqv, m_seqnos);
        }
    };

    template <std::size_t Nout>
    tbb::flow::receiver<tagged_msg_t>* // the input to the output
    build_hydra_output(tbb::flow::graph& graph,
                       WireCell::INode::pointer wcnode,
                       sender_port_vector& senders,
                       std::vector<tbb::flow::graph_node*>& nodes,
                       NodeInfo& info)                       
    {
        auto mfn = new mfunc_node_type<Nout>(graph, 1, HydraOutputBody<Nout>(wcnode, info));
        nodes.push_back(mfn);

        auto spv = outdexer_ports(*mfn, std::make_index_sequence<Nout>{});

        for (size_t ind=0; ind<Nout; ++ind) {
            auto qn = new seq_node(graph, [](const msg_t& m) {return m.first;});
            nodes.push_back(qn);
            tbb::flow::make_edge(*spv[ind], *qn);
            senders.push_back(dynamic_cast<sender_type*>(qn));
        }
        return mfn;
    }


    // Wrap TBB around WCT hydra node
    class HydraWrapper : public NodeWrapper {
        std::vector<tbb::flow::graph_node*> m_nodes;
        receiver_port_vector m_rx;
        sender_port_vector m_tx;
      public:
        HydraWrapper(tbb::flow::graph& graph, WireCell::INode::pointer wcnode)
        {
            m_info.set(wcnode);
            const size_t nin = wcnode->input_types().size();
            const size_t nout = wcnode->output_types().size();

            tbb::flow::sender<tagged_msg_t>* fe;
            if(nin < 1 || nin > 3) {
                THROW(WireCell::ValueError() << WireCell::errmsg{"unsupported size for hydra input"});
            }
            if (1 == nin) fe = build_hydra_input<1>(graph, m_rx, m_nodes);
            if (2 == nin) fe = build_hydra_input<2>(graph, m_rx, m_nodes);
            if (3 == nin) fe = build_hydra_input<3>(graph, m_rx, m_nodes);
            // ... more if needed


            tbb::flow::receiver<tagged_msg_t>* be;
            if(nout < 1 || nout > 3) {
                THROW(WireCell::ValueError() << WireCell::errmsg{"unsupported size for hydra output"});
            }
            if (1 == nout) be = build_hydra_output<1>(graph, wcnode, m_tx, m_nodes, m_info);
            if (2 == nout) be = build_hydra_output<2>(graph, wcnode, m_tx, m_nodes, m_info);
            if (3 == nout) be = build_hydra_output<3>(graph, wcnode, m_tx, m_nodes, m_info);
            // ... more if needed

            // Join input front-end half to output back-end half
            tbb::flow::make_edge(*fe, *be);
        }

        virtual sender_port_vector sender_ports()
        {
            return m_tx;
        }

        virtual receiver_port_vector receiver_ports()
        {
            return m_rx;
        }          
    };
}

#endif
