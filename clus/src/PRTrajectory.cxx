#include "WireCellClus/PRTrajectory.h"

namespace WireCell::Clus::PR {

    Trajectory::Trajectory(graph_type& sg)
        : m_graph(sg)
    {
    }
    
    Trajectory::~Trajectory()
    {
    }


    void Trajectory::post_construct_init()
    {
        auto& gb = get_graph_bundle(m_graph);
        gb.trajectory = std::weak_ptr<Trajectory>(shared_from_this());
    }


    TrajectoryPtr get_trajectory(graph_type& sg)
    {
        auto& gb = get_graph_bundle(sg);
        // Return shared version of weak pointer.  This can be nullptr.
        return gb.trajectory.lock();
    }

}
