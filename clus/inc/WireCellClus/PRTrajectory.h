/** Basic classes and functions to build and work with trajectory graphs
 */
#ifndef WIRECELL_CLUS_PR_TRAJECTORY
#define WIRECELL_CLUS_PR_TRAJECTORY

#include "WireCellClus/PRTrajectoryView.h"

namespace WireCell::Clus::PR {

    class TrajectoryView;
    using TrajectoryViewPtr = std::shared_ptr<TrajectoryView>;

    /** Manage a trajectory graph

        A PR::Trajectory provides a collection of PR::TrajectoryView which help
        to construct and query a user's trajectory graphs.

        A trajectory graph is owned by the user and is borrowed by reference by
        the Trajectory (and TrajectoryView).

        A PR::Trajectory (and views) provide a mechanism to associate additional
        methods and state to trajectory graph (or subgraph).

        Application code may create one or more PR::Trajectory (or view)
        instances on the user's graph.  An application class is recommended to
        keep a base Trajectory as a data member through the application class
        may extend PR::Trajectory (and PR::TrajectoryView) through inheritance.

     */
    class Trajectory {

    public:

        Trajectory(Graph& graph);

        /** Make and store a trajectory view of a given type.
         *
         * See PR::TrajectoryView.
         */
        template<typename ViewType>
        std::shared_ptr<ViewType> make_view() = {
            auto ptr = PR::make_view<ViewType>(m_graph);
            m_views.push_back(ptr);
            return ptr;
        };
            
        using ViewVector = std::vector<TrajectoryViewPtr>;

        /** Access all known views */
        const ViewVector& views() const { return m_views; }

        /** Number of views */
        size_t size() const { return m_views.size(); }

        /** Get an existing view as a particular type.
         *
         * Unlike STL at(), no exception is thrown and instead a nullptr is
         * returned on out-of-bounds index or type mismatch error.
         */
        template<typename ViewType>
        std::shared_ptr<ViewType> at(size_t index) = {
            if (index >= m_views.size()) return nullptr;
            auto base = m_views.at(index);
            return dynamic_pointer_cast<ViewType>(base);
        }

        
    private:
        Graph& m_graph;

        std::vector<TrajectoryViewPtr> m_views;
    };

}

#endif

