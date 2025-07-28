#ifndef WIRECELL_CLUS_PR_SHOWER
#define WIRECELL_CLUS_PR_SHOWER

#include "WireCellClus/PRTrajectory.h"

namespace WireCell::Clus::PR {


    /** A shower is a type of trajectory with "shower like" appearance."
     */
    class Shower : public Trajectory {
    private:

        // FIXME: for real, initial dummy data just to distinquish this from base Trajectory.
        std::string m_shower_specific_data;

    protected:
        // Shower Root constructor (calls Trajectory's root constructor)
        Shower() : Trajectory(), m_shower_specific_data("Default Shower Data") {
            std::cout << "Shower (Root) constructor called." << std::endl;
        }

        // Shower Child constructor (calls Trajectory's child constructor)
        explicit Shower(graph_type&& graph)
            : Trajectory(std::move(graph)), m_shower_specific_data("Child Shower Data") {
            std::cout << "Shower (Child) constructor called." << std::endl;
        }

        // Custom Shower constructor (calls Trajectory's root constructor)
        Shower(const std::string& data) : Trajectory(), m_shower_specific_data(data) {
            std::cout << "Shower (Custom Root) constructor called with data: " << data << std::endl;
        }

        // Custom Shower constructor for child with specific data
        // Calls Trajectory's child constructor
        Shower(graph_type&& graph, const std::string& data)
            : Trajectory(std::move(graph)), m_shower_specific_data(data) {
            std::cout << "Shower (Custom Child) constructor called with data: " << data << std::endl;
        }

    public:
        virtual ~Shower() {
            std::cout << "Shower destructor called." << std::endl;
        }

        const std::string& get_shower_data() const {
            return m_shower_specific_data;
        }

        void print_type() const override {
            std::cout << "Type: Shower (Data: " << m_shower_specific_data << ")" << std::endl;
        }
    };

}
#endif
