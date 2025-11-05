#ifndef WIRECELL_CLUS_PR_SHOWER_FUNCTIONS
#define WIRECELL_CLUS_PR_SHOWER_FUNCTIONS

namespace WireCell::Clus::PR {

    /** Modify shower assuming shower kinematics.
     *
     * This free function is is equivalent to the method of WCP's
     * WCShower::calculate_kinematics().
     */
    void shower_kinematics(ShowerPtr shower);

    void update_particle_type(ShowerPtr shower);
    
    void longmuon_kinematics(ShowerPtr shower);
}

#endif
