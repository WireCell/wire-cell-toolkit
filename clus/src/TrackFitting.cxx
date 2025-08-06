#include "WireCellClus/TrackFitting.h"
#include "WireCellUtil/Logging.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

using geo_point_t = WireCell::Point;


TrackFitting::TrackFitting(FittingType fitting_type) 
    : m_fitting_type(fitting_type) 
{

}

geo_point_t TrackFitting::adjust_rough_path(PR::Segment& segment){
    return geo_point_t(0,0,0);
}

void TrackFitting::collect_charge(double dis_cut, double range_cut){
    // collect charge within dis_cut and range_cut
}