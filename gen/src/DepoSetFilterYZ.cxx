#include "WireCellGen/DepoSetFilterYZ.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellIface/IDepo.h"
#include "WireCellUtil/Persist.h"
#include "WireCellIface/IAnodePlane.h"

WIRECELL_FACTORY(DepoSetFilterYZ, WireCell::Gen::DepoSetFilterYZ, WireCell::INamed, WireCell::IDepoSetFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Gen;

DepoSetFilterYZ::DepoSetFilterYZ()
  : Aux::Logger("DepoSetFilterYZ", "gen")
{
}
DepoSetFilterYZ::~DepoSetFilterYZ() {}

WireCell::Configuration DepoSetFilterYZ::default_configuration() const
{
    Configuration cfg;
    //../../util/src/Response.cxx

    cfg["yzmap_filename"] = "YZMap_filename";
    cfg["bin_width"]      = "BinWidth";
    cfg["tpc_width"]      = "TPCWidth";
    cfg["bin_height"]     = "BinHeight";
    cfg["resp"]           = "Response";
    cfg["anode"]          = "AnodePlane";
    cfg["plane"]          = "WirePlane";
    return cfg;
}

void DepoSetFilterYZ::configure(const WireCell::Configuration& cfg)
{

  //  log->debug("Conifg File Name");
  const std::string filename = cfg["yzmap_filename"].asString();
    if (filename.empty()) {
        THROW(ValueError() << errmsg{"DepoSetFilterYZ requires an YZ region"});
    }

    //    log->debug("Conifg Anode TN");    
    const std::string anode_tn = cfg["anode"].asString();
    if (anode_tn.empty()) {
        THROW(ValueError() << errmsg{"DepoSetFilterYZ requires an anode plane"});
    }
    //    log->debug("Config Andoe");
    WireCell::IAnodePlane::pointer anode = Factory::find_tn<IAnodePlane>(anode_tn);
    if (anode == nullptr) {
        THROW(ValueError() << errmsg{"Input anode is a nullptr"});
    }
    //    log->debug("Conifg Andode Face");
    IAnodeFace::vector abode_faces = anode->faces();
    for (auto face : abode_faces) {
        m_boxes.push_back(face->sensitive());
    }
    //    log->debug("Rest...");
    bin_width =  get<double>(cfg, "bin_width") * units::cm;
    tpc_width =  get<double>(cfg, "tpc_width") * units::cm;
    bin_height = get<double>(cfg, "bin_height") * units::cm;
    resp =       get<int>   (cfg, "resp");
    plane =      get<int>   (cfg, "plane");
    //std::cout << "Here is my issue? " << std::endl;
    anode_name = get<std::string>(cfg, "anode");
    //std::cout << "Nope!" << std::endl;
    jmap = WireCell::Persist::load(filename);

}

bool DepoSetFilterYZ::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("DepoSetFilterYZ fail with no input on call = {}", m_count);
        return true;
    }
    IDepo::vector output_depos;
    
    for (auto idepo : *(in->depos())) {
      log->debug("I GOT A DEPO IN DepoSetFilterYZ!");
        bool pass_resp = false;
	bool pass_anod = false;

	//double depo_x = idepo->pos().x();
	double depo_y = idepo->pos().y();
	double depo_z = idepo->pos().z();

	//int depo_bin_x = std::round(depo_x/tpc_width);
	int depo_bin_y = std::round(depo_y/bin_height);
	int depo_bin_z = std::round(depo_z/bin_width);

	log->debug(" DepoSetFilterYZ depo_bin_y: ", depo_bin_y);
	log->debug(" DepoSetFilterYZ depo_bin_z: ", depo_bin_z);

	for (const auto& node : jmap) {
	  
	  std::cout << "Next string : " << node["anode"].asString() << std::endl; 
	  if(node["anode"].asString() == anode_name &&
	     node["plane"].asInt() == plane){
	
	    auto map_resp = node["map"][depo_bin_y][depo_bin_z];
	    
	    if (map_resp.asInt() == resp) {
	      pass_resp = true;
	      break;
	    }
	  }
	}
	
	for (auto box : m_boxes) {
	  WireCell::Ray r = box.bounds();
	  
	  if (box.inside(idepo->pos())) {
	    pass_anod = true;
	    break;
	  }
        }	

        if (pass_resp && pass_anod) {
	  output_depos.push_back(idepo);
        }
    }
    log->debug("call={} Number of Depos for a give APA={}", m_count, output_depos.size());
    out = std::make_shared<WireCell::Aux::SimpleDepoSet>(m_count, output_depos);
    ++m_count;
    return true;
}
