// this file collects destructors in cases where the class is too
// simple to warrant its own .cxx file.

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellAux/SimplePhotonHit.h"
#include "WireCellAux/SimplePhotonHitSet.h"
#include "WireCellAux/SimplePrimaryVertexSet.h"
#include "WireCellAux/SimpleSimTruth.h"
#include "WireCellAux/SimpleTrackSegmentSet.h"
#include "WireCellAux/SimpleTrajectory.h"
#include "WireCellAux/SimpleTrajectorySet.h"
#include "WireCellAux/SimpleWire.h"

using namespace WireCell::Aux;

SimpleBlob::~SimpleBlob() {}
SimpleBlobSet::~SimpleBlobSet() {}
SimpleDepoSet::~SimpleDepoSet() {}
SimplePhotonHit::~SimplePhotonHit() {}
SimplePhotonHitSet::~SimplePhotonHitSet() {}
SimplePrimaryVertexSet::~SimplePrimaryVertexSet() {}
SimpleSimTruth::~SimpleSimTruth() {}
SimpleTrackSegmentSet::~SimpleTrackSegmentSet() {}
SimpleTrajectory::~SimpleTrajectory() {}
SimpleTrajectorySet::~SimpleTrajectorySet() {}
SimpleWire::~SimpleWire() {}
