// this file collects destructors in cases where the class is too
// simple to warrant its own .cxx file.

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellAux/SimpleTrackSegmentSet.h"
#include "WireCellAux/SimpleWire.h"

using namespace WireCell::Aux;

SimpleBlob::~SimpleBlob() {}
SimpleBlobSet::~SimpleBlobSet() {}
SimpleDepoSet::~SimpleDepoSet() {}
SimpleTrackSegmentSet::~SimpleTrackSegmentSet() {}
SimpleWire::~SimpleWire() {}
