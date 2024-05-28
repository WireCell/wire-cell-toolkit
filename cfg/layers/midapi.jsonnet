// This file defines the "mid" layer API.
//
// It provides a "base class" with dummy methods that will assert if not
// overridden.
//
// Caution: due to laziness it is likely woefully incomplete.
//
// Each detector variant mid API object will override what is here.
// Thus, this object simply gives a place to define an "abstract base
// class".  Any method not provided by a concrete implementation and
// which gets called will result in an error.

{
    // Return a drifter configured for a detector.  For optimal
    // runtime performance, the implementation should return a
    // DepoSetDrifter and NOT a chain with singular Drifter followed
    // by Bagger.  Mid implementer should use low.drifter  
    drifter() :: std.assertEqual("drifter" == "implemented"),

    
}
