/** This file provides the "omni detector map".

It aggregates known omni objects by name.

See omni.jsonnet for the required API to produce your own.

Generally, experiment-specific omni objects are expected to be defined in
experiment/<name>/omni.jsonnet.  Each omni.jsonnet file may produce a single
omni object or an array of omni objects.

*/

local make_api(obj) =
    if std.type(obj) == "array"
    then {[o.name]:o for o in obj}
    else {[obj.name]:obj};
    
make_api(import "omni.jsonnet")  
+ make_api(import "omni/pdsp/omni.jsonnet")
+ make_api(import "omni/pdhd/omni.jsonnet")
