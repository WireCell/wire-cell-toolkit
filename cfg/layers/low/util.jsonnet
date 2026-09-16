{
    // Various configuration objects, especially those for AnodePlanes
    // have a ident data attribute which is useful to use as a string.
    idents :: function(obj) std.toString(obj.data.ident),

    // Extract and return unique plane trio objects aka xregions from
    // an array of "drifts".
    //
    // A volume's faces array may contain null placeholders for insensitive
    // faces (e.g. single-drift VD CRPs).  Filter those explicitly: WCT's
    // bundled jsonnet std.prune raises on null array elements.
    driftsToXregions :: function(vols)
        std.set(std.filter(function(f) f != null,
                           std.flattenArrays([v.faces for v in vols])),
            function(o) std.toString(o)),

    assure_name :: function(name, obj) 
        if std.type(name) == "null"
        then $.idents(obj)
        else name,


}
