// Return object/set difference: big - small
function(small, big) 
    {[k]:big[k] for k in std.setDiff(std.objectFields(big), std.objectFields(small))}

          
                            
