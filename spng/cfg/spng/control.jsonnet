// A "control" object collects various component config objects and basic
// parameters that are used by many components and which are independent of
// detector and particular job makeup.  Most of the other configs have a
// top-level function with an argument "control" that accepts the .config
// attribute.  Give the .uses attribute to pg.main().

local wc = import "wirecell.jsonnet";

function(device="cpu", seeds=[0,1,2,3,4], verbosity=0, semaphore=1)
{
    objects: {
        rng: {
            type: 'Random',
            name: 'rng' + std.join("",std.map(std.toString, seeds)),
            data: {
                generator: "default",
                seeds: seeds,
            }
        },
        /// Fixme: need a way to specify this.
        dft: {
            type: "FftwDFT",
            name: "",
            data: {}
        },
        
        semaphore: {
            type: 'Semaphore',
            name: "",               // for now, only allow one
            data: { concurrency: semaphore },
        },

    },
    uses: std.objectValues($.objects),

    /// The control config object
    config: {
        [key]: wc.tn($.objects[key])
        for key in std.objectFields($.objects)
    } + {
        device: device,
        verbosity: verbosity,
    },

}
