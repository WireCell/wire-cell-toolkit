// A "control" object collects various component config objects and basic parameter.

{
    random(seeds=[0,1,2,3,4]):: {
        type: 'Random',
        name: 'rng' + std.join("",std.map(std.toString, seeds)),
        data: {
            generator: "default",
            seeds: seeds,
        }
    },

    fftwdft: {
        type: "FftwDFT",
        name: "",
        data: {}
    },

    bundle(device="cpu", rng=$.random(), dft=$.fftwdft, verbosity=0):: {
        device: device, rng:rng, dft:dft, verbosity:verbosity
    },

}
