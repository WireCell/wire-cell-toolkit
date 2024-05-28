local nominal = import "nominal.jsonnet";
nominal {

    // The ductor transforms drifted depos to currents
    ductor: {

        # field_file: "dune-garfield-1d565.json.bz2",
        field_file: $.detector.fields,

        // The distance from the anode centerline to where the field
        // response calculation begins drifting.  Take care that field
        // response calculations are not always done with the
        // thickness of the anode in mind and may report their
        // "response plane" measured relative to a number of other
        // monuments.  This value is almost certainly required in the
        // "geometry" section.
        response_plane : 10*wc.cm,

        // The "absolute" time (ie, in G4 time) that the lower edge of
        // of final readout tick #0 should correspond to.  This is a
        // "fixed" notion.
        local tzero = -250*wc.us,

        // Simulate extra time so we can see the tail of response
        // from activity from depos that fall outside the readout
        // window.  This amount of extra is somewhat arbitrary but
        // the duration of the field response calculation is a
        // good measure.
        local pre = self.response_plane / $.lar.drift_speed,

        // The time bin where the readout should be considered to
        // actually start given the pre-signal simulated.
        tbin : wc.roundToInt(pre / self.binning.tick),
        binning: {
            // Ductor ticks at same speed as ADC
            tick : $.binning.tick,
            // Over the somewhat enlarged domain
            nticks : $.ductor.tbin + $.binning.nticks,
        },
        start_time : tzero - pre,
        readout_time : self.binning.nticks * self.binning.tick,
    },

}
