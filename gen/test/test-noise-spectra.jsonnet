// A contrived noise model suitable for GroupNoiseModel
// See aux/docs/noise.org for description.

local wc = import "wirecell.jsonnet";

// The Rayleigh PDF is shaped similar to real world noise
local rayleigh = function(x, sigma=0.1)
    x / (sigma*sigma) * std.exp(-0.5*std.pow(x/sigma, 2));

local gen = function(grpnum, nsamples, period, fpeak, rms, nsave) {
    local fsample = 1.0/period,
    local fnyquist = 0.5*fsample,
    local frayleigh = fsample/nsamples,
    
    local fsigma = fpeak*fnyquist,
    local nstep = std.floor(nsamples/(2*nsave)),

    // only up to Nyquist frequency, and then double.  Note, should -1 the range max
    local half_freqs = [ind*frayleigh for ind in std.range(0, nsamples/2)],
    local half_Rk = [rayleigh(freq, fsigma) for freq in half_freqs],

    local num = 3.141592*nsamples*nsamples*rms*rms,
    local den = 4*2*std.foldl(function(s,r) s + r*r, half_Rk, 0),
    local S = std.sqrt(num/den),
    nsamples: nsamples,
    period: period,
    group: grpnum,
    freqs: half_freqs[0:nsamples:nstep],
    amps: [S*R for R in half_Rk[0:nsamples:nstep]],
};
function(ngrps=10, nsamples=4096, nsave=64, period=0.5*wc.us, fpeak=0.1, rms=1*wc.mV)
[
    gen(grp, nsamples, period, fpeak, rms*grp/ngrps, nsave) for grp in std.range(1, ngrps)
]
