'''
A proxy for the supernova pointing resolution of a TA selection.

The direction to a supernova is estimated from the electrons of the
selected SNB interactions.  Each reaction carries the neutrino direction
differently: elastic scattering on electrons (ES) gives forward electrons,
CC on 40Ar nearly isotropic ones (slightly backward).  The proxy combines
them as an optimally weighted mean of the measured electron directions u_i
around the true neutrino direction n:

- group the selected interactions by reaction and electron kinetic energy;
- per group g, from the measured directions, the mean cosine mu_g and the
  transverse variance per axis s_g = (1 - <cos^2>_g) / 2;
- with N_g selected interactions expected in a burst, the information is
  I = sum_g N_g mu_g^2 / s_g (the weights w_g = mu_g / s_g are optimal) and
  the angular error per transverse axis sigma = 1/sqrt(I);
- mu_g^2 is estimated without bias from the k_g simulated events, as
  mean^2 - variance/k_g (CC electrons are so nearly isotropic that the
  plain mean^2 would measure only the sample's statistical fluctuation);
- the proxy is the cone containing 68% of the estimates of a 2D Gaussian,
  theta68 = 1.515 sigma.

The measured direction is the true electron direction (MARLEY) smeared by a
reconstruction model: a 2D Gaussian kick of sigma(T) = sigma10 (10 MeV /
T)^power and a head-tail flip of probability p_flip.  The selection enters
only through which interactions are found and so it respects the
simulation: it uses the electron kinematics from MARLEY, not the toy
transport.

Requires only numpy.
'''

import os

import numpy

CC, ES = 0, 1
theta68_per_sigma = numpy.sqrt(-2 * numpy.log(1 - 0.68))  # 1.5096


class Reco:
    '''Reconstruction model of the electron direction.'''

    def __init__(self, sigma10_deg=0.0, power=1.0, p_flip=0.0, nrep=16, seed=1):
        self.sigma10 = numpy.radians(sigma10_deg)
        self.power = power
        self.p_flip = p_flip
        self.nrep = nrep if (sigma10_deg > 0 or p_flip > 0) else 1
        self.seed = seed

    def measured_moments(self, cos, te):
        '''Per event, the mean of cos and cos^2 of the measured direction to
        the neutrino over nrep smearings.'''
        cos = numpy.clip(numpy.asarray(cos, dtype=float), -1, 1)
        if self.nrep == 1:
            return cos, cos**2
        rng = numpy.random.default_rng(self.seed)
        n = len(cos)
        sig = self.sigma10 * (10.0 / numpy.maximum(te, 0.1))**self.power
        delta = sig[:, None] * numpy.sqrt(-2 * numpy.log(rng.uniform(size=(n, self.nrep)) + 1e-300))
        phi = rng.uniform(0, 2 * numpy.pi, size=(n, self.nrep))
        sin = numpy.sqrt(1 - cos**2)
        c = cos[:, None] * numpy.cos(delta) + sin[:, None] * numpy.sin(delta) * numpy.cos(phi)
        flip = rng.uniform(size=(n, self.nrep)) < self.p_flip
        c = numpy.where(flip, -c, c)
        return c.mean(axis=1), (c**2).mean(axis=1)


def summary_path(path):
    '''The MARLEY summary of a run output, from run-snb-ta-study.sh.'''
    base = path[:-4] if path.endswith(".npz") else path
    return base + "-marley.npz"


class Kinematics:
    '''The MARLEY kinematics of the signal interactions of a Study.

    LowEDepoSource takes the MARLEY events of a job in order, so the k-th
    signal interaction of a file (by id) is event k of its summary.'''

    def __init__(self, study):
        ints = study.ints
        n = len(ints["id"])
        self.cos = numpy.full(n, numpy.nan)
        self.te = numpy.full(n, numpy.nan)
        self.enu = numpy.full(n, numpy.nan)
        self.reaction = numpy.full(n, -1, dtype=int)
        sig = ints["signal"] == 1
        ifile = ints["file"].astype(int)
        for f, path in enumerate(study.paths):
            m = numpy.flatnonzero(sig & (ifile == f))
            if len(m) == 0:
                continue
            spath = summary_path(path)
            if not os.path.exists(spath):
                raise FileNotFoundError(f"no MARLEY summary {spath}")
            summ = numpy.load(spath)
            m = m[numpy.argsort(ints["id"][m])]
            ev = numpy.arange(len(m)) % len(summ["enu"])
            if not numpy.allclose(summ["enu"][ev], ints["energy"][m]):
                raise ValueError(f"MARLEY summary {spath} does not match the interactions")
            self.cos[m] = summ["cos"][ev]
            self.te[m] = summ["te"][ev]
            self.enu[m] = summ["enu"][ev]
            self.reaction[m] = summ["reaction"][ev]


class Pointing:
    '''Pointing information of selections of one Study.

    count: the expected number of interactions of the study's reaction(s)
    in the burst for the detector considered (eg a module), all of which
    arrive in the readout; the found fraction scales it.'''

    te_edges = (0, 5, 10, 15, 20, 30, 200)

    def __init__(self, study, kin=None, reco=None):
        self.study = study
        self.kin = kin if kin is not None else Kinematics(study)
        self.reco = reco if reco is not None else Reco()
        self.pool = study.inwin & (self.kin.reaction >= 0)
        c1, c2 = self.reco.measured_moments(numpy.nan_to_num(self.kin.cos[self.pool]),
                                            numpy.nan_to_num(self.kin.te[self.pool]))
        self.c1 = numpy.zeros(len(self.pool))
        self.c2 = numpy.zeros(len(self.pool))
        self.c1[self.pool] = c1
        self.c2[self.pool] = c2

    def information(self, found, count, reaction, weights=None):
        '''Information (1/rad^2) of the found interactions of a reaction,
        for count interactions in the burst.  weights (per interaction, eg
        Poisson for a bootstrap) default to one.'''
        w = numpy.ones(len(self.pool)) if weights is None else numpy.asarray(weights, dtype=float)
        m = self.pool & (self.kin.reaction == reaction)
        ntot = w[m].sum()
        if ntot == 0:
            return 0.0
        info = 0.0
        te = self.kin.te
        for lo, hi in zip(self.te_edges[:-1], self.te_edges[1:]):
            g = m & found & (te >= lo) & (te < hi)
            wg = w[g]
            k = wg.sum()
            if k < 2:
                continue
            c1 = self.c1[g]
            mean = (wg * c1).sum() / k
            var = (wg * (c1 - mean)**2).sum() / (k - 1)
            mu2 = mean**2 - var / k
            sv = 0.5 * (1 - (wg * self.c2[g]).sum() / k)
            info += count * (k / ntot) * mu2 / max(sv, 1e-6)
        return max(info, 0.0)

    def mean_cos(self, found, reaction):
        '''Mean measured cosine to the neutrino of the found interactions.'''
        m = self.pool & found & (self.kin.reaction == reaction)
        return float(self.c1[m].mean()) if m.any() else float("nan")

    def efficiency(self, found, reaction):
        m = self.pool & (self.kin.reaction == reaction)
        return (found & m).sum() / m.sum() if m.any() else float("nan")


def theta68_deg(information):
    '''The 68% cone of the pointing estimate given its information.'''
    if information <= 0:
        return float("inf")
    return float(numpy.degrees(theta68_per_sigma / numpy.sqrt(information)))


def es_per_cc(cc_xs, es):
    '''Number of ES (all flavors) per nu_e CC interaction.

    cc_xs: flux averaged nu_e CC cross section per atom; es: dict flavor ->
    (flux averaged ES cross section per atom, mean energy, multiplicity),
    with nu_e's mean energy first.  Assumes equal luminosity per flavor, so
    the flux of a flavor is proportional to 1 / mean energy.'''
    emean_e = es["ve"][1]
    return sum(xs * mult * emean_e / emean for xs, emean, mult in es.values()) / cc_xs
