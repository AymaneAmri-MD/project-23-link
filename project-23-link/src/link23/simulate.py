"""Meiosis and offspring simulation under the Haldane map function.

Model
-----
Both parents are heterozygous (Aa) at every locus, and the question is the
number X of loci at which an offspring shows the dominant phenotype
(carries at least one A). Without linkage X ~ Binomial(n, 0.75).

Along a chromosome, the parental haplotype transmitted in a gamete switches
between two adjacent loci with probability r = (1 - exp(-2 d)) / 2, where d is
their distance in Morgans (Haldane: crossovers are Poisson, no interference).
Different chromosomes always assort independently (r = 0.5).
"""
from collections import namedtuple

import numpy as np

Loci = namedtuple("Loci", ["chrom", "pos_cM"])


def haldane_r(d_cM):
    """Recombination fraction for a distance in cM (Haldane map function)."""
    d = np.asarray(d_cM, dtype=float) / 100.0
    return 0.5 * (1.0 - np.exp(-2.0 * d))


def make_loci(chrom, pos_cM):
    """Build a Loci tuple, sorted by (chromosome, position)."""
    chrom = np.asarray(chrom)
    pos = np.asarray(pos_cM, dtype=float)
    order = np.lexsort((pos, chrom))
    return Loci(chrom[order], pos[order])


def simulate_gametes(loci, n, rng, phase="coupling"):
    """Simulate n gametes from one Aa/Aa/... parent.

    Returns an (n, L) uint8 array: 1 if the gamete carries allele A at that locus.

    phase="coupling": one parental haplotype carries all A, the other all a
                      (strongest possible linkage signal).
    phase="random"  : each locus gets an independent random phase for this parent.
    """
    L = len(loci.chrom)
    out = np.empty((n, L), dtype=np.uint8)

    if phase == "coupling":
        phase_vec = np.zeros(L, dtype=np.uint8)
    elif phase == "random":
        phase_vec = rng.integers(0, 2, size=L, dtype=np.uint8)
    else:
        raise ValueError("phase must be 'coupling' or 'random'")

    for c in np.unique(loci.chrom):
        idx = np.where(loci.chrom == c)[0]
        start = rng.integers(0, 2, size=n)
        if len(idx) > 1:
            r = haldane_r(np.diff(loci.pos_cM[idx]))
            switches = rng.random((n, len(r))) < r
            steps = np.concatenate(
                [np.zeros((n, 1), dtype=int), np.cumsum(switches, axis=1)], axis=1)
            hap = (start[:, None] + steps) % 2
        else:
            hap = start[:, None]
        # hap == 0 -> parental haplotype 1 (carries A under coupling)
        out[:, idx] = (hap ^ phase_vec[idx][None, :]) ^ 1
    return out


def simulate_offspring_counts(loci, n, rng=None, phase="coupling"):
    """Number of dominant phenotypes per offspring for n offspring of Aa x Aa."""
    rng = np.random.default_rng(rng)
    g1 = simulate_gametes(loci, n, rng, phase)
    g2 = simulate_gametes(loci, n, rng, phase)
    return (g1 | g2).sum(axis=1)
