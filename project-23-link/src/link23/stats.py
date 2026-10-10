"""Summary statistics to compare a simulated distribution with the binomial."""
import numpy as np
from scipy import stats as sst


def binomial_pmf(n, p=0.75):
    k = np.arange(n + 1)
    return sst.binom.pmf(k, n, p)


def empirical_pmf(counts, n):
    return np.bincount(counts, minlength=n + 1) / len(counts)


def shannon_entropy(pmf, base=2):
    pmf = np.asarray(pmf)
    pmf = pmf[pmf > 0]
    return float(-(pmf * np.log(pmf)).sum() / np.log(base))


def kl_divergence(p, q):
    """KL(p || q) in bits, summed over the support of p."""
    p, q = np.asarray(p), np.asarray(q)
    m = p > 0
    return float((p[m] * np.log2(p[m] / q[m])).sum())


def total_variation(p, q):
    return 0.5 * float(np.abs(np.asarray(p) - np.asarray(q)).sum())


def bimodality_coefficient(counts):
    """Sarle's bimodality coefficient; values above ~0.555 hint at bimodality."""
    x = np.asarray(counts, dtype=float)
    n = len(x)
    g = sst.skew(x)
    k = sst.kurtosis(x)  # excess kurtosis
    return float((g**2 + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))))
