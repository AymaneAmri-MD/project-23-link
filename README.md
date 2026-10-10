# Project 23-Link

**Boundary conditions of the (3 + 1)ⁿ Mendelian model**: what happens to the ideal polyhybrid distribution when loci sit on real chromosomes and recombine at realistic rates.

> Author: Aymane Amri Amine, pre-med candidate, biomedical engineering and AI.
> The original note is in [`paper/`](paper/Project_23_Link_v1.pdf).

## Question

For *n* independent loci in an *Aa × Aa* cross, the number of dominant phenotypes is Binomial(*n*, 0.75). With *n* = 23 this gives μ = 17.25 and σ ≈ 2.08. How far does this ideal description hold once loci are physically linked on shared chromosomes?

## Open in Colab

| Notebook | Content | |
|---|---|---|
| `01_binomial_vs_gaussian` | Checks of the analytical numbers, Gaussian approximation and its tails, dominance counting vs additive polygenic model | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/USER/project-23-link/blob/main/notebooks/01_binomial_vs_gaussian.ipynb) |
| `02_linkage_simulation` | Meiosis simulation with the Haldane map function, linkage sweeps, effective number of independent blocks, LD decay | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/USER/project-23-link/blob/main/notebooks/02_linkage_simulation.ipynb) |

Replace `USER` with your GitHub username in the badges and in the `REPO_URL` variable of the first notebook cell.

## Main results

All numbers come from 200 000 simulated offspring per scenario (see the notebooks).

- **Control (1 locus per chromosome):** the simulation reproduces the binomial (variance 4.31, KL ≈ 5·10⁻⁵ bits, i.e. sampling noise). Independent assortment means no linkage by construction.
- **23 loci on one chromosome:** variance rises from 4.3 to about 98 at 0 cM spacing and still about 29 at 10 cM. At complete linkage only 0 or 23 dominant phenotypes occur.
- **Effective number of independent blocks:** with 23 loci 1 cM apart, spreading them over 1, 2, 3, 5, 12 chromosomes gives variances of about 82, 45, 31, 20 and 8. The chromosome acts as the unit of independence.
- **Gaussian tails:** the normal approximation underestimates P(X = 0) = 2⁻⁴⁶ ≈ 1.42·10⁻¹⁴ by a factor of about 40.

![Linkage histograms](figures/02_linkage_histograms.png)
![Variance vs number of chromosomes](figures/02_effective_blocks.png)

## Repository layout

```
README.md
requirements.txt
notebooks/   01_binomial_vs_gaussian.ipynb, 02_linkage_simulation.ipynb
src/link23/  genome.py (map lengths), simulate.py (meiosis, offspring), stats.py (KL, entropy, TV, bimodality)
figures/     figures exported by the notebooks
paper/       original note (PDF)
```

## Run locally

```bash
pip install -r requirements.txt
jupyter lab notebooks/
```

## Model and limitations

- Haldane map function r = ½(1 − e^(−2d)); no crossover interference; sex-averaged, rounded chromosome lengths (23 = 22 autosomes + X).
- Both parents are heterozygous at every locus; default phase is coupling (strongest linkage signal), `phase="random"` is available.
- This is a **family-cross** model. Population-level linkage disequilibrium (D′, r²) needs forward simulation (`msprime`, `stdpopsim`, SLiM) and real data such as 1000 Genomes. The LD curve in notebook 02 is analytic, not simulated.
- Pharmacogenomic response clusters are mostly driven by specific variants and diplotypes (CYP2D6, CYP2C19...), a natural next use case.

## Roadmap

1. Real genetic map (HapMap or deCODE) instead of rounded lengths.
2. Population LD with `msprime` / `stdpopsim`.
3. 1000 Genomes case study on a pharmacogene.

## License

MIT, see `LICENSE`.
