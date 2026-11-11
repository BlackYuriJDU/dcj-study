# dcj-study — Independent verification of prion-disease data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22164909.svg)](https://doi.org/10.5281/zenodo.22164909)

Open reanalysis of publicly deposited sporadic Creutzfeldt–Jakob disease (sCJD) datasets, run by a non-laboratory independent initiative in Brazil. The statistical core is written in standard-library Python; reading the data, drawing figures, building the PDF and running the tests need the pinned packages in `requirements.txt`, and the GWAS file (188 MB, not versioned) is downloaded by script with an MD5 check.

> **Correction, November 2026.** Earlier versions of this repository and of the preprint (up to v0.4) stated that the STX6 GWAS signal does not colocalize with the brain STX6 eQTL. That was the product of a formula error in our own code. The corrected analysis, checked variant by variant against the R package `coloc`, supports a shared causal variant (posterior 0.98 under the default priors; 0.83 with p12 = 1e-6), which reproduces already published work. Colocalization is not mediation. Every change, with before/after values, is listed in [`CORRECOES.md`](CORRECOES.md).

## Headline results

| Result | Number |
|---|---|
| GWAS GCST90001389 variants parsed | 6,314,492 (0 malformed) |
| Published loci independently replicated | **3/3** — PRNP p=1.6e-15 · GAL3ST1 p=6.2e-10 · STX6 p=7.5e-9 |
| Genomic inflation λ | 1.0587 (MAF-stratum gradient 0.016) |
| rs3747957 (*Brain* 2025 index) already in 2020 sumstats | p=9.7e-9, same β direction |
| STX6 lead-cluster posterior mass (r²≥0.8) | 90.5% |
| STX6: GWAS × brain eQTL colocalization (5 cohorts, 2,156 samples) | PP.H4 = 0.98 with default `coloc` priors; 0.83 with p12 = 1e-6 (single-causal-variant model; reproduces published work) |
| Brain expression GSE160208 vs original paper | published count of 184 DEGs reproduced with the magnitude filter applied before the FDR (183 with the FDR over the whole panel); r = 1.000 |
| Brain expression GSE160208, FDR < 0.05 | 441 of 742 estimable genes (437 of 800 when the 58 zero-variance genes are counted) |
| Blood miRNA GSE140069 after age/sex/RIN adjustment | 89 → 1 significant among 865 estimable miRNAs (84 → 1 when the 74 constant miRNAs are counted) |

## Repository map

```
ARQUIVO_COMPLETO.md      ← EVERYTHING: all reports + code inline (start here)
preprint/                ← English manuscript draft
research/                ← state of the art, dataset catalogs, ecosystem map
caso_referencia/               ← simulated reference case (training/didactic; not a real patient)
pipeline/lib/            ← statistical core (standard library only): colocalization, harmonisation, OLS, GWAS/eQTL access
pipeline/scripts/        ← analysis scripts (dependencies pinned in requirements.txt)
tests/                   ← tests that recompute every central result against independent references
CORRECOES.md             ← what changed in the November 2026 correction, with before/after values
pipeline/reports/        ← every generated report with intermediate numbers
colaboracao/             ← family guide (PT-BR), engagement letters, audit reports
memory/                  ← mistakes, decisions, successful patterns
```

## Reproducibility

```bash
make install             # virtual environment (.venv) with the pinned versions in requirements.txt
make test                # test suite; no test may be skipped
make dados-gwas          # downloads GCST90001389 (188 MB) from the GWAS Catalog and checks its MD5
make relatorios-offline  # brain expression, blood miRNA, blood × brain: versioned data only
make relatorios-gwas     # GWAS QC, clumping, colocalization (needs the GWAS file)
make relatorios-rede     # LD fine-mapping, offset check (needs Ensembl REST and NCBI E-utilities)
make relatorios          # all three
```

The GEO matrices and the eQTL Catalogue extracts are versioned; release, sample sizes and SHA-256
checksums of the eQTL extracts are in `pipeline/data/proveniencia_qtl.tsv`.

The statistics are implemented from first principles (Welch t-test, Benjamini–Hochberg FDR,
Wakefield Bayes factors, colocalization, least squares) and each central result is recomputed in
the test suite against an independent reference: `coloc.abf` output for colocalization, SciPy,
statsmodels and NumPy for the rest. Before November 2026 the colocalization "validation" compared
two copies of the same wrong formula; that is the failure the current tests are built to catch.

## Ethics

Public, de-identified data only. The "Caso Referência" dossier is a **simulated**
generic training profile, not any real person. No clinical claims are made about
any individual.

## License & citation

MIT — see `LICENSE`. Cite via `CITATION.cff` (DOI minted on first Zenodo release).
