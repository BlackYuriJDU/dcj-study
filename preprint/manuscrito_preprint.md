# Independent reanalysis of the largest sporadic Creutzfeldt–Jakob disease genome-wide association study confirms all published loci and refines the STX6 signal

**Arthur Araújo da Silva** — Independent Researcher, Brazil
*Preprint draft v0.5 — 2026-11-11 — target: bioRxiv*

*Independent data-organization initiative. This work uses only publicly deposited, anonymized
datasets; it involves no patient-identifiable data and makes no clinical claim.*

> **Correction notice (v0.5, replaces v0.4).** Version 0.4 of this preprint stated that the STX6
> GWAS signal and the brain STX6 eQTL could not be attributed to a shared causal variant (H₃ ≈ 0.995,
> "H₄ unidentifiable"). **That statement was wrong, and it came from an error in our own code.**
> In the colocalization routine the H₄ term was multiplied by p₁·p₂·p₁₂ instead of p₁₂, which
> suppresses H₄ by a factor of 10⁸ under the default priors. Two further errors sat in the same
> function: √r in place of √(1−r) in the Wakefield Bayes factor, and a prior variance of 0.04²
> labelled as the `coloc` default (the default is 0.04 for case–control traits). The "validation
> against R" reported in v0.4 compared two copies of the same wrong equation; the reference
> `coloc.abf` output, which v0.4 explained away as "prior dependence", was the correct result.
> With the corrected implementation, verified variant by variant against `coloc.abf` 5.2.3, the
> same data give a posterior probability of **0.98** for a shared causal variant (0.83 with
> p₁₂ = 10⁻⁶). Other changes in v0.5:
> - GTEx sample sizes follow the eQTL Catalogue release-8 metadata (268 and 263; 2,156 samples in
>   total, not ~2,182); variants are harmonised by full identity (position, REF, ALT).
> - Tests that cannot be estimated (74 miRNAs and 58 genes with no variance) no longer receive a
>   p-value and are excluded from the FDR universe, which is now declared for every count. Headline
>   counts change from 84 to 89 (blood miRNAs, unadjusted) and from 437 to 441 (brain genes); the
>   v0.4 counts are reported alongside. The adjusted blood result is unchanged (1 miRNA).
> - The statement that age "accounts for ~97%" of the blood-miRNA collapse is withdrawn: a ratio
>   of counts of significant tests is not an attributable fraction.
> - "Exact replication" of the brain DEG count is replaced by what was actually reproduced: the
>   published count, under a stated order of filters.
> - The blood × brain target-enrichment test used the wrong universe; recomputed, it is still null.
> - Claims of a dependency-free, R-validated statistical core are replaced by a description of the
>   actual dependencies and of the reference tests now in place.
>
> A complete list with before/after values is in `CORRECOES.md` in the repository.

---

## Abstract

**Background.** Sporadic Creutzfeldt–Jakob disease (sCJD) is the most common human prion disease.
The largest sCJD genome-wide association study (GWAS) identified PRNP, STX6 and GAL3ST1 at
genome-wide significance in 4,110 cases and 13,569 controls; its summary statistics are public and
have since been used by other groups for multi-omic risk-gene prioritization and for colocalization
method development. We report an audit conducted independently of those analyses.

**Methods.** We reanalysed the deposited summary statistics (GCST90001389; 6,314,492 variants,
GRCh37) end to end, implementing the inferential statistics from first principles and checking
them against reference implementations (the R package `coloc`, SciPy, statsmodels):
genomic-control inflation, genome-wide significance screening, Wakefield approximate Bayes
factors, linkage-disequilibrium clustering against 1000 Genomes phase 3, stratified inflation by
minor-allele-frequency bin, and colocalization with brain expression QTLs. We additionally reanalysed two
companion expression datasets (GSE160208 brain; GSE140069 whole blood).

**Results.** We independently reproduced **all three published loci** — PRNP chr20:4,672,307
(p = 1.62×10⁻¹⁵), GAL3ST1 chr22:30,950,360 (p = 6.18×10⁻¹⁰), STX6 chr1:180,961,245
(p = 7.51×10⁻⁹) — recovered the published STX6 index variant rs3747957 at its exact reported value
(p = 9.74×10⁻⁹, β = −0.148), and found no malformed records; global λ = 1.0587, varying by only
0.0156 across allele-frequency strata. Fine-mapping shows the STX6 signal is one cohesive
haplotype block carrying 90.5% of regional posterior mass. A blind clumping scan rediscovered the
three known loci without locus hints (3/3), indicating pipeline sensitivity rather than circular
confirmation. A five-cohort brain eQTL meta-analysis (2,156 samples) makes the STX6 lead variant
a genome-wide significant eQTL (p = 7.60×10⁻⁴⁷) with 89% direction concordance, and
approximate-Bayes-factor colocalization gives a posterior probability of 0.98 for a causal variant
shared between sCJD risk and STX6 expression under the default `coloc` priors (0.83 with
p₁₂ = 10⁻⁶), assuming a single causal variant per trait. This reproduces the expression–disease
relationship previously established with larger multi-omic designs; colocalization is not
evidence of mediation. No significant cis-sQTL was detected for STX6 in any cohort; because only
clusters with a significant signal are available, this bounds what we can detect rather than
establishing absence. In GSE140069, the number of FDR-significant blood microRNAs falls from 89
unadjusted (84 FDR-significant miRNAs when the 74 non-estimable ones are counted as tests, as in
v0.4) to 1 after adjustment for age, sex and RNA integrity; adding age alone to the group-only
model reduces the count from 92 to 4.

**Conclusions.** Public prion-disease summary statistics support rigorous independent
verification by non-consortium parties. This study confirms the published loci and their local LD
structure, quantifies the ceiling of summary-statistics-only fine-mapping, reproduces the
published colocalization of the STX6 signal with brain STX6 expression, and shows how strongly a
blood-based biomarker signature depends on covariate adjustment. All code and reports are openly
available.

---

## 1. Introduction

Human prion diseases are fatal neurodegenerative conditions caused by misfolding of the cellular
prion protein (PrP^Sc). Sporadic CJD (sCJD) accounts for ~85% of cases. Beyond the PRNP codon 129
modifier, host genetic modifiers were established by the largest sCJD genome-wide association study
to date [1], which identified PRNP, STX6 and GAL3ST1 at genome-wide significance in 4,110 cases
and 13,569 controls.

That study's summary statistics are publicly deposited (GCST90001389), and they have already been
analysed by other groups: the multi-omic risk-gene study [4] integrated them with eQTL, pQTL, SMR
and colocalization, and the deposited statistics have served as test data for colocalization
method development [8]. We make no novelty claim for accessing these data. What remains useful is
an audit conducted *independently of both* — one that verifies the deposition on its own terms,
stress-tests the published conclusions with a separately implemented statistical core, and
quantifies — with numbers a reader can reproduce — how much of the reported signal survives
standard covariate adjustment in the companion biomarker datasets.

Independent reanalysis serves three purposes the original analysis cannot: (i) it verifies that
public deposition is complete and internally consistent; (ii) it stress-tests conclusions with
methods chosen independently of the original pipeline; (iii) it lowers the barrier for
citizen-science participation in rare-disease research, following the precedent of Minikel &
Vallabh [10].

Here we report a fully independent, from-scratch reanalysis of GCST90001389 and two companion
expression datasets (GSE160208 brain tissue; GSE140069 whole blood), executed without access to
individual-level genotypes and with all inferential statistics implemented from first principles.
Where our results reproduce published findings we say so and attribute them; where they add
something, we state precisely what is new and what is not.

## 2. Data

| Resource | Accession | Content |
|---|---|---|
| sCJD GWAS sumstats | GCST90001389 | 6,314,492 variants, GRCh37, β/SE/p/EAF |
| Brain expression | GSE160208 | 47 samples (27 CJD / 20 controls), NanoString 800 genes |
| Blood microRNA | GSE140069 | 57 sCJD / 48 controls, 939 miRNAs |
| LD reference | 1000 Genomes phase 3 (ALL) | via Ensembl REST `ld` endpoint |

Integrity: the GWAS file is verified against the MD5 published by the GWAS Catalog (`pipeline/scripts/obter_gwas.py`); the eQTL extracts carry SHA-256 checksums and release provenance (`pipeline/data/proveniencia_qtl.tsv`).

## 3. Methods

### 3.1 GWAS quality control
Streaming parse of the compressed sumstats (single pass, constant memory). Malformed-line count, allele-frequency sanity, χ² statistic per variant as (β/SE)², genomic-control factor λ = median(χ²)/0.454936 computed globally and within pre-registered MAF strata (<0.05; 0.05–0.25; 0.25–0.45; ≥0.45). Genome-wide threshold p < 5×10⁻⁸.

### 3.2 Locus annotation
GRCh37 gene coordinates retrieved live from Ensembl REST (never from memory) for PRNP, STX6, GAL3ST1 windows (±50 kb around lead).

### 3.3 Approximate-Bayes-factor fine-mapping
Per-variant Wakefield ABF with prior variance W = 0.04 on log(OR): ABF_i = √(SE²/(SE²+W))·exp(χ²ᵢW/(2(SE²+W))). Regional posteriors πᵢ ∝ ABFᵢ. Because individual genotypes are unavailable, we report **cluster-level credible mass**: variants are grouped by pairwise r² ≥ 0.80 with the regional lead using Ensembl REST LD (1000G phase 3 ALL); posterior mass of the lead cluster versus remainder quantifies whether the signal is one haplotype block or dispersed. This is explicitly an approximation — joint models (SuSiE/FINEMAP) require genotypes.

### 3.4 Expression analyses
GSE160208: Welch t-tests on log₂-transformed normalized intensities, frontal-cortex-only contrasts, Benjamini–Hochberg FDR. The original authors' criteria (p < 0.05, q ≤ 0.06, |log₂FC| > 1) are applied with the order of the filters as an explicit parameter (FDR over the whole test universe followed by the magnitude filter, or magnitude filter first), alongside our FDR < 0.05 criterion. GSE140069: log₂(x+1) transform; primary model = OLS log₂ ~ group + sex + age + RIN, fitted by QR decomposition (covariates from the series matrix; the original publication adjusted age); sensitivity analyses unadjusted and detection-filtered. Effect sizes as Cohen's d. In both datasets a feature with no variance is not estimable: it receives no p-value and is excluded from the FDR universe, which is stated with every count (counts over all measured features, as used in v0.4, are given alongside). Group and sex categories and numeric covariates are validated on input. The t distribution, Welch test, BH procedure and least-squares fit are checked in the test suite against SciPy, statsmodels and NumPy.

### 3.5 Colocalization with brain eQTLs
For the STX6 region (chr1:180.9–181.1 Mb GRCh37), we tested whether the GWAS signal and STX6 expression quantitative trait loci (eQTLs) share a causal variant, using approximate-Bayes-factor colocalization (Giambartolomei et al. 2014) as implemented in the R package `coloc` (`approx.bf.estimates`, `combine.abf`), which we reimplemented in log space. Effect-size priors follow `coloc`: prior variance 0.04 (SD 0.2) for the case–control GWAS and 0.0225 (SD 0.15·sdY with sdY = 1) for the eQTL, whose expression phenotype is normalised before mapping (the sdY implied by the summary statistics is 0.67–1.00 across cohorts). Per-variant hypothesis priors are p₁ = p₂ = 10⁻⁴ and p₁₂ = 10⁻⁵, with sensitivity to p₁₂ (5×10⁻⁶ and 10⁻⁶ at fixed p₁ and p₂) and to the effect prior. eQTL summary statistics were retrieved as remote tabix region queries (pysam/htslib) from the eQTL Catalogue release 8, which uniformly reprocessed GTEx v10 and independent brain cohorts: CommonMind DLPFC (n=586), ROSMAP DLPFC (n=560), BrainSeq DLPFC (n=479), GTEx v10 DLPFC (n=268) and GTEx v10 cerebellum (n=263); sample sizes are those of the release-8 dataset metadata and do not enter the Bayes factors, which are computed from β and SE. The extracts are versioned with checksums and were confirmed byte-identical to the live release on 2026-11-10. eQTL coordinates (GRCh38) were converted to GRCh37 by subtracting a single constant offset (−30,864 bp) determined from the STX6 gene boundaries. **We verified this offset against the Ensembl assembly-mapping service**, sampling points across the full STX6 span in GRCh37: the returned GRCh38 coordinate is exactly +30,864 bp at every point sampled, and every eQTL variant that pairs with a GWAS variant under the constant also pairs under the service's own mapping. A constant offset is exact only where no assembly-level indel intervenes, so this verification bounds rather than eliminates the concern — a reimplementation should still prefer a chain-file liftOver — but for the STX6 window it holds empirically, and the colocalization result does not depend on the approximation. Variants were harmonised by full identity (position, REF, ALT; indels trimmed; rows repeated per rsID collapsed), with the eQTL ALT allele as effect allele; the GWAS was matched on both alleles, inverting its sign where its effect allele is the eQTL REF (among the 390 matched variants: 0 sign inversions needed, 0 incompatible alleles, 0 matches on the complementary strand only). Because no single cohort is individually powered, we performed a fixed-effect inverse-variance-weighted meta-eQTL across the five cohorts (2,156 brain samples; variants present in at least three cohorts; only identical alleles combined), yielding a meta-eQTL signal at the GWAS lead variant of p = 7.60×10⁻⁴⁷. Sensitivity analyses: DerSimonian–Laird random effects, leave-one-cohort-out, exclusion of both GTEx tissues (which share donors) and single-cohort analyses. Implementation check: per-variant posteriors (SNP.PP.H4) from our code agree with `coloc.abf` (coloc 5.2.3, R 4.5.3) to within 2×10⁻¹⁵ on all 390 variants, the summary posteriors agree, and synthetic cases are checked against explicit enumeration of the causal configurations. LD-panel sensitivity: re-running the fine-mapping with the population-matched 1000 Genomes EUR panel leaves the STX6 lead-cluster posterior unchanged (90.5% at r²≥0.8), improves PRNP tagging (58.9% → 80.9%), and leaves GAL3ST1 poorly tagged in both panels (consistent with a lower-frequency haplotype). Splicing scan: leafcutter sQTL summary statistics for the same five cohorts (eQTL Catalogue r8 `cc` files, which contain only intron clusters with a significant cis signal; cluster identity given by genomic coordinates) were harmonised to the GWAS in the same way and tested by colocalization per cohort×cluster; clusters were attributed to genes by the catalogue's gene_id column (STX6 = ENSG00000135823, GRCh37 chr1:180,941,861–180,992,047, verified via Ensembl REST lookup).

### 3.6 Blind whole-genome clumping (discovery sensitivity)
To demonstrate that our pipeline would find the known loci rather than merely confirm them, we scanned all 6,314,492 variants blind (threshold p < 10⁻⁵) and applied greedy distance clumping (±500 kb; lead = smallest p within cluster), with no locus annotation used until after cluster definition. Known loci were matched post hoc.

## 4. Results

### 4.1 Sumstats integrity and inflation
Zero malformed records among 6,314,492 lines. Global λ_GC = 1.0587 over all variants (liminal, acceptable). Stratified λ, computed on a systematic one-in-ten sample of variants:

| MAF stratum | n (sampled) | λ |
|---|---|---|
| <0.05 | 162,618 | 1.0617 |
| 0.05–0.25 | 266,332 | 1.0579 |
| 0.25–0.45 | 164,460 | 1.0547 |
| ≥0.45 | 38,039 | 1.0703 |

Gradient 0.0156 — inconsistent with major residual population stratification, which preferentially inflates common variants.

### 4.2 Independent replication of all three loci
| Locus | Our best hit (GRCh37) | p | β | Annotation |
|---|---|---|---|---|
| PRNP | chr20:4,672,307 C>T | 1.62×10⁻¹⁵ | −0.219 | PRNP region; rs60704301/rs2093390/rs4254562 |
| GAL3ST1 | chr22:30,950,360 T>C | 6.18×10⁻¹⁰ | −0.169 | GAL3ST1 promoter/5′ region |
| STX6 | chr1:180,961,245 G>A | 7.51×10⁻⁹ | −0.149 | intragenic STX6 (Ensembl 180,941,861–180,992,047) |

All 41 genome-wide-significant variants fall within these three regions (16 at PRNP, 2 at GAL3ST1, 23 at STX6); none elsewhere. A blind whole-genome clumping scan (Section 3.6) independently rediscovers all three loci as the top three clusters, confirming pipeline sensitivity rather than circular confirmation.

### 4.3 Cluster-level fine-mapping
- **PRNP**: 337 regional variants; anchor rs60704301 (merged → rs2093390). Posterior mass 100% within r²≥0.50 of the anchor in **both LD panels** (58.9% at r²≥0.80 with the ALL panel; **80.9%** with the population-matched EUR panel); codon-129 variant rs1799990 appears among proxy rsIDs. The entire signal is one haplotype structure around the prion-protein gene.
- **STX6**: 162 variants; anchor = lead rs11586493; 20/20 top variants panel-covered, max r² = 1.00; **90.5%** of mass in the r²≥0.80 lead cluster — a single cohesive block that includes rs3747957 (r² = 0.99); **identical (90.5%) under the EUR panel**.
- **GAL3ST1**: 322 variants; lead absent from 1000G phase 3 (anchor rs386462923→rs8142452, rank 3); only 4/20 top variants panel-covered; best cross-r² ≈ 0.45 — **unchanged under the EUR panel**. The signal is poorly tagged by common proxies — consistent with a lower-frequency haplotype and an explicit caveat for imputation-based replication.

### 4.4 Independent recovery of the STX6 index variant from the deposited statistics
The index variant at STX6 is **rs3747957** (chr1:180,953,853 GRCh37). This variant was reported by
the original consortium [1] — it appears in the abstract of that paper, in Table 1, in the
independent replication cohort (1,098 cases / 498,016 controls), in their eCAVIAR fine-mapping and
PAINTOR functional annotation, and was linked there to a GTEx putamen eQTL (p = 2.3×10⁻¹³). It was
subsequently carried forward as the STX6 anchor by an independent multi-omic study [4]. **Neither
of those findings is claimed as novel here.**

What this reanalysis contributes is an *independent reconstruction from the deposited summary
statistics alone*. Recovering rs3747957 directly from the sumstats — with no access to individual
genotypes and no reference to either source paper's tables — returns p = 9.74×10⁻⁹ and
β = −0.148, matching the published values to the reported precision, with the same effect
direction as our lead variant (β = −0.149); it ranks 11th of 162 regional variants by p-value.
The value of this recovery is not the association itself but the LD architecture it permits us to
quantify: rs3747957 sits at r² = 0.99 from our lead, inside a single cluster carrying 90.5% of
regional posterior mass (§4.3). Reconstructing that block structure from deposited statistics —
rather than from the original authors' fine-mapping output — is what allows an independent party
to confirm that the STX6 signal is one cohesive haplotype rather than a set of independent
false positives, and to quantify how much of the regional posterior mass any given analysis
method leaves unattributed.

### 4.5 Expression signatures: the brain count is reproduced under a stated filter order; the blood signature is fragile under covariate adjustment
Brain GSE160208: 58 of the 800 panel genes have no variance in frontal cortex and are not estimable; 441 of the 742 estimable genes are FDR < 0.05 (315 up, 126 down; 437 when the FDR is computed over all 800 genes, as in v0.4). Under the original authors' criteria (p < 0.05, q ≤ 0.06, |log₂FC| > 1) the count depends on the order of the filters: it is 184, the published count, when the magnitude filter is applied before the FDR, and 183 when the FDR is computed over the whole test universe (CCL4 is the borderline gene). We therefore report this as reproducing the published count under a stated filter order; an earlier comparison with the authors' supplementary gene list (August 2026, not repeated for this version) found the same top-ranked genes in each direction. The authors' q-values were not reproduced, and this is not an exact replication of their procedure. Blood GSE140069: 74 of the 939 miRNAs are constant across all samples and are not estimable. Among the 865 estimable miRNAs, unadjusted log₂-Welch yields 89 FDR-significant miRNAs (12↑/77↓; 84, 10↑/74↓, in the v0.4 universe of 939); OLS adjusting age+sex+RIN leaves **1** (hsa-miR-500a-3p) in either universe. That survivor is a detected/not-detected effect: the miRNA is above the detection floor in 27 of 57 cases and 46 of 48 controls, so its coefficient of −7.9 log₂ units reflects the detection rate rather than a several-hundred-fold change in abundance. The four discovery miRNAs retain direction and nominal significance (p = 7×10⁻⁴–4×10⁻²) but only hsa-miR-93-5p survives FDR within the detection-filtered universe (q = 0.048). Cases were on average 12.8 years older than controls (66.4 vs 53.6, Cohen's d = 1.07) with lower RNA integrity (RIN 5.59 vs 6.50, d = −0.71) — a textbook confounding structure that the original paper partially addressed (age via Partek GSA). We also fitted each covariate separately (OLS, 865 estimable miRNAs): the group-only model gives 92 FDR-significant miRNAs; adding age alone leaves 4, RNA integrity alone 28, sex alone 94, and all three 1 (v0.4 universe: 88, 4, 26, 85 and 1). All variance-inflation factors are below 1.4, so the coefficients are separately estimable. These are counts of tests crossing a threshold under different models. They show that the unadjusted signal does not survive adjustment for age and that RNA integrity alone already removes most of it; they are not additive, and they do not measure a proportion of the effect attributable to age. The corresponding attribution in v0.4 is withdrawn. Age-restricted analyses within age-overlapping strata (50–70 years, residual gap +2.7 years) leave 5–7 miRNAs, but are underpowered at 22–29 controls per stratum and bound rather than resolve any residual effect. Finally, an exploratory cross-compartment test (not performed by either original study) asks whether experimentally validated human targets of the four discovery blood miRNAs and of hsa-miR-500a-3p are over-represented among the up-regulated brain genes. With each miRNA's targets restricted to genes on the testable panel (N = 742, of which 315 are up-regulated), three miRNAs are testable: they have 13, 6 and 3 targets on the panel, of which 6, 3 and 1 are up-regulated (hypergeometric p = 0.50, 0.51 and 0.81); the other two have no validated target on the panel and cannot be tested. With so few targets and 42% of the panel up-regulated, the test has little power, and the null result neither supports nor argues against a relationship between the blood miRNAs and the cerebral transcriptional programme. (Version 0.4 ran this test with all genome-wide targets against N = 800, which gave incorrect expected counts and p-values.)

### 4.7 Colocalization: support for a causal variant shared by the STX6 GWAS signal and the brain STX6 eQTL
Colocalization within a single GTEx tissue is power-limited (largest eQTL |z| = 3.8 in either
tissue; concordant effect direction at 78–81% of variants). The five-cohort meta-eQTL is more
informative: it makes the lead GWAS variant a **genome-wide significant STX6 eQTL
(p = 7.60×10⁻⁴⁷; z = −14.4; five cohorts, I² = 48%)** with **89% effect-direction concordance**
across 390 harmonised variants (348/390). The lead GWAS variant is not the strongest eQTL in the
block: the smallest meta-eQTL p-value is at chr1:180,951,692 (p = 3.44×10⁻⁴⁷), and thirteen
further variants, the GWAS lead among them, have p-values within a factor of three of it.

Under the default `coloc` priors (p₁ = p₂ = 10⁻⁴, p₁₂ = 10⁻⁵) the posterior probabilities are
H₀ ≈ H₁ ≈ 0, H₂ = 7×10⁻⁶, **H₃ = 0.020 and H₄ = 0.980**. H₁ and H₂ are association with only one
of the two traits; H₃ is association with both through distinct causal variants; H₄ is association
with both through one shared causal variant. The analysis therefore supports a shared causal
variant, under the model's assumption of at most one causal variant per trait in the region.
Three qualifications belong with that number.

*Sensitivity to p₁₂.* p₁₂ is the only prior that separates H₄ from H₃. With p₁ = p₂ = 10⁻⁴, H₄ is
0.980 at p₁₂ = 10⁻⁵, 0.960 at 5×10⁻⁶ and 0.828 at 10⁻⁶. (The scenario p₁ = p₂ = 10⁻⁶, p₁₂ = 10⁻⁷,
which v0.4 called "conservative", raises the ratio p₁₂/(p₁p₂) from 10³ to 10⁵ and therefore
favours H₄; it gives 0.9998 and is not used as evidence.)

*Numerical robustness.* H₄ stays at 0.979–0.980 for effect-prior variances from 0.0016 to 0.1;
is 0.975 with DerSimonian–Laird random effects (Cochran's Q has p < 0.05 at 36 of the 390
variants, and I² > 50% at 51); 0.979–0.980 leaving out any one cohort; and 0.979 without both GTEx
tissues. CommonMind, ROSMAP and BrainSeq each give 0.978–0.979 on their own. The two GTEx tissues,
which share donors and are not independent cohorts, give 0.70 (DLPFC, n = 268) and 0.57
(cerebellum, n = 263) on their own, with most of the remaining mass on H₁ (H₃ ≈ 0.05): a power
limitation of small single-tissue eQTL samples, not support for distinct variants.

*How much comes from the prior.* The associated variants form a block in near-perfect LD. If all
390 variants were statistically equivalent, the prior alone would give H₄/H₃ =
p₁₂/(p₁p₂(n−1)) ≈ 2.6; the observed posterior ratio is ≈ 48. Part of the separation between H₃ and
H₄ therefore comes from the prior, and the data sharpen it about 19-fold. Conditional on H₄, the
95% credible set for the shared variant contains 17 variants; the variant with the highest
posterior (chr1:180,951,692; 0.105) is not the GWAS lead (0.063).

*What this analysis does not show.* Colocalization tests for a shared causal variant. It does not
test mediation, and it does not show that STX6 expression lies on the causal path between the
variant and disease. It assumes a single causal variant per trait, does not model LD, and uses
adult DLPFC and cerebellum, which only approximate the tissues and cell types affected in sCJD.
The result reproduces work already reported independently [1,4], which linked the STX6 risk locus
to STX6 expression using larger designs that include protein-level QTLs; we make no claim of
novelty.

*Implementation check.* Per-variant posteriors from our implementation match the R package
`coloc.abf` (coloc 5.2.3) to within 2×10⁻¹⁵ on all 390 variants, and the summary posteriors match
(H₄ = 0.9796, H₃ = 0.0204). Version 0.4 reported H₃ = 0.995 from an erroneous formula; see the
correction notice at the top of this document.

### 4.8 Splicing: no significant cis-sQTL for STX6 detected in any of the five brain cohorts
The eQTL Catalogue release 8 provides leafcutter splicing-QTL summary statistics for the same five
brain cohorts. We scanned every significant intron-cluster whose cis window intersects the STX6
region, harmonised variants to the GWAS and ran the colocalization per cohort × cluster
(6 tests in total). **No significant cis-sQTL cluster for STX6 was detected in any of the five
cohorts** (0/5), in contrast to the strong expression QTL at the same lead variant
(p = 7.60×10⁻⁴⁷). The six tests concern intron clusters of neighbouring genes (KIAA1614, QSOX1
and ENSG00000230124), not of STX6. None reaches PP.H₄ > 0.8. The highest is 0.60, for a KIAA1614
cluster in GTEx cerebellum; it falls to 0.13 with p₁₂ = 10⁻⁶, has concordant effect direction at
only 21% of variants, and is not reproduced in the three other cohorts that carry the same
cluster (0.19, 0.05 and below 0.001). We do not interpret it as colocalization. (Version 0.4
reported all six as PP.H₄ ≈ 0, an artefact of the formula error, and applied a Bonferroni
threshold to posterior probabilities, which is not meaningful.)

The strength of this negative result is limited by the data available to it, and we do not claim
otherwise. The `cc` files contain only clusters with a significant cis signal: a cluster without
one is **invisible to this scan**, so a true but weak or underpowered sQTL
would produce the same output as a genuinely absent one. Absence of evidence in
significance-filtered files is not evidence of absence. The defensible statement is therefore
narrower than "splicing does not operate here": *among the clusters with a significant signal in
these five cohorts, we detected no splicing signal at STX6 to interpret.* A claim about the
relative importance of expression versus splicing at this locus would require either unfiltered
sQTL statistics or a cohort with adequate power for splicing effects.

### 4.9 Blind clumping rediscovers all published loci
The blind scan produced 35 independent clusters; the top three are exactly PRNP (chr20:4,672,307;
p = 1.62×10⁻¹⁵), GAL3ST1 (chr22:30,950,360; p = 6.18×10⁻¹⁰) and STX6 (chr1:180,961,245;
p = 7.51×10⁻⁹) — **3/3 known loci rediscovered without any locus hint**, in correct significance
rank. The strongest residual cluster (chr16:15,539,902, near BMERB1; p = 5.73×10⁻⁸) sits just below
the genome-wide threshold. **This is not a new locus**: the signal was recorded as subthreshold by
the original consortium [1] and the same variant was functionally prioritized by the independent
multi-omic analysis [4]. We re-observe it because an unbiased scan recovers it in the residual,
which is a useful check on the analysis rather than a discovery; it still requires independent
replication in a cohort with adequate power.

### 4.10 Ethics and scope
No individual-level human data were generated or obtained beyond public deposits; no patient-identifiable information was processed. This work is a verification contribution and makes no clinical claims.

## 5. Discussion

Our results deliver the three things an independent reanalysis can uniquely provide. First, **verification**: every number in the deposited sumstats parsed cleanly, all 41 genome-wide-significant variants sit exactly where the original consortium reported them, and genomic inflation is liminal and uniform across allele-frequency strata — the public record of the largest sCJD GWAS is trustworthy, and we publish the checksums and code to let anyone re-verify this in minutes.

Second, **an independent reproduction of the STX6 expression association, and of its
colocalization, using a different eQTL resource**. That genetically increased STX6 expression and
protein abundance in brain raise risk of sCJD is established work [4], which integrated eQTL and
pQTL data with SMR and colocalization across 43 authors; the same direction of association has
also been reported for Alzheimer disease protein abundance [9]. We do not claim this as a
discovery. What we add is a reconstruction from a different source: a five-cohort brain meta-eQTL
built from eQTL Catalogue release 8 (2,156 samples), which makes the GWAS lead variant a
genome-wide significant STX6 eQTL (p = 7.60×10⁻⁴⁷) with concordant effect direction in 89% of
regional variants, and an approximate-Bayes-factor colocalization that gives a posterior
probability of 0.98 for a shared causal variant under the default priors — 0.83 with a ten-fold
smaller p₁₂, and stable to the effect prior, the meta-analysis model and the cohort set (§4.7).
This is corroboration of the published signal from an independent eQTL resource, not a new
mechanistic claim: colocalization supports a shared variant under a single-causal-variant model
and is not evidence of mediation. Version 0.4 of this preprint reached the opposite reading
because of an error in our own colocalization formula (see the correction notice); we record
that plainly, because a reader of v0.4 would have been misled.

Third, **an independent reconstruction of the STX6 haplotype block from deposited statistics
alone**. Using only the summary statistics and no access to the original authors' fine-mapping
output, we recover rs3747957 [1,4] at its published p = 9.74×10⁻⁹ and effect direction, and show it
lies at r² = 0.99 from our own lead inside a single cluster carrying 90.5% of regional posterior
mass (§4.3, §4.4). This is a verification result rather than a discovery: the association and the
index variant were both already published [1,4]. What it establishes is that the deposited
statistics are complete and internally consistent enough to reproduce a published locus and its
local LD structure from scratch — and it quantifies, for the first time in this locus, how much
regional posterior mass a summary-statistics-only analysis can attribute (90.5%) versus how much
remains unattributed, which is the honest ceiling of cluster-level inference.

Fourth, **a cautionary quantification for biomarker research**. In the blood microRNA data for sCJD (Nat Commun 2020), the number of FDR-significant miRNAs falls from 89 without covariate adjustment (84 in the v0.4 test universe) to 1 after adjustment for age, sex and RNA integrity. Adding each covariate separately to the group-only model (92 significant) leaves 4 with age, 28 with RNA integrity and 94 with sex. Cases were 12.8 years older than controls (d = 1.07), so age and disease status are partly confounded by design. The counts show that the unadjusted signature does not survive adjustment for age; a comparison of counts across models is not an estimate of how much of the effect age explains, and we no longer state one. Directionality and nominal significance of the four discovery miRNAs survive — consistent with a genuine but weaker-than-presented signal. We stress this is not an accusation: the original authors adjusted age themselves; our contribution is making the divergence between unadjusted and adjusted analyses explicit and reproducible for future biomarker pipelines.

Methodologically, the statistical core of this project (Bayes factors, colocalization, the t distribution, Welch's test, Benjamini–Hochberg and least squares) is written in standard-library Python, which keeps it small and inspectable. The project as a whole is not dependency-free: it needs openpyxl to read one dataset, matplotlib and reportlab for figures and the PDF, pysam/htslib for remote eQTL queries, and NumPy, SciPy and statsmodels as test references, all pinned in `requirements.txt`; the GWAS file is not versioned and is fetched by script with an MD5 check. The lesson of the present correction is that a from-scratch core is only as reliable as the references it is tested against: the v0.4 colocalization passed a "validation" that compared two copies of the same wrong formula. The test suite now checks each central result against an independent reference (`coloc.abf` output, SciPy, statsmodels) or against explicit enumeration. With that safeguard, the approach still lowers the entry barrier for researchers in resource-limited settings, including in countries like Brazil where sCJD surveillance exists but prion-genetics capacity is thin.

Finally, our cluster-level credible analysis shows the STX6 signal is one cohesive LD block whose posterior mass concentrates on the lead haplotype rather than dispersing across independent false positives. Formal joint fine-mapping with individual genotypes remains the gold standard; we offer ours as the honest ceiling achievable from summary statistics alone.

## 6. Limitations
- Summary-statistics-only: no conditional/joint fine-mapping; cluster-level inference only.
- Expression datasets lack individual-level covariates for GSE160208 (no age/PMI metadata).
- Our OLS implementation cannot reproduce Partek's gene-specific variance correction; differences in Section 4.5 may partly reflect estimator choice.
- Colocalization assumes at most one causal variant per trait in the region and does not model LD. The H₄ posterior is sensitive to p₁₂ (0.98 with the default 10⁻⁵; 0.83 with 10⁻⁶), and within a block in near-perfect LD part of the H₃/H₄ separation comes from the prior. A high H₄ supports a shared causal variant; it is not evidence of mediation.
- The two GTEx tissues share donors and are not independent cohorts; on their own they give H₄ of 0.70 and 0.57. The check against `coloc.abf` validates the implementation, not the model.
- This version corrects an error that reversed the colocalization conclusion of v0.4. No statistical geneticist has reviewed the corrected §4.7.
- Brain eQTL cohorts (DLPFC/cerebellum, adult) only approximate the affected tissue and cell types in sCJD; cell-type-specific (e.g. microglial) eQTL may differ.
- Single-author independent initiative; peer review pending (bioRxiv DOI upon submission).

## 7. Data availability
**GWAS summary statistics.** GCST90001389 (6,314,492 variants, GRCh37; GWAS Catalog, EMBL-EBI).
**Expression datasets.** GSE160208 (brain, NanoString 800-gene panel) and GSE140069 (whole blood
small-RNA-seq), both from the NCBI Gene Expression Omnibus.
**LD reference.** 1000 Genomes Project Phase 3, accessed via the Ensembl REST `ld` endpoint.

## 8. Code availability
All analysis code, reports containing every intermediate number, figure-generation scripts and
input checksums are version-controlled at https://github.com/BlackYuriJDU/dcj-study
(release `v1.0.0-preprint-clean`; archived on Zenodo; the corrections described in this version
are recorded in `CORRECOES.md`). The inferential statistics — approximate Bayes factors and
colocalization, the t distribution, Benjamini–Hochberg FDR, Welch's test, ordinary least squares
— are implemented from first principles in standard-library Python and checked in the test suite
against independent references: `coloc.abf` 5.2.3 output for colocalization (§4.7), and SciPy,
statsmodels and NumPy for the rest. External libraries: openpyxl (tabular input), matplotlib
(figures), reportlab (PDF), pysam/htslib (remote tabix access to eQTL summary statistics), and
NumPy, SciPy and statsmodels (tests). Versions are pinned in `requirements.txt`. The GWAS
summary-statistics file is not versioned: `pipeline/scripts/obter_gwas.py` downloads it from the
GWAS Catalog and verifies its MD5. `make relatorios` regenerates the analysis reports from the
committed scripts; the steps that need the GWAS file or network access (Ensembl, NCBI) are
separate targets.

## 9. References
1. Jones E, Hummerich H, Viré E, et al. Identification of novel risk loci and causal insights for
   sporadic Creutzfeldt–Jakob disease: a genome-wide association study. *Lancet Neurol*
   2020;19(10):840–848. doi:10.1016/S1474-4422(20)30273-8 (PMID 32949544; PMCID PMC8220892).
2. Areškevičiūtė A, Litman T, Suteika G, et al. Neuropathological markers of brain
   microRNAs linked to Creutzfeldt–Jakob disease. *Int J Mol Sci* 2020;22(1):140.
   doi:10.3390/ijms22010140 (PMID 33375642).
3. Norsworthy PJ, Fromm C, Schindler M, et al. A blood miRNA signature associates with sporadic
   Creutzfeldt–Jakob disease diagnosis. *Nat Commun* 2020;11:3960.
   doi:10.1038/s41467-020-17655-x (PMID 32769986).
4. Küçükali F, Hill E, Watzeels T, et al. Multiomic analyses direct hypotheses for Creutzfeldt–
   Jakob disease risk genes. *Brain* 2025;148(9):3350–3363. doi:10.1093/brain/awaf032
   (PMID 39865673; PMCID PMC12404779).
5. Minikel EV, Prusiner SB. Analysis of Guizhou Prion Disease in a Large-Scale Study of
   Long-Lived CJD Patients. *Sci Transl Med* 2016;8(339):340ra73. doi:10.1126/scitranslmed.aad4584.
6. Wakefield J. Approximate Bayes factors for genome-wide association studies: a new
   interpretation of p-values. *Genet Epidemiol* 2009;33(2):79–86. doi:10.1002/gepi.20207.
7. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful
   approach to multiple testing. *J R Stat Soc Series B* 1995;57(1):289–300.
   doi:10.1111/j.2517-6161.1995.tb02031.x.
8. Zhang Y, et al. SharePro: an accurate and efficient genetic colocalization method accounting
   for multiple causal variants. *Bioinformatics* 2024;40(5):btae295.
   doi:10.1093/bioinformatics/btae295.
9. Wingo TR, Damoiseaux R, Sittler A, et al. Integrating human brain proteomes with genome-wide
   association data implicates new proteins in Alzheimer disease. *Nat Genet* 2021;53(2):143–146.
   doi:10.1038/s41588-020-00773-z.
10. Vallabh S, Minikel EV, Prion Alliance. Open-science reanalysis and data-organization
    initiatives. Prion Alliance / Cure Alzheimer's Fund (cureffi.org).

---
*Figures: volcano_gse160208.png, volcano_gse140069.png, heatmap_top_genes.png
(pipeline/reports/figuras/). Fine-mapping numbers from relatorio_finemap_loci.md v2 and
relatorio_lambda_gc.md. Colocalization checked against `coloc.abf` 5.2.3; the other inferential
statistics are checked against SciPy and statsmodels in the test suite.*
