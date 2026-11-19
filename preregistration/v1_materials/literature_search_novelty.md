# Literature search on prior procedures (plan A §6)

Run on 2026-11-19, before Zenodo v1, by the author with AI assistance (Claude). Purpose: decide whether the manuscript may say that no comparable procedure was found (plan A §6). This search is about **methods only**. It is not part of the positive-control search of plan B, no gene was recorded from it, and no result file of this study was opened.

Raw outputs are in `literature_search_records/`: the query strings (`queries.json`, `queries_extra.json`, `queries_supplement.json`), every retrieved PubMed and Europe PMC record (PMID, year, journal, title) and the result lists of the three web searches (`web_searches_2026-11-19.json`), as JSON. A supplement (L8, L9, W3, the tool documentation and the methods of the consortium papers) was added on the same day, after a pre-release audit (`audits/独立核对意见_冻结前复查.md`).

## Question

Is there a published procedure for **genome-wide lineage-dependency testing in cancer cell-line loss-of-function screens** that (i) stratifies or adjusts for **culture format** (suspension / adherent), and (ii) evaluates its **type I error and power** (e.g. by simulation)?

## Databases and queries

PubMed (NCBI E-utilities `esearch`, title/abstract fields, no date or language limits). Full query strings are in `literature_search_records/`.

| ID | Topic | PubMed hits |
|---|---|---|
| L1 | DepMap / CRISPR screen / gene dependency × lineage / tissue- / context-specific × stratified / covariate / confounding | 9 |
| L2 | DepMap / CRISPR screen / gene dependency × culture condition / growth pattern / suspension / culture format / medium | 19 |
| L3 | Van Elteren / stratified Wilcoxon / stratified rank × gene / genomic / CRISPR / expression / cancer | 8 |
| L4 | lineage-selective / lineage-specific / selective dependency × method / statistical / framework / benchmark | 29 |
| L5 | DepMap / CRISPR screen / gene effect × confounding / bias × lineage / cancer type / tissue × cell line | 1 |
| L6 | blood / haematopoietic / leukaemia × suspension / adherent × DepMap / CRISPR screen / gene dependency | 0 |
| L7 | Hodges–Lehmann / stratified Wilcoxon / Van Elteren / Cochran–Mantel–Haenszel × CRISPR / DepMap / gene dependency / functional genomic | 0 |
| L8 | DepMap / CRISPR screen / gene effect / gene dependency × limma / linear model / linear regression / mixed model × lineage / tissue / cancer type | 2 |
| L9 | DepMap / CRISPR screen / gene effect / gene dependency × growth rate / doubling time / proliferation rate × confounding / covariate / adjustment / bias | 0 |

Europe PMC (abstract field, includes preprints):
- E1: DepMap or CRISPR screens × culture condition / growth pattern / suspension × lineage / confound: 0 hits.
- E2: Van Elteren / stratified Wilcoxon × CRISPR / dependency / cell lines: 0 hits.

Europe PMC all-field counts for L1–L4 (2,101; 4,827; 468; 1,011) were recorded as a check only. All-field matching includes full text, so these sets were too unspecific to screen and were not screened.

Three general web searches (2026-11-19) were also run, to catch preprints and portal documentation; the result lists are archived:
- W1: "DepMap lineage-selective dependencies confounding culture condition suspension adherent growth pattern statistical adjustment"
- W2: "bioRxiv 2025 2026 CRISPR screen cancer cell lines 'culture conditions' confound gene essentiality tissue-specific dependencies stratified analysis"
- W3: "DepMap portal Context Explorer methods in-context out-of-context t-test comparison group selective dependency documentation"

Tool documentation and consortium papers were checked directly, because their lineage tests are described in the methods or in the documentation rather than in titles or abstracts (see "Tools and consortium methods" below).

## Screening

All 68 PubMed records of L1–L9 (66 distinct; PMID 31039782 appears in L1, L2 and L5) were screened by title. Records that might be relevant were also screened by abstract or full text. No record met both (i) and (ii). The closest records are listed below (one of them, Shi 2024, was found by web search).

| Reference | What it does | (i) culture format | (ii) error / power |
|---|---|---|---|
| Liang S et al. 2021, IEEE/ACM TCBB, doi:10.1109/tcbb.2021.3094650 (PMID 34232885) | Van Elteren test for differential expression in single-cell data, stratified by batch / participant; simulation and ROC evaluation | No: strata are batches, not culture format; the setting is not dependency screens | Yes |
| Lagziel S et al. 2019, BMC Biol, doi:10.1186/s12915-019-0654-4 (PMID 31039782) | Shows culture medium confounds cell-line-specific essentiality of metabolic genes; predictive models include medium | Medium, not format; predictive models for metabolic genes, not lineage tests | No |
| Neiswender JV et al. 2026, Nature 657:789–799, doi:10.1038/s41586-026-10843-7 (PMID 42557315) | 3D / NextGen DepMap. Linear regression in a 250-gene minipool screen separates growth format (dome vs plastic) from medium. Reports that growth format had little effect under serum-free conditions, unlike an earlier serum-containing study | 2D vs 3D, not suspension vs adherent; not genome-wide lineage testing | No |
| Shi X et al. 2024, Nat Cancer, doi:10.1038/s43018-024-00789-y (PMID 39009815; found by web search W1, not by L1–L9; checked in the full text, PMC11358024) | TCGA translational dependency map. Prediction models are compared with models built on the DepMap confounders dataset (sex, Cas9 activity, age, lineage, primary or metastasis, growth pattern, library, screen quality, cancer type) as a null distribution | Covariate in a null model for prediction, not a lineage test | No |
| Rossiter NJ et al. 2021, Cell Metab, doi:10.1016/j.cmet.2021.02.005 (PMID 33651980) | CRISPR screens in a medium resembling human plasma; conditionally essential genes in human cancer cell lines (abstract wording) | Medium, experimental | No |
| Metz P et al. 2026, BMC Genomics, doi:10.1186/s12864-026-12658-2 (PMID 41715004) | Guide library representation biases gene effect scores | Library, not format (relevant to the plan A library sensitivity strata) | No |

## Tools and consortium methods

Checked on 2026-11-19 in the full text (PMC, or the publisher's page) or the documentation, by searching the text for the terms suspension, adherent, culture condition, culture medium, growth pattern, t-test, linear model, limma, Wilcoxon, Mann–Whitney, covariate, doubling and growth rate, and reading the sections found.

| Source | Lineage / context test | Culture format |
|---|---|---|
| Meyers RM et al. 2017, Nat Genet (CERES), PMC5709193 | Lineage vs rest of the dataset: difference in means, two-tailed Student's t-test, FDR | Culture medium is listed per line in a table; not used in the test |
| Tsherniak A et al. 2017, Cell (DEMETER), PMC5667678 | Culture medium and conditions are recorded in the sample sheet; no culture adjustment of lineage tests found | Recorded only |
| Behan FM et al. 2019, Nature (Project Score), doi:10.1038/s41586-019-1103-9 (publisher's page) | Cancer-type-specific core fitness genes; ANOVA for biomarker associations within cancer types | None of the search terms for culture found |
| Pacini C et al. 2021, Nat Commun (integrated Broad–Sanger screens), PMC7955067 | — | None of the search terms found |
| Dharia NV et al. 2021, Nat Genet (paediatric cancers), PMC8049517 | Tumour type vs all other lines: limma `lmFit`/`eBayes` and two-sided t-tests with BH correction; Fisher's exact test on dependency probability | Not adjusted |
| Dempster JM et al. 2021, Genome Biol (Chronos), PMC8686573 | Not a lineage test; Chronos models cell growth rate and screen quality per line when estimating gene effects | Growth rate, not culture format |
| DepMap portal "enriched in" lineage / Context Explorer (DepMap forum, staff answer, 2021; Context Explorer introduced in 25Q3) | For each gene and context, in-group vs out-group t-test; p values not corrected for multiple testing | Not adjusted. Whether the out-group can be chosen (for example other blood or other solid lines) could not be checked: the portal page returned a browser verification page |
| gdsctools ANOVA (GDSC, documentation `anova_parttwo`) | Drug response ~ tissue + medium + MSI + feature, sequential (type I) sums of squares; medium and MSI included by default | Culture **medium** as covariate; drug response, not lineage dependency; no simulation of error rates reported |

Lineage-Selective Dependencies in Pediatric Cancers (Cold Spring Harb Perspect Med 2025, doi:10.1101/cshperspect.a041573) was found in the web search. It is a review and was not counted as a procedure. Its full text could not be opened (browser verification page), so its content was not checked.

## Conclusion and permitted wording

- **What was found**:
  - A stratified Wilcoxon (Van Elteren) test with simulation evaluation exists in genomics, for batch effects in single-cell expression data (Liang 2021).
  - Culture medium and growth format have been handled as covariates or regression terms in other settings (Lagziel 2019, Neiswender 2026, Shi 2024; culture medium in the GDSC ANOVA).
  - The lineage tests of the DepMap and Project Score consortium papers and of the DepMap portal compare a lineage with all other lines (t-test, limma or ANOVA) without stratifying by culture format.
  - No published procedure was found that does both (i) and (ii) for genome-wide lineage dependencies.
- **Wording the manuscript may use** (plan A §6): "In a documented literature search (2026-11-19; queries archived with the preregistration), we found no published genome-wide procedure for lineage-dependency testing in cancer cell-line screens that stratifies by culture format and reports type I error and power; the Van Elteren test itself is standard and has been applied to batch-stratified differential expression in single-cell data [Liang 2021]. Culture conditions have been handled as covariates in pharmacogenomic association testing (culture medium in the GDSC ANOVA) and in prediction models [Lagziel 2019; Shi 2024], and lineage tests in dependency studies compare a lineage with all other lines [Meyers 2017; Dharia 2021]; none of these stratifies lineage tests by culture format or evaluates error rates by simulation."
- **Not permitted**: "first", "novel test", or any statement that no such procedure exists.
- **The manuscript must cite** Liang 2021 (the closest methodological precedent), Lagziel 2019, Neiswender 2026, the GDSC ANOVA (gdsctools), Meyers 2017 and Dharia 2021.
  - Neiswender 2026 reports a small growth-format effect under serum-free conditions. This is relevant to the interpretation of plan A; the DepMap 2D lines are mostly in serum-containing medium.

## Limitations

- Two databases (PubMed, Europe PMC), three web searches, six consortium papers and two tools. Google Scholar, Web of Science, Scopus and CRAN/Bioconductor package documentation were not searched.
- The consortium papers were checked by searching their text for fixed terms and reading the sections found, not by reading every supplementary file; Behan 2019 was read on the publisher's page without its supplementary methods.
- Neiswender 2026 and Shi 2024 appeared after, or near, the knowledge cutoff of the AI model that wrote the pre-release audit; their descriptions here are from the publisher's page (Neiswender 2026) and the full text in PMC (Shi 2024), and were not checked by that audit.
- Screening was done by one reviewer (AI-assisted, then checked by the author), not in duplicate.
- Methods described only in the supplementary material of a biological paper may be missed by title/abstract queries.
- This is a documented scoping search, not a systematic review. The wording above is limited accordingly ("we found no").
