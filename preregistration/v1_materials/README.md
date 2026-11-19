# Zenodo v1 materials (plans A and B)

Prepared 2026-11-19. Deposited together with `plan_A_stratified_procedure.md` and `plan_B_positive_controls.md` as Zenodo v1, before screening starts and before any simulation or real-data analysis of plan A is run. SHA-256 of every file in the deposit is listed in `SHA256SUMS.txt` at the top level of the deposit, generated at deposit.

| File | Content | Used by |
|---|---|---|
| `screener_prompt.md` | Prompts for step 1 (extraction) and step 2 (judgement), given to both AI screeners | plan B §2 |
| `search_queries.md` | 26 lineage queries (Q01–Q26) and the review query (Q27); two pre-freeze count checks | plan B §2 |
| `lineage_mapping.csv` | Cancer types as named in papers → OncotreeLineage, one row per case | plan B §1 criteria 1, 3, 4 and 6 |
| `screen_project_institutions.csv` | Institutions of the genome-scale screen projects, for the screen-project flag | plan B §1 |
| `requirements_lock.txt` | `pip freeze` of the analysis environment (Python 3.13.12) | plan A §7 |
| `literature_search_novelty.md` | Documented search for prior culture-stratified lineage-dependency procedures (2026-11-19), with permitted wording | plan A §6 |
| `literature_search_records/` | Query strings and all retrieved records (JSON) of that search | plan A §6 |
| `../audits/` | Seven audit reports on the drafts and materials (in Chinese), produced by an AI model (Claude, Anthropic) in a session separate from the one that drafted the plans; they are not independent human review. The final drafts show how each point was handled | plans A and B |
| `README.md` | This file, including the rules below | plan B §1–2 |

## Release tags

"Zenodo v1" is the repository release `v2.5.0-prereg`; v1b, v2 and v3 will be `v2.5.1-prereg-v1b`, `v2.5.2-prereg-v2` and `v2.5.3-prereg-v3` (plan A, Section 7). Each has its own Zenodo version DOI within the repository's concept record. The deposit is dated 2026-11-19; if the release had been made on a later day, the dates in the status lines would have been changed before tagging.

## Lineage-assignment rule (plan B, criteria 1, 3, 4 and 6)

The same rule is used for the target cell lines and for the comparison cell lines. Rules 1–3 assign a lineage to each cell line: between the steps for the cell lines listed in step 1 (to fix the target lineage of each unit), and by script after step 2 for the cell lines each screener listed in step 2 (to compute criteria 1, 3 and 4 against the fixed target lineage). Rule 4 is used only between the steps.

1. Each named cell line is looked up in DepMap 26Q1 `Model.csv` (annotations only, no screening results) by exact match after converting to upper case and removing all characters other than letters and digits, against `StrippedCellLineName`, then `CellLineName`, then the part of `CCLEName` before the first underscore. A line with no match, or with matches to more than one model of different lineages, is "unmatched". The lineage of a matched line is its `OncotreeLineage`.
2. A line whose matched model has `OncotreePrimaryDisease` = "Non-Cancerous" is a non-cancer line: it does not count as a target or a comparison cancer line.
3. An unmatched line is assigned by the cancer type or tissue the paper gives for it, using `lineage_mapping.csv`; a line given as non-malignant or immortalised is a non-cancer line.
4. Target lineage of a unit (used only between the steps, rule 3 of the next section): the lineage holding most of the matched target cancer lines listed in step 1 by either screener (unmatched lines are used only if no target line is matched); a tie makes the unit "unclear". If no target line is named, or the evidence is in vivo only, the target cancer type given in step 1 is mapped with `lineage_mapping.csv`; such units cannot meet criterion 1 unless step 2 finds target cell lines.
5. Cancer types not in the table, or marked "not mapped", make the unit "unclear" (listed, not analysed); a comparison line of such a type does not count for criterion 4.
6. Only the 26 lineages with ≥ 5 cancer models in DepMap 26Q1 can receive a positive control; pairs mapping to other lineages are listed, not analysed.
7. Matching of cancer-type text to rows of `lineage_mapping.csv` is done by the author; doubtful cases are recorded with the reason.

## Rules between screening steps 1 and 2

1. **Identifier verification**: every source PMID or DOI is checked programmatically against PubMed (and Crossref for DOIs) for title and year. DOIs are converted to PMIDs where PubMed has the record; otherwise the DOI is kept. Sources that cannot be verified are removed and listed in v2.
2. **Gene normalisation**: gene names as written are matched, after converting to upper case and removing spaces, against the HGNC complete set (downloaded on the day of normalisation; its date and SHA-256 are recorded in v2) in this order: approved symbol; previous symbols; alias symbols. A name is normalised only if it matches exactly one approved gene at the first level where it matches; otherwise the author assigns it by reading the paper and records the reason. The approved symbol and its Entrez (NCBI Gene) ID are used.
3. **Provisional lineage for merging**: the cell lines that either screener listed for the unit in step 1 (union of both lists) are assigned by the lineage rule above; where the screeners' lists conflict in a way that changes the lineage, the author reads the paper and records the decision.
4. **Merge key**: (Entrez ID, OncotreeLineage, source PMID, or DOI where there is no PMID). Several cancer types of the same lineage in one paper form one unit. The union of both screeners' step-1 candidates under this key is the fixed list for step 2; the target lineage of each unit is fixed here and not re-determined later.
5. **Step-2 input** for each unit: `unit_id`, approved gene symbol (with the name as written), target lineage, the cancer type(s) as written in step 1, and the source PMID or DOI. The cell lines listed in step 1 are not given, so that step 2 is judged afresh from the paper.

## Screen-project flag

An eligible paper is flagged if the printed affiliation of its corresponding or last author contains, case-insensitively, any of the terms in the `match_terms` column of `screen_project_institutions.csv` (terms separated by "; "). This is a coarse flag used only for a sensitivity analysis and is applied by the author, not by the screeners. Several of the institutions are large, so many papers unrelated to screen projects will be flagged; for a sensitivity analysis that removes flagged papers this errs on the conservative side.

## Environment note

`requirements_lock.txt` is a `pip freeze` of the global Python environment used for the analysis; it contains packages that are not used. Only the packages in the repository's `requirements.txt` (numpy, pandas, scipy, statsmodels, lifelines, matplotlib, adjustText, openpyxl, requests, Pillow) are imported, at the versions shown in the lock file.

## What was done before freezing

- The 26 lineages and their primary-disease composition were read from `Model.csv` (annotation only).
- The queries were checked for syntax with PubMed `esearch` counts only, twice (see `search_queries.md`); no record was viewed.
- Garraway & Sellers 2006, Table 2 was checked against the publisher's page (doi:10.1038/nrc1947, table 2) on 2026-11-19.
- The plans and materials were revised after seven rounds of audit by an AI model (Claude) in a session separate from the one that drafted the plans (not independent human review); all seven reports are deposited in `audits/`. Two earlier rounds of review, also by Claude in another separate session, were given in that session and are not deposited as files.
- A documented literature search on prior procedures was run (`literature_search_novelty.md`); it concerned methods only, and no gene was recorded from it.
- No positive-control candidate was looked up in any result file.
