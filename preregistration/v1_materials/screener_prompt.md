# Screener prompts (plan B, Section 2)

Frozen in Zenodo v1. Screening has two steps with separate prompts. Both screeners receive exactly the text of each prompt (between the lines), followed by the input of that step. Nothing else is added (no gene names, no hints, no lists of expected results). Screener 1 runs in an environment without access to the authors' analysis files; screener 2 is a different AI model run by the author. The model name and version, the environment and the date are recorded for each run.

- **Step 1 (extraction)**: each screener receives the PubMed records and lists candidate pairs.
- Between the steps (rules in `README.md`): source identifiers are verified programmatically and unverifiable sources removed (listed in v2); gene names are normalised to HGNC symbols and Entrez IDs; the union of both screeners' candidates is formed with the merge key (Entrez ID, OncotreeLineage, source PMID, or DOI where there is no PMID), and the target lineage of each unit is fixed. The result is a fixed list of "candidate pair × source paper" units.
- **Step 2 (judgement)**: each screener receives the fixed list and judges every unit. Screeners judge the experimental facts; whether cell lines belong to the target lineage or to a different lineage is decided afterwards by script, with the lineage rule in `README.md`, and criteria 3, 4 and "eligible" are computed from the screeners' structured fields. Agreement is computed on step 2, with "candidate pair × source paper" as the unit, for c1, c2, c5, the computed c3 and c4, and "eligible": raw agreement, counts, the proportion "met" for each screener, unweighted Cohen's kappa (three categories) and Gwet's AC1 (kappa is depressed when almost all units are "not met"). Agreement of step 1 is reported as the overlap of the two candidate lists.

---

## Prompt for step 1 (extraction)

You are screening published papers for a literature benchmark. In this step you only list candidates; you do not judge them.

Rules:
- Do not use or look up any large-scale dependency resource or its results (for example DepMap, Project Achilles, Project DRIVE, Project Score, DEMETER, cancer dependency portals) or any table or manuscript about lineage-selective dependencies. Do not use your memory of such resources. Do not search for, open or use any repository, archive, preprint or preregistration of the study for which this screening is done (for example GitHub, Zenodo, OSF or bioRxiv records about lineage-selective dependencies in DepMap), or any document that lists expected results. Use web access only for PubMed/PMC, DOI resolution and publishers' full-text pages.
- Do not invent papers, PMIDs, DOIs, genes or cell lines. If you cannot open a paper, say so in the notes.
- List every candidate the records contain; do not try to produce a balanced or a short list.

Input: a batch of PubMed records (PMID, title, abstract; full text where available).
- For a primary research paper: list every gene for which the paper reports loss-of-function experiments in cancer cells (shRNA, siRNA, CRISPR knockout or interference, dominant-negative construct, or targeted degrader), with each cancer type in which this was done. Each gene × cancer type is one row.
- For a review: list every gene the review presents as a lineage dependency or lineage-survival gene in a specific cancer type, and, in separate rows, each primary paper the review cites for it (PMID or DOI as given in the review's reference list). Do not follow the references of those primary papers further.

Output: comma-separated values with this header, and nothing else except a final line "END":

record_PMID,record_type,source_PMID_or_DOI,gene_as_written,cancer_type_as_written,target_cell_lines,notes

record_type is "primary" or "review". For a primary paper, source_PMID_or_DOI is the paper itself. Put any text containing commas in double quotes; separate cell line names with semicolons.

Example of the format (invented, not a real paper):
12345678,primary,12345678,GENEX,"example carcinoma","LINE-A; LINE-B",

---

## Prompt for step 2 (judgement)

You are screening published papers for a literature benchmark. You receive a fixed list of units; each unit is one gene, one target lineage (with the cancer type(s) as written) and one source paper. For each unit you receive: unit_id, gene symbol (and the name as written), target lineage, cancer type(s) as written, and source PMID or DOI. Judge every unit in the list, from that source paper alone. You report what the paper shows; you do not decide which lineage a cell line belongs to (this is done afterwards from your list of cell lines). "Target cancer type" below means the cancer type(s) given for the unit.

Rules:
- Do not use or look up any large-scale dependency resource or its results (for example DepMap, Project Achilles, Project DRIVE, Project Score, DEMETER, cancer dependency portals) or any table or manuscript about lineage-selective dependencies. Do not use your memory of such resources. Do not search for, open or use any repository, archive, preprint or preregistration of the study for which this screening is done (for example GitHub, Zenodo, OSF or bioRxiv records about lineage-selective dependencies in DepMap), or any document that lists expected results. Use web access only for PubMed/PMC, DOI resolution and publishers' full-text pages. Judge only what the paper reports.
- Only experiments performed and reported in the paper itself count. Data reproduced from dependency databases or other published screens (for example a figure based on DepMap) do not count.
- Do not invent papers, PMIDs, DOIs, cell lines or results. If you cannot open the paper at all, mark every item "unclear" and say so.
- If you have the full text: information that is not reported → "not met"; information that is reported but ambiguous → "unclear".
- If you have only the abstract: items the abstract does not address → "unclear" (not assessable), not "not met". Mark the unit "abstract only".
- Judge each unit independently.

Items to judge:
1. c1 Perturbation: genetic loss of function (shRNA, siRNA, CRISPR knockout or interference, dominant-negative construct, or targeted degrader) of the gene was performed in cancer cells of the target cancer type. Drug inhibition alone does not count; record it in the notes.
2. c2 Reliability: at least 2 independent reagents (e.g. 2 different shRNAs or sgRNAs) or a rescue experiment.
3. Target cell lines: list every cancer cell line of the target cancer type in which the perturbation was tested, each with "met" if loss of function reduced proliferation or viability, "not met" if it did not, "unclear" otherwise.
4. Comparison cell lines: list every other cell line in which the same perturbation was tested, each with the cancer type or tissue the paper gives for it and with "met" if the effect was absent or at most half as large as in the target cell lines, "not met" if it was larger, "unclear" otherwise. Include non-cancer lines; they are excluded afterwards.
5. In vivo: whether loss of function reduced tumour growth in vivo in a model of the target cancer type (met / not met / unclear / not done).
6. c5 Origin of the gene: the paper's own stated reason for studying the gene, followed back one level of citation, is not a genome-scale loss-of-function screen (one targeting at least 1,000 genes; CRISPR, RNAi or similar, by any group). If the gene was first proposed by such a screen, c5 is "not met", even if later work confirmed it. Smaller, focused library screens do not make c5 "not met". A paper that cites a screen only as supporting evidence, after choosing the gene for another reason, can meet c5. If it cannot be determined whether the gene came from a genome-scale screen, c5 is "unclear".

Also record (these do not decide eligibility):
- Origin category, one of: DepMap-derived (DepMap, Project Achilles, Project DRIVE, Project Score, DEMETER2 component screens) / other genome-scale screen / biologically motivated candidate / genomics (amplification, mutation) / clinical / unclear. Here "unclear" means: not from a genome-scale screen, but the reason for choosing the gene fits none of the other categories.
- Subgroup restriction: whether the dependency is reported only in a molecular subgroup of the cancer type (yes / no; which subgroup).
- Affiliations of the corresponding and last authors, as printed in the paper.
- Publication year.
- Full text or abstract only.

Output: comma-separated values with this header, one row per unit in the order received, and nothing else except a final line "END":

unit_id,c1_perturbation,c2_reliability,target_lines_with_effect,comparison_lines_with_effect,in_vivo,c5_origin,origin_category,subgroup_restricted,subgroup,senior_author_affiliations,year,full_text_or_abstract,notes

Values: c1, c2 and c5 are "met", "not met" or "unclear". In target_lines_with_effect write each line as NAME:met, NAME:not met or NAME:unclear, separated by semicolons. In comparison_lines_with_effect write each line as NAME|cancer type or tissue as given:met (or not met, unclear), separated by semicolons. Put any field containing commas in double quotes.

Example of the format (invented, not a real paper):
U0001,met,met,"LINE-A:met; LINE-B:met","LINE-C|example sarcoma:met; LINE-D|example lymphoma:not met",not done,met,biologically motivated candidate,no,,"Example University",2010,full text,

---
