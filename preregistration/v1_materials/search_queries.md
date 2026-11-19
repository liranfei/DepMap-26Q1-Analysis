# PubMed search queries (plan B, Section 2)

Frozen in Zenodo v1. Each query is run once in PubMed (https://pubmed.ncbi.nlm.nih.gov/), sort = Best Match, and the first 100 records are taken (fewer if the query returns fewer). The search date, total hit count and the PMIDs of the records taken are recorded at the time of search. No date limit, no language limit.

Template: `(<cancer type terms>)` + common block:

```
(lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q01. Ampulla of Vater

```
("ampullary carcinoma"[tiab] OR "ampullary cancer"[tiab] OR "ampullary adenocarcinoma"[tiab] OR "ampulla of Vater"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q02. Biliary Tract

```
(cholangiocarcinoma[tiab] OR "biliary tract cancer"[tiab] OR "biliary tract carcinoma"[tiab] OR "bile duct cancer"[tiab] OR "gallbladder cancer"[tiab] OR "gallbladder carcinoma"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q03. Bladder/Urinary Tract

```
("bladder cancer"[tiab] OR "bladder carcinoma"[tiab] OR "urothelial carcinoma"[tiab] OR "urothelial cancer"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q04. Bone

```
("Ewing sarcoma"[tiab] OR "Ewing's sarcoma"[tiab] OR osteosarcoma[tiab] OR chondrosarcoma[tiab] OR chordoma[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q05. Bowel

```
("colorectal cancer"[tiab] OR "colorectal carcinoma"[tiab] OR "colorectal adenocarcinoma"[tiab] OR "colon cancer"[tiab] OR "colon carcinoma"[tiab] OR "rectal cancer"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q06. Breast

```
("breast cancer"[tiab] OR "breast carcinoma"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q07. CNS/Brain

```
(glioma[tiab] OR glioblastoma[tiab] OR astrocytoma[tiab] OR medulloblastoma[tiab] OR "atypical teratoid"[tiab] OR meningioma[tiab] OR oligodendroglioma[tiab] OR ependymoma[tiab] OR "embryonal tumor"[tiab] OR "embryonal tumour"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q08. Cervix

```
("cervical cancer"[tiab] OR "cervical carcinoma"[tiab] OR "cervical squamous cell carcinoma"[tiab] OR "cervical adenocarcinoma"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q09. Esophagus/Stomach

```
("esophageal cancer"[tiab] OR "oesophageal cancer"[tiab] OR "esophageal carcinoma"[tiab] OR "oesophageal carcinoma"[tiab] OR "esophageal squamous cell carcinoma"[tiab] OR "esophageal adenocarcinoma"[tiab] OR "gastric cancer"[tiab] OR "gastric carcinoma"[tiab] OR "gastric adenocarcinoma"[tiab] OR "stomach cancer"[tiab] OR esophagogastric[tiab] OR "gastroesophageal junction"[tiab] OR "gastro-oesophageal junction"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q10. Eye

```
("uveal melanoma"[tiab] OR "ocular melanoma"[tiab] OR retinoblastoma[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q11. Head and Neck

```
("head and neck squamous cell carcinoma"[tiab] OR "head and neck cancer"[tiab] OR HNSCC[tiab] OR "oral squamous cell carcinoma"[tiab] OR "laryngeal carcinoma"[tiab] OR "pharyngeal carcinoma"[tiab] OR "salivary gland carcinoma"[tiab] OR "nasopharyngeal carcinoma"[tiab] OR "nasopharyngeal cancer"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q12. Kidney

```
("renal cell carcinoma"[tiab] OR "kidney cancer"[tiab] OR "renal cancer"[tiab] OR "rhabdoid tumor of the kidney"[tiab] OR "rhabdoid tumour of the kidney"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q13. Liver

```
("hepatocellular carcinoma"[tiab] OR "liver cancer"[tiab] OR hepatoblastoma[tiab] OR hepatoma[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q14. Lung

```
("lung cancer"[tiab] OR "lung adenocarcinoma"[tiab] OR "lung squamous cell carcinoma"[tiab] OR "non-small cell lung"[tiab] OR NSCLC[tiab] OR "small cell lung cancer"[tiab] OR SCLC[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q15. Lymphoid

```
(lymphoma[tiab] OR "lymphoblastic leukemia"[tiab] OR "lymphoblastic leukaemia"[tiab] OR "lymphocytic leukemia"[tiab] OR "lymphocytic leukaemia"[tiab] OR "multiple myeloma"[tiab] OR myeloma[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q16. Myeloid

```
("acute myeloid leukemia"[tiab] OR "acute myeloid leukaemia"[tiab] OR AML[tiab] OR "chronic myeloid leukemia"[tiab] OR "chronic myeloid leukaemia"[tiab] OR "myeloproliferative neoplasm"[tiab] OR "acute myelogenous leukemia"[tiab] OR "acute myelogenous leukaemia"[tiab] OR "myeloproliferative neoplasms"[tiab] OR myelodysplastic[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q17. Ovary/Fallopian Tube

```
("ovarian cancer"[tiab] OR "ovarian carcinoma"[tiab] OR "high-grade serous"[tiab] OR "fallopian tube carcinoma"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q18. Pancreas

```
("pancreatic cancer"[tiab] OR "pancreatic ductal adenocarcinoma"[tiab] OR "pancreatic adenocarcinoma"[tiab] OR "pancreatic carcinoma"[tiab] OR PDAC[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q19. Peripheral Nervous System

```
(neuroblastoma[tiab] OR "malignant peripheral nerve sheath tumor"[tiab] OR "malignant peripheral nerve sheath tumour"[tiab] OR MPNST[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q20. Pleura

```
(mesothelioma[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q21. Prostate

```
("prostate cancer"[tiab] OR "prostate carcinoma"[tiab] OR "prostate adenocarcinoma"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q22. Skin

```
(melanoma[tiab] OR "cutaneous squamous cell carcinoma"[tiab] OR "Merkel cell carcinoma"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q23. Soft Tissue

```
(rhabdomyosarcoma[tiab] OR liposarcoma[tiab] OR "synovial sarcoma"[tiab] OR leiomyosarcoma[tiab] OR "soft tissue sarcoma"[tiab] OR "undifferentiated pleomorphic sarcoma"[tiab] OR "rhabdoid tumor"[tiab] OR "rhabdoid tumour"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q24. Testis

```
("testicular germ cell tumor"[tiab] OR "testicular germ cell tumour"[tiab] OR "testicular cancer"[tiab] OR "embryonal carcinoma"[tiab] OR nonseminoma[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q25. Thyroid

```
("thyroid cancer"[tiab] OR "thyroid carcinoma"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q26. Uterus

```
("endometrial cancer"[tiab] OR "endometrial carcinoma"[tiab] OR "uterine cancer"[tiab] OR "uterine carcinoma"[tiab] OR "endometrial adenocarcinoma"[tiab]) AND (lineage[tiab] OR "master regulator"[tiab] OR "lineage-survival"[tiab] OR "lineage-specific"[tiab] OR "lineage dependency"[tiab]) AND (shRNA[tiab] OR siRNA[tiab] OR CRISPR[tiab] OR knockdown[tiab] OR knockout[tiab]) AND (proliferation[tiab] OR viability[tiab] OR survival[tiab] OR dependency[tiab])
```

## Q27. Reviews on lineage-survival oncogenes and lineage transcription-factor dependencies

```
("lineage-survival oncogene"[tiab] OR "lineage-survival oncogenes"[tiab] OR "lineage survival oncogene"[tiab] OR "lineage survival oncogenes"[tiab] OR "lineage dependency"[tiab] OR "lineage dependencies"[tiab] OR "lineage-specific dependency"[tiab] OR "lineage-specific dependencies"[tiab] OR "lineage transcription factor"[tiab] OR "lineage transcription factors"[tiab] OR "master transcription factor"[tiab] OR "master transcription factors"[tiab]) AND (cancer[tiab] OR tumor[tiab] OR tumour[tiab] OR carcinoma[tiab] OR leukemia[tiab] OR leukaemia[tiab] OR lymphoma[tiab] OR sarcoma[tiab]) AND review[pt]
```

From each review taken, every gene presented as a lineage dependency or lineage-survival oncogene in a specific cancer type is listed, and the primary papers the review cites for it are screened against the criteria (one level of snowballing only: references of those primary papers are not followed).
## Pre-freeze syntax check

On 2026-11-19, before freezing, each query was sent to the PubMed E-utilities `esearch` service with `rettype=count` only, to check syntax and that no query returns zero records; no record, title or PMID was retrieved or viewed. Counts on that date: Q01 2, Q02 5, Q03 6, Q04 20, Q05 39, Q06 80, Q07 57, Q08 6, Q09 31, Q10 19, Q11 17, Q12 13, Q13 43, Q14 68, Q15 100, Q16 85, Q17 26, Q18 30, Q19 25, Q20 3, Q21 57, Q22 51, Q23 3, Q24 2, Q25 6, Q26 4, Q27 113. After this first check, synonyms were added to Q07, Q08, Q09, Q11, Q13, Q16 and Q26 following an audit of the queries by an AI model in a separate session (no record had been viewed), and the count check was repeated (below).

Second count check, 2026-11-19, after adding the synonyms (counts only, no record viewed): Q01 2, Q02 5, Q03 6, Q04 20, Q05 39, Q06 80, Q07 59, Q08 6, Q09 31, Q10 19, Q11 20, Q12 13, Q13 44, Q14 68, Q15 100, Q16 95, Q17 26, Q18 30, Q19 25, Q20 3, Q21 57, Q22 51, Q23 3, Q24 2, Q25 6, Q26 4, Q27 113. The queries were not changed after this second check. Counts at the actual search date may differ and are recorded then, together with the full list of PMIDs returned by each query (esearch with retmax set to the hit count), so that the effect of the Best Match cut-off can be checked later; Best Match ordering is not reproducible over time.
