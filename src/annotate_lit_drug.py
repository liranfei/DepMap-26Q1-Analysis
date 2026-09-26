"""Coarse external annotation of the candidate pairs: (i) PubMed record counts for gene + lineage-specific cancer terms + perturbation/dependency terms; (ii) DGIdb drug-gene interactions and druggable-genome categories.
Both are automated proxies (record counts do not establish that a dependency was reported); queried 2026-09-27."""
import json, time, re, urllib.parse, urllib.request, pandas as pd, os
O = os.environ.get("RESULTS_DIR", "results") + "/"
F = pd.read_csv(O + "final_targets.csv"); F["Symbol"] = F.Gene.str.replace(r" \(.*", "", regex=True)
LT = {"ampulla of vater": "ampullary OR \"ampulla of Vater\"", "bladder/urinary tract": "bladder OR urothelial", "bowel": "colorectal OR colon OR \"colon cancer\"", "breast": "\"breast cancer\"", "cervix": "cervical cancer",
      "cns/brain": "glioma OR glioblastoma OR medulloblastoma OR neuroblastoma OR astrocytoma", "esophagus/stomach": "esophageal OR gastric", "eye": "\"uveal melanoma\" OR retinoblastoma OR \"ocular melanoma\"",
      "head and neck": "\"head and neck\"", "kidney": "\"renal cell\" OR \"kidney cancer\" OR \"rhabdoid\"", "liver": "hepatocellular OR \"liver cancer\"", "lymphoid": "lymphoma OR \"lymphoblastic leukemia\" OR myeloma OR \"lymphoid\"",
      "myeloid": "\"myeloid leukemia\" OR AML OR \"myeloid\"", "ovary/fallopian tube": "ovarian", "pancreas": "pancreatic", "peripheral nervous system": "neuroblastoma OR \"nervous system\"", "pleura": "mesothelioma",
      "prostate": "\"prostate cancer\"", "skin": "melanoma OR \"skin cancer\"", "soft tissue": "sarcoma OR rhabdomyosarcoma", "thyroid": "thyroid cancer"}
PERT = "(CRISPR OR knockout OR knockdown OR siRNA OR shRNA OR dependency OR essential OR vulnerability OR inhibitor)"
def pm(term):
    u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term=" + urllib.parse.quote(term)
    for _ in range(4):
        try: return int(json.load(urllib.request.urlopen(u, timeout=30))["esearchresult"]["count"])
        except Exception: time.sleep(2)
    return None
rows = []
for _, r in F.iterrows():
    q = f"({r.Symbol}[Title/Abstract]) AND ({LT[r.Lineage]}) AND {PERT}"; n = pm(q); time.sleep(0.4)
    rows.append((r.Symbol, r.Lineage, n))
L = pd.DataFrame(rows, columns=["Symbol", "Lineage", "pubmed_records"]); L["pubmed_category"] = pd.cut(L.pubmed_records, [-1, 0, 4, 10**9], labels=["none", "1-4", ">=5"])
syms = sorted(F.Symbol.unique()); D = {}
def dg(names):
    q = "{ genes(names: %s) { nodes { name geneCategories { name } interactions { drug { name approved } interactionScore } } } }" % json.dumps(names)
    req = urllib.request.Request("https://dgidb.org/api/graphql", data=json.dumps({"query": q}).encode(), headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=60))["data"]["genes"]["nodes"]
for i in range(0, len(syms), 20):
    for nd in dg(syms[i:i + 20]):
        ap = sorted({x["drug"]["name"] for x in nd["interactions"] if x["drug"]["approved"]}); allr = nd["interactions"]
        top = [x["drug"]["name"] for x in sorted(allr, key=lambda x: -(x["interactionScore"] or 0))[:3]]
        D[nd["name"]] = dict(n_interactions=len(allr), n_approved_drugs=len(ap), top_drugs=top, categories=sorted({c["name"] for c in nd["geneCategories"]}))
G = pd.DataFrame([dict(Symbol=s, **D.get(s, dict(n_interactions=0, n_approved_drugs=0, top_drugs=[], categories=[]))) for s in syms])
G["druggable_genome"] = G.categories.apply(lambda c: any("DRUGGABLE" in x.upper() or "KINASE" in x.upper() or "CLINICALLY" in x.upper() for x in c))
L.to_csv(O + "annotation_pubmed.csv", index=False); G.to_csv(O + "annotation_dgidb.csv", index=False)
S = dict(n_pairs=len(L), pubmed_none=int((L.pubmed_category == "none").sum()), pubmed_1_4=int((L.pubmed_category == "1-4").sum()), pubmed_ge5=int((L.pubmed_category == ">=5").sum()), pubmed_failed=int(L.pubmed_records.isna().sum()),
         n_genes=len(G), genes_with_any_interaction=int((G.n_interactions > 0).sum()), genes_with_approved_drug=int((G.n_approved_drugs > 0).sum()), genes_druggable_category=int(G.druggable_genome.sum()),
         approved_drug_genes=G[G.n_approved_drugs > 0].Symbol.tolist(), no_interaction_genes=G[G.n_interactions == 0].Symbol.tolist(), no_pubmed_pairs=L[L.pubmed_category == "none"][["Symbol", "Lineage"]].values.tolist())
json.dump(S, open(O + "annotation_summary.json", "w"), indent=1); print(json.dumps(S, indent=1))
