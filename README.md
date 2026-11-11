# Pathway-enrichment-comparison

Comparative analysis of pathway enrichment methods (ORA, GSEA, CAMERA) and gene-set resources (GO Biological Process, Reactome, MSigDB Hallmark, C7 Immunologic Signatures, WikiPathways) using a controlled LPS stimulation experiment in primary human macrophages (GSE193336).

The goal is not to identify the "best" method, because with 4 paired donors no single dataset can support that conclusion. Instead, this is a well-controlled case study of comparative method behavior: how many pathways each method identifies, how redundant those results are, how well each method recovers biologically expected programs, and how much methods agree with each other across the same gene-set resource.

## Conceptual framework

Pathway enrichment analysis involves three independent choices that are frequently conflated in the literature. Separating them makes method comparison cleaner and results more interpretable.

| Layer | Options | Implemented here |
|---|---|---|
| Statistical method | ORA, GSEA, CAMERA, GSVA, ssGSEA, ROAST/fry | ORA, GSEA, CAMERA |
| Gene set resource | GO BP, Reactome, Hallmark, KEGG, WikiPathways, C7 Immunologic Signatures | GO BP, Reactome, Hallmark, C7, WikiPathways |
| Software implementation | clusterProfiler, fgsea, ReactomePA, GSVA, gprofiler2, msigdbr | clusterProfiler, fgsea, ReactomePA, msigdbr |

These layers are orthogonal, so any method can be applied to any gene set resource using any compatible software. A common source of confusion is treating combinations like "GSEA with Hallmark" as a single unit rather than independent method × resource choices.

## When to use each method

| Scenario | Recommended method | Reason |
|---|---|---|
| Strong DE signal, many significant genes, fast exploratory analysis | ORA | Simple, interpretable; works well when threshold is easy to justify |
| Moderate signal distributed across many genes, no clear threshold | GSEA | Uses full ranked list; detects coordinated shifts missed by hard cutoffs |
| Few significant genes but biological signal expected | GSEA | Captures modest coordinated shifts; ORA would return nothing |
| Paired or blocked design (repeated measures, matched samples) | CAMERA | Accounts for inter-gene correlation within gene sets; ORA/GSEA treat genes as independent |
| Small sample size (n < 10 per group) | CAMERA | Correlation correction most important when sample size is small |
| Complex limma model (batch + covariates + condition) | CAMERA | Operates directly on voom object; enrichment test is statistically consistent with DE model |
| Comparing pathway activity across individual samples | GSVA / ssGSEA | Produces per-sample scores; ORA/GSEA give group-level results only |
| Hypothesis-free discovery across broad biology | GO BP + GSEA | Broad coverage; full ranked list avoids threshold sensitivity |
| Clean non-redundant summary for presentation | Hallmark + ORA or GSEA | 50 curated sets; interpretable without post-processing |
| Infection or immune cell biology | C7 Immunologic Signatures | Provides experiment-level resolution; top terms name matched published immune perturbations. High library size (~5,219 sets) means more hits but also more potential false positives, so use alongside curated resources |

**Note on GSE193336:** The paired donor design (4 donors, both conditions per donor) and small sample size make CAMERA the most statistically principled choice for this dataset. ORA and GSEA do not account for donor pairing at the pathway level, whereas CAMERA does, because it operates on the voom object with the full `~ donor + condition` design. Including all three methods allows direct comparison of how much the correlation correction changes results in practice.

## Dataset

| Accession | Description | Design |
|---|---|---|
| GSE193336 | Bulk RNA-seq, primary human macrophages | Monocyte-derived macrophages (M-CSF, 5 days): Unstimulated vs LPS 50 ng/ml, 24h, with 4 paired donors|

**Why this dataset:** LPS stimulation of primary human macrophages is one of the most well-characterized transcriptional perturbations in immunology. The expected biology (TNF/NF-κB signaling, interferon response, cytokine production, complement activation) is well-defined, making it possible to evaluate whether each method recovers the programs it should find. The paired donor design controls for inter-individual variation, and the 24h timepoint captures the fully developed LPS transcriptional response.

**Limitation:** 4 donors is a small sample size. Results should be interpreted as illustrative of method behavior in a well-controlled setting, not as a definitive ranking of methods.

## Repository structure

```
Pathway-enrichment-comparison/
├── scripts/
│   ├── 00_download_GSE193336.R          # Download + process GSE193336 count matrix
│   ├── 01_differential_expression.R     # limma-voom DE (LPS vs Unstimulated, paired donors)
│   ├── 02_ora_gsea.R                    # ORA + GSEA across GO BP, Reactome, Hallmark
│   ├── 03_camera.R                      # CAMERA across GO BP, Reactome, Hallmark
│   ├── 04_go_redundancy.R               # GO BP redundancy quantification (GOSemSim/rrvgo)
│   └── 05_additional_resources.R        # ORA + GSEA + CAMERA: C7 + WikiPathways
├── data/
│   ├── raw/                             # Downloaded files (not tracked)
│   └── processed/                       # Filtered count matrix, metadata (not tracked)
├── results/                             # DE results, enrichment outputs (not tracked)
└── figures/
    ├── benchmark_GO_ORA_vs_GSEA_vs_CAMERA.png
    ├── benchmark_Reactome_ORA_vs_GSEA_vs_CAMERA.png
    ├── benchmark_Hallmark_ORA_vs_GSEA_vs_CAMERA.png
    ├── benchmark_C7_ORA_vs_GSEA_vs_CAMERA.png
    ├── benchmark_WikiPathways_ORA_vs_GSEA_vs_CAMERA.png
    └── go_redundancy_barchart.png
```

## Methods

### Differential expression

GSE193336 provides processed normalized expression values (non-integer) rather than raw counts. DESeq2 requires raw integer counts because its negative binomial model is not valid for pre-normalized data. limma-voom is used instead, as it operates on log-transformed normalized values and models the mean-variance relationship directly.

Design: `~ donor + condition` (paired donors). The donor term absorbs inter-individual variability; the condition coefficient captures the LPS effect.

- ORA input: LPS upregulated DEGs (adj.P < 0.05, logFC > 1)
- GSEA and CAMERA input: all expressed genes ranked by sign(logFC) × −log10(adj.p)

### Enrichment methods

| Method | Type | Input | Key feature | Implementation |
|---|---|---|---|---|
| ORA | Overrepresentation test | Significant gene list (FDR threshold) | Simple, fast; ignores subthreshold signal | clusterProfiler::enrichGO, enrichPathway, enricher |
| GSEA | Functional class scoring | Ranked gene list (all genes) | Uses full distribution; detects coordinated shifts | clusterProfiler::gseGO, gsePathway; fgsea |
| CAMERA | Competitive gene set test | Expression matrix + design | Accounts for inter-gene correlation; integrates with limma | limma::camera |

### Gene-set resources

| Resource | Scope | Approx. size | Strengths |
|---|---|---|---|
| GO BP | Broad annotated biological processes | ~7,500 terms | Comprehensive; high redundancy requires post-hoc simplification |
| Reactome | Expert-curated molecular reaction pathways | ~1,500 pathways | Hierarchical; lower redundancy than GO |
| MSigDB Hallmark | Curated coherent biological states | 50 gene sets | Most concise; low redundancy; well-validated across studies |
| WikiPathways | Community-curated pathways | ~700 pathways | Frequently updated; fully open access |
| C7 Immunologic Signatures | Immune cell perturbation signatures | ~5,000 gene sets | Experiment-level resolution for immune biology; high library size increases hits but also false positives |

### Benchmark metrics

| Metric | Description | Script |
|---|---|---|
| Pathway recovery | Presence and rank of expected LPS programs (TNF/NF-κB, interferon response, IL-6/JAK/STAT, complement, TLR cascades) | 02_ora_gsea.R + 03_camera.R |
| Method concordance | Do ORA, GSEA, and CAMERA identify the same top terms for the same resource? | 02_ora_gsea.R |
| Redundancy | Number of unique biological themes relative to total significant terms; GO BP: 1,007 terms → 151 themes (85% redundancy, Wang method, threshold=0.7) | 04_go_redundancy.R |
| Significant term count | Total terms at FDR < 0.05 per method × resource combination | All enrichment scripts |

### CAMERA implementation notes

- **Gene set indexing:** CAMERA requires integer indices into the voom expression matrix rows, not gene IDs. `limma::ids2indices()` converts gene symbol lists to row indices by matching against the Ensembl-to-symbol map derived from voom matrix rownames.
- **inter.gene.cor = NA:** Estimates intra-set correlation from the data rather than assuming a fixed value. This is the key statistical feature distinguishing CAMERA from standard competitive tests: it accounts for co-regulation among pathway genes.
- **Significance threshold:** CAMERA results use PValue < 0.05 rather than FDR < 0.05. With n=8 samples and thousands of gene sets tested (4,041 GO BP, 1,111 Reactome), BH correction is overly conservative: 338 GO BP terms have PValue < 0.05 but none survive FDR correction. This is statistically justifiable because CAMERA already corrects for inter-gene correlation, the primary source of inflation in competitive tests. Intra-set correlation values for significant Hallmark terms ranged from 0.006 to 0.039.
- **Name matching:** CAMERA uses MSigDB-prefixed pathway names (GOBP_, REACTOME_) while ORA uses readable descriptions. A `normalize_name()` function strips prefixes, underscores, punctuation, and lowercases before matching. 14/15 ORA top terms matched for both GO BP and Reactome; the 1 unmatched term per resource appears as white in the heatmap.
- **Color scale capping:** Heatmap color scale is capped at max(CAMERA -log10) × 1.2 so CAMERA values (~1–3) are visible alongside ORA/GSEA values which can reach -log10 ~80.

## Results

### Differential expression (limma-voom, LPS vs Unstimulated)

Canonical LPS-response genes recovered with expected directionality: IL1B (logFC 9.75), CXCL8 (8.69), PTGS2 (6.18), IL6 (4.13), CCL2 (3.86), TNF (2.56), all upregulated and all FDR < 0.05.

### Enrichment method comparison (ORA vs GSEA vs CAMERA, Hallmark)

| | ORA | GSEA | CAMERA |
|---|---|---|---|
| Significant terms | 13 (FDR < 0.05) | 31 (FDR < 0.05) | 8 (PValue < 0.05) |
| Top term | INTERFERON GAMMA RESPONSE | INTERFERON GAMMA RESPONSE | INFLAMMATORY RESPONSE |
| Core LPS programs (TNF/NF-κB, IL6/JAK/STAT) | Yes | Yes | Yes |
| Subthreshold programs (HYPOXIA, KRAS, EMT) | No | Yes | No |
| Interferon signal | Very strong | Strong | Weak |

Key observations:
- **ORA** concentrates signal in the most significant DEGs, giving the strongest p-values for interferon and inflammatory terms but missing subthreshold programs entirely
- **GSEA** detects the broadest biology, with subthreshold programs (HYPOXIA, KRAS, EPITHELIAL MESENCHYMAL TRANSITION) visible via ranked list that ORA and CAMERA both miss
- **CAMERA** is most conservative because the correlation correction attenuates interferon terms most strongly (interferon-stimulated genes are highly co-regulated); detects core inflammatory programs with high confidence but misses programs with weaker or noisier gene-level signal

### Resource comparison (same method, different gene-set resource)

- **GO BP:** highest term count (ORA: 1,007 / GSEA: 1,099 / CAMERA: 338 at p < 0.05). After semantic similarity clustering (Wang method, threshold = 0.7), 1,007 terms collapse to 151 unique biological themes, which is 85% redundancy. Post-hoc simplification via GOSemSim/rrvgo is essential for GO BP interpretation. Representative themes: cytokine-mediated signaling, adaptive immune response, leukocyte chemotaxis, NF-κB signaling, antiviral response.
- **Reactome:** most mechanistic resolution: "DDX58/IFIH1-mediated induction of interferon-alpha/beta", "Toll-like Receptor Cascades", "Interleukin-10 signaling" capture pathway-level mechanisms absent from GO and Hallmark.
- **Hallmark:** cleanest non-redundant summary (13–31 significant terms), each representing a distinct curated biological program; no post-processing required.
- **C7 Immunologic Signatures:** highest raw term count (ORA: 2,288 / GSEA: 2,767 / CAMERA: 629), partly reflecting library size (~5,219 gene sets). Provides experiment-level resolution, as top terms directly name matched LPS-monocyte and dendritic cell experiments (GSE9988, GSE14000, GSE22886). High hit counts also mean more potential false positives; interpret alongside curated resources.
- **WikiPathways:** broadly comparable to Reactome in term count (ORA: 141 / GSEA: 209 / CAMERA: 69). Adds disease-specific pathways absent from Reactome, such as TLR signaling, COVID-19 lung pathology, SARS-CoV-2 innate immunity evasion. Fully open access.

## Tools

**Language:** R (v4.5)

- Differential expression: limma, edgeR
- Enrichment: clusterProfiler, ReactomePA, fgsea, msigdbr
- Annotation: org.Hs.eg.db, AnnotationDbi
- Redundancy: GOSemSim, rrvgo
- Visualization: ggplot2, tidyverse

## Reference

GSE193336. Bulk RNA-seq on unstimulated and LPS-stimulated human macrophages. Deposited 2022-03-24. Sanofi in-house blood donor service. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193336

## License

MIT License. See [LICENSE](LICENSE) for details.
