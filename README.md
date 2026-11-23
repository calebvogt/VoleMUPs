# VoleMUPs: the single MUP gene of voles

Code behind the project page **<https://calebvogt.github.io/projects/vole-mups/>**.

House mice carry a cluster of 22 major urinary protein (*Mup*) genes and put large amounts of MUP protein in their
urine. Most other mammals keep a single ancestral *Mup* gene (Logan et al. 2008), and voles are known to use
odorant-binding proteins (OBPs) as their main urinary lipocalins (Stopková et al. 2010, 2021; Loxley et al. 2017).
This pipeline checks the vole *Mup* gene directly in open genome data: whether it is a real MUP, whether it is the
only one, and whether it is expressed in the liver.

**Result (first pass, 2026-11-13):** every one of 12 vole assemblies (10 species) has exactly one intact MUP gene, in
the ancestral position between *Zfp37* and *Slc46a2* (the genes flanking the 2.2 Mb mouse cluster); the gap-free 2026
long-read prairie vole assembly (both haplotypes) finds no second MUP locus anywhere. The vole protein is 69% identical
to mouse MUP4. In liver RNA-seq from one animal per species (n = 1), the gene gives 734 reads per million in a male
prairie vole, 4 in a male meadow vole (same lab and protocol) and 0 in a female bank vole; the mouse MUP family gives
42,260. Voles instead carry an expanded X-linked OBP1 / aphrodisin family (15 intact + 11 broken in the 2026 prairie
vole assembly; 9-25 intact in the well-assembled annotated voles; mouse 5).

## Pipeline

| Step | Script | What it does |
|------|--------|--------------|
| 01 | `scripts/step01_download.py` | RefSeq proteins + GFF3 + assembly reports for the 17 annotated genomes in `config.yml`; reference lipocalins from UniProt; the Pfam PF00061 HMM |
| 02 | `scripts/step02_annotated.py` | every gene of every annotation; HMMER `hmmsearch` PF00061 on the longest protein per gene; name matches add pseudogenes |
| 03 | `scripts/step03_regions.py --define --fetch --search` | the 8 flanking genes on each side of three mouse lipocalin clusters (MUP chr 4, OBP1 chr X, LCN2 chr 2) located by symbol in each genome; region DNA from NCBI; assembly gaps; miniprot of the query proteins to each region |
| 04 | `scripts/step04_genomes.py` | whole-genome miniprot for the 7 assemblies without an annotation |
| 05 | `scripts/step05_tree.py` | blastp labels; lipocalin tree (MAFFT, trimAl, IQ-TREE), rooted on ApoD; MUP / OBP / OBP2 clades; focused tree with ModelFinder and 1,000 ultrafast bootstraps |
| 06 | `scripts/step06_expression.py --download`, then `--count` | first 120 MB of eight liver RNA-seq runs (ENA); 25-mer read counts per lipocalin transcript, reads shared within the MUP or OBP family pooled; prairie vole MUP/OBP from its 2026 gene models, montane vole from its own MUP model and the prairie vole OBP models |
| 06b | `scripts/step06b_egapx.py` | sanity check against the Dennis lab EGAPx annotation of the 2026 prairie and meadow vole assemblies (Abuelanin et al. 2026; fetched from their GitHub, used for analysis only): gene-by-gene match to the DNA-search models, and prairie vole liver reads counted against their MUP and OBP transcripts |
| 07 | `scripts/step07_collect.py` | summary tables: genes, species, regions, CDS spans, QC |
| 08 | `scripts/step08_structures.py --alphafold --esmfold` | 3D models for the page: AlphaFold database where it has the protein, ESMFold (run locally; `pip install torch transformers`, ~8 GB of weights cached in `models/` next to `data_root`) for the prairie vole MUP; signal peptides removed |
| web | `python web/build_web_data.py`, then `python web/make_tile.py` | JSON and tile for the project page (written into a `calebvogt.github.io` checkout next to this one) |

`config.yml` sets the data folder (`data_root`, default `data/` in this repo), the species and the thread budget.
Put machine-specific changes in an untracked `config.local.yml`. Step 01 also fetches the seven assemblies without an annotation (about 6 GB) into
`<data_root>/raw/genome_<accession>/`, unless they are already in a `shared_roots` folder (`--skip-genomes` to skip).

Tools: `environment.yml` (conda/mamba: Python 3.12, Biopython, HMMER, seqkit, miniprot, MAFFT, trimAl, IQ-TREE 3,
BLAST+). Steps 02-07 need Linux or macOS; on Windows run them in WSL:

```bash
bash scripts/run_all.sh                # steps 02-07
bash scripts/run_all.sh step05 step07  # some steps
```

## Liver RNA-seq samples (one animal each)

| Species | Animal | Source | Run |
|---|---|---|---|
| prairie vole | adult male | UC Davis vole genome project, PRJEB108798 (Abuelanin et al. 2026) | ERR17322872 |
| meadow vole | adult male | same project | ERR17322876 |
| house mouse C57BL/6 | male, ZT0 | circadian study, University of Tokyo, PRJDB7789 | DRR188168 |
| bank vole | adult female | tissue-expression study, Lund University, PRJNA1119683 | SRR29282985 |

## Key references

- Abuelanin M et al. (2026) Single-library chromosome-scale diploid assemblies of vole genomes resolve a species-specific duplication implicated in pair bonding. Cell Genomics. https://doi.org/10.1016/j.xgen.2026.101336
- Logan DW, Marton TF, Stowers L (2008) Species specificity in major urinary proteins by parallel evolution. PLoS ONE 3:e3280. https://doi.org/10.1371/journal.pone.0003280
- Loxley GM et al. (2017) Glareosin: a novel sexually dimorphic urinary lipocalin in the bank vole. Open Biol 7:170135. https://doi.org/10.1098/rsob.170135
- Mudge JM et al. (2008) Dynamic instability of the major urinary protein gene family revealed by genomic and phenotypic comparisons between C57 and 129 strain mice. Genome Biol 9:R91. https://doi.org/10.1186/gb-2008-9-5-r91
- Stopková R et al. (2010) Novel OBP genes similar to hamster aphrodisin in the bank vole, *Myodes glareolus*. BMC Genomics 11:45. https://doi.org/10.1186/1471-2164-11-45
- Stopková R et al. (2021) Biological roles of lipocalins in chemical communication, reproduction, and regulation of microbiota. Front Physiol 12:740006. https://doi.org/10.3389/fphys.2021.740006
- Li H (2023) Protein-to-genome alignment with miniprot. Bioinformatics 39:btad014. https://doi.org/10.1093/bioinformatics/btad014

## License

The code is MIT-licensed (`LICENSE`). The genomes and RNA-seq are public data from NCBI and ENA, under their own terms.
