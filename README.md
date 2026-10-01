# Molecular Phylogenetics Tool: MSA + Tree Construction

A command-line pipeline that takes a multi-FASTA file of homologous
sequences and produces a real, properly-methodology'd phylogenetic tree:
Clustal Omega multiple sequence alignment → distance matrix → UPGMA *and*
Neighbor-Joining tree construction → bootstrap branch-support analysis →
publication-style tree figures.

## Why this matters

"Draw a phylogenetic tree" sounds like a diagram exercise, but doing it
properly — real alignment, a defensible distance/tree-building method, and
bootstrap support so you know which branches to actually trust — is core
methodology across evolutionary biology, molecular epidemiology (tracking
viral variant spread), and comparative genomics. This tool implements that
full, correct pipeline rather than just prettifying a hand-drawn tree.

## Data

Ships with **real GenBank sequence data**: chloroplast *rpl16* gene
sequences from 7 *Opuntia* (prickly pear cactus) species — the standard
teaching dataset used throughout Biopython's own official tutorial for
alignment and phylogenetics (`data/opuntia.fasta`, via
biopython/biopython, Doc/examples/opuntia.fasta; original data from
Wallace & Cota, 1996, *Southwestern Naturalist*, GenBank accessions
AF191658–AF191665).

**To run this on a biomedical gene family instead** (e.g. comparing
hemoglobin, insulin, or a disease gene across species): use
`src/fetch_sequences.py` to pull real sequences from UniProt or NCBI at
runtime — see the docstring in that file for a ready-to-run hemoglobin-beta
example. This needs a normal internet connection (NCBI/UniProt aren't
reachable from every sandboxed environment, but work fine on a regular
laptop).

## Usage

```bash
sudo apt install clustalo   # or: conda install -c bioconda clustalo
pip install -r requirements.txt

python src/build_tree.py --fasta data/opuntia.fasta --out results/ \
    --bootstrap-replicates 100
```

Works with any multi-FASTA of homologous sequences — just point `--fasta`
at your own file.

## Method

1. **Alignment**: real Clustal Omega (Sievers et al. 2011), not a
   simplified reimplementation — this is the actual tool used in
   published research
2. **Distance matrix**: identity-based pairwise distances from the
   alignment (`Bio.Phylo.TreeConstruction.DistanceCalculator`)
3. **Tree building**: both **UPGMA** (assumes a constant molecular clock —
   simpler, but can mislead if evolutionary rates differ across lineages)
   and **Neighbor-Joining** (clock-free, generally preferred) — run both,
   see where they agree/disagree
4. **Bootstrap support**: 100 resampled alignments → 100 trees → majority-
   rule consensus, with support values on internal branches. A branch with
   <70% bootstrap support is conventionally treated as unresolved — the
   tool reports this so you know which parts of the topology to trust.

## Outputs (`results/`)

- `alignment.fasta` — the Clustal Omega multiple sequence alignment
- `tree_upgma.png` / `tree_nj.png` — the two tree-building methods
- `tree_bootstrap_consensus.png` (+ `.nwk`) — final tree with branch
  support values; the `.nwk` (Newick format) files open directly in
  iTOL, FigTree, or Dendroscope for further exploration/publication figures

## Possible extensions

- Swap identity-distance for a proper substitution model (Kimura 2-
  parameter for DNA, or a BLOSUM/PAM-based model for protein) via
  `DistanceCalculator("blosum62")` etc. — one-line change, more accurate
  for divergent sequences
- Maximum-likelihood tree building (e.g. via IQ-TREE or RAxML) as a more
  rigorous alternative to the distance-based methods here
- Ancestral sequence reconstruction at internal nodes
