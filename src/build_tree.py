"""
build_tree.py
--------------
End-to-end molecular phylogenetics pipeline:

  1. Multiple sequence alignment (MSA) via Clustal Omega (real alignment
     tool, not a toy reimplementation — install with `apt install clustalo`
     or `conda install -c bioconda clustalo`)
  2. Pairwise distance matrix from the alignment (identity-based)
  3. Tree construction: both UPGMA and Neighbor-Joining (Bio.Phylo), so you
     can compare a clock-assuming vs. clock-free method on the same data
  4. Bootstrap resampling for branch-support values (standard practice —
     an unsupported tree topology shouldn't be trusted)
  5. Tree visualization (saved as PNG + Newick for use in other tools like
     iTOL or FigTree)

Usage:
    python src/build_tree.py --fasta data/opuntia.fasta --out results/
"""
import argparse
import shutil
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_trees, majority_consensus


def relabel_with_species_names(fasta_path: str, out_dir: Path) -> str:
    """Rewrites sequence IDs from raw GenBank accessions (e.g.
    'gi|6273291|gb|AF191665.1|AF191665') to a short readable label pulled
    from the FASTA description (e.g. 'Opuntia_marenae'), purely for
    cleaner tree plots. Falls back to the original ID if no clean genus+
    species pattern is found. The alignment itself is unaffected — only
    display labels change."""
    import re

    records = list(SeqIO.parse(fasta_path, "fasta"))
    used = set()
    for rec in records:
        m = re.search(r"\b([A-Z][a-z]+)\s+([a-z][a-z\-]+)\b", rec.description)
        label = f"{m.group(1)}_{m.group(2)}" if m else rec.id
        base, i = label, 2
        while label in used:  # keep labels unique
            label = f"{base}_{i}"
            i += 1
        used.add(label)
        rec.id = label
        rec.description = ""

    relabeled_path = out_dir / "input_relabeled.fasta"
    SeqIO.write(records, relabeled_path, "fasta")
    return str(relabeled_path)


def run_clustalo(fasta_path: str, out_dir: Path) -> Path:
    if shutil.which("clustalo") is None:
        raise RuntimeError(
            "clustalo not found on PATH. Install it first:\n"
            "  Ubuntu/Debian: sudo apt install clustalo\n"
            "  conda:         conda install -c bioconda clustalo"
        )
    aln_path = out_dir / "alignment.fasta"
    subprocess.run(
        ["clustalo", "-i", fasta_path, "-o", str(aln_path), "--force", "--outfmt=fasta"],
        check=True, capture_output=True,
    )
    print(f"Wrote alignment to {aln_path}")
    return aln_path


def build_and_draw_tree(alignment, method: str, out_dir: Path, label: str):
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, method=method)  # "upgma" or "nj"
    tree = constructor.build_tree(alignment)
    tree.ladderize()

    newick_path = out_dir / f"tree_{label}.nwk"
    Phylo.write(tree, str(newick_path), "newick")

    fig, ax = plt.subplots(figsize=(9, max(4, 0.5 * len(alignment))))
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.set_title(f"{label.upper()} phylogenetic tree ({len(alignment)} taxa)")
    fig.tight_layout()
    fig.savefig(out_dir / f"tree_{label}.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {out_dir / f'tree_{label}.png'} and {newick_path}")
    return tree


def bootstrap_analysis(alignment, out_dir: Path, n_replicates: int = 100):
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, method="nj")

    trees = list(bootstrap_trees(alignment, n_replicates, constructor))
    consensus = majority_consensus(trees, cutoff=0.5)
    consensus.ladderize()

    fig, ax = plt.subplots(figsize=(9, max(4, 0.5 * len(alignment))))
    Phylo.draw(consensus, axes=ax, do_show=False)
    ax.set_title(f"Bootstrap majority-rule consensus tree ({n_replicates} replicates)")
    fig.tight_layout()
    fig.savefig(out_dir / "tree_bootstrap_consensus.png", dpi=150)
    plt.close(fig)

    Phylo.write(consensus, str(out_dir / "tree_bootstrap_consensus.nwk"), "newick")
    print(f"Wrote {out_dir / 'tree_bootstrap_consensus.png'} "
          f"(bootstrap support labels are on internal branches)")
    return consensus


def run(fasta_path: str, out_dir: str, bootstrap_replicates: int):
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"Loaded {len(records)} sequences from {fasta_path}")
    for r in records:
        print(f"  {r.id}  ({len(r.seq)} bp)  {r.description}")

    relabeled_fasta = relabel_with_species_names(fasta_path, out_path)

    print("\nRunning Clustal Omega multiple sequence alignment...")
    aln_path = run_clustalo(relabeled_fasta, out_path)
    alignment = AlignIO.read(aln_path, "fasta")
    print(f"Alignment length: {alignment.get_alignment_length()} columns")

    print("\nBuilding UPGMA tree (assumes a molecular clock)...")
    build_and_draw_tree(alignment, "upgma", out_path, "upgma")

    print("\nBuilding Neighbor-Joining tree (no clock assumption)...")
    build_and_draw_tree(alignment, "nj", out_path, "nj")

    print(f"\nRunning bootstrap analysis ({bootstrap_replicates} replicates)...")
    bootstrap_analysis(alignment, out_path, bootstrap_replicates)

    print("\nDone. See results/ for trees (.png + .nwk) and the alignment.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", default=str(Path(__file__).resolve().parent.parent / "data" / "opuntia.fasta"))
    parser.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "results"))
    parser.add_argument("--bootstrap-replicates", type=int, default=100)
    args = parser.parse_args()
    run(args.fasta, args.out, args.bootstrap_replicates)
