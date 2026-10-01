"""
fetch_sequences.py
--------------------
The bundled data/opuntia.fasta (chloroplast rpl16 gene, 7 Opuntia cactus
species, real GenBank data via the Biopython tutorial corpus) is included so
this tool has real, working example data with zero external downloads.

For a biomedical gene family instead, use this script to pull real
sequences from NCBI or UniProt at runtime (needs an internet connection —
NCBI/UniProt are not reachable from a fully offline environment, but work
fine from a normal laptop/server).

Example — human hemoglobin beta subunit (HBB), a classic comparative-
genomics teaching example, across several species:

    python src/fetch_sequences.py --db uniprot \\
        --accessions P68871,P02067,P02070,P02112 \\
        --out data/hemoglobin_beta.fasta

    # P68871 Homo sapiens (human), P02067 Sus scrofa (pig),
    # P02070 Bos taurus (cow), P02112 Gallus gallus (chicken)
    # (verify/update accessions at https://www.uniprot.org before use)

Or from NCBI, by accession or gene name search, via Biopython's Entrez:

    python src/fetch_sequences.py --db ncbi \\
        --query "HBB[gene] AND Mammalia[organism]" --retmax 10 \\
        --out data/hemoglobin_ncbi.fasta

Always set --email to a real address — NCBI requires this to identify API
traffic per their usage policy.
"""
import argparse
import urllib.request


def fetch_uniprot(accessions: list[str], out_path: str):
    with open(out_path, "w") as out:
        for acc in accessions:
            url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
            print(f"Fetching {url}")
            with urllib.request.urlopen(url) as resp:
                out.write(resp.read().decode())
    print(f"Wrote sequences for {len(accessions)} accessions to {out_path}")


def fetch_ncbi(query: str, retmax: int, email: str, out_path: str):
    from Bio import Entrez, SeqIO

    Entrez.email = email
    handle = Entrez.esearch(db="protein", term=query, retmax=retmax)
    record = Entrez.read(handle)
    ids = record["IdList"]
    print(f"Found {len(ids)} matching records for query: {query}")

    fetch_handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
    records = list(SeqIO.parse(fetch_handle, "fasta"))
    SeqIO.write(records, out_path, "fasta")
    print(f"Wrote {len(records)} sequences to {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", choices=["uniprot", "ncbi"], required=True)
    parser.add_argument("--accessions", help="Comma-separated UniProt accessions (--db uniprot)")
    parser.add_argument("--query", help="Entrez search query (--db ncbi)")
    parser.add_argument("--retmax", type=int, default=10)
    parser.add_argument("--email", default="", help="Required by NCBI for --db ncbi")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if args.db == "uniprot":
        fetch_uniprot(args.accessions.split(","), args.out)
    else:
        if not args.email:
            raise SystemExit("--email is required for NCBI Entrez queries")
        fetch_ncbi(args.query, args.retmax, args.email, args.out)
