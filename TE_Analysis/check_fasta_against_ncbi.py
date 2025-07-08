import os
import time
import argparse
from Bio import SeqIO
from Bio.Blast import NCBIWWW

def blast_sequences(fasta_file, output_folder, program='blastn', database='nt', entrez_query=None, delay=5):
    os.makedirs(output_folder, exist_ok=True)

    print(f"🔍 Starting BLAST for sequences in {fasta_file}...")

    total_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

    if total_sequences == 0:
        print(f"❌ No sequences found in {fasta_file}. Exiting.")
        return

    print(f"🔍 Found {total_sequences} sequences in {fasta_file}. Starting BLAST...")

    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta"), 1):
        print(f"🔄 [{i}/{total_sequences}] BLASTing sequence: {record.id}")
        
        output_path = os.path.join(output_folder, f"{record.id}_blast.txt")
        
        # Check if output already exists
        if os.path.exists(output_path):
            print(f"⚠️ Output already exists for {record.id}, skipping...")
            continue

        # Check for empty sequence
        if len(record.seq) == 0:
            warning = f"⚠️ Skipping {record.id}: Empty sequence.\n"
            with open(output_path, "w") as out_handle:
                out_handle.write(warning)
            print(warning.strip())
            continue

        # Check for single-letter sequences
        unique_bases = set(str(record.seq).upper())
        if len(unique_bases) == 1:
            warning = f"⚠️ Skipping {record.id}: Sequence consists of a single nucleotide: {unique_bases.pop()}.\n"
            with open(output_path, "w") as out_handle:
                out_handle.write(warning)
            print(warning.strip())
            continue

        try:
            result_handle = NCBIWWW.qblast(
                program=program,
                database=database,
                sequence=str(record.seq),
                entrez_query=entrez_query,
                format_type="Text",
                expect=0.001,
                hitlist_size=50,
                megablast=True
            )

            with open(output_path, "w") as out_handle:
                out_handle.write(result_handle.read())

            print(f"✅ Saved: {output_path}")
            result_handle.close()

            time.sleep(delay)  # Respect NCBI rate limits

        except Exception as e:
            print(f"❌ Error with sequence {record.id}: {e}")

    print(f"✅ All BLASTs completed for {fasta_file}.")


def main():
    parser = argparse.ArgumentParser(description="BLAST sequences in a FASTA file against NCBI web database.")
    parser.add_argument("fasta_file", help="Input FASTA file with sequences to BLAST")
    parser.add_argument("output_folder", help="Directory to save BLAST results (XML format)")
    parser.add_argument("--program", default="blastn", choices=["blastn", "blastp", "blastx", "tblastn", "tblastx"],
                        help="BLAST program to use (default: blastn)")
    parser.add_argument("--database", default="nt", help="NCBI database to search (default: nt)")
    parser.add_argument("--entrez_query", default=None, help="Optional Entrez query to restrict search")
    parser.add_argument("--delay", type=int, default=5, help="Delay in seconds between BLAST requests (default: 5)")

    args = parser.parse_args()

    blast_sequences(
        fasta_file=args.fasta_file,
        output_folder=args.output_folder,
        program=args.program,
        database=args.database,
        entrez_query=args.entrez_query,
        delay=args.delay
    )


if __name__ == "__main__":
    main()
