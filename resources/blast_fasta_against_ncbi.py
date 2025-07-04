import argparse
import os
import time
from Bio.Blast import NCBIWWW
from Bio import SeqIO


def blast_sequences(fasta_file, output_folder, program='blastn', database='nt', entrez_query=None, delay=5):
   
    os.makedirs(output_folder, exist_ok=True)

    #get the total number of sequences for progress tracking
    total_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta"), 1):
        print(f"🔄 [{i}/{total_sequences}] BLASTing sequence: {record.id}")

        # check if output file already exists
        output_path = os.path.join(output_folder, f"{record.id}_blast.txt")
        if os.path.exists(output_path):
            print(f"⚠️ Output already exists for {record.id}, skipping...")
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

            time.sleep(delay)  # Delay to respect NCBI rate limits

        except Exception as e:
            print(f"❌ Error with sequence {record.id}: {e}")


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
