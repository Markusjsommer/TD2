from Bio import SeqIO
import os

def partition_fasta(input_fasta, output_dir):
    # Read all records from the input FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))
    total_records = len(records)

    # Define the percentage splits
    percentages = [20, 40, 60, 80, 100]

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Write partitions to separate files
    for percentage in percentages:
        num_records = int((percentage / 100.0) * total_records)
        output_file = os.path.join(output_dir, f"partition_{percentage}.fasta")
        with open(output_file, "w") as output_handle:
            SeqIO.write(records[:num_records], output_handle, "fasta")

    print(f"Partitions saved in {output_dir}")

if __name__ == "__main__":
    import argparse

    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Partition a FASTA file into different percentages.")
    parser.add_argument("input_fasta", type=str, help="Path to the input FASTA file")
    parser.add_argument("output_dir", type=str, help="Directory to save the partitioned FASTA files")

    args = parser.parse_args()

    # Call the partition function
    partition_fasta(args.input_fasta, args.output_dir)