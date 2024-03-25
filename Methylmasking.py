#!/usr/bin/env python3
import sys
import pysam

def process_and_write_bam(input_bam_path, output_bam_path):
    # Open input BAM file for reading
    with pysam.AlignmentFile(input_bam_path, "rb") as input_bam:
        # Create a new BAM file for writing
        header = input_bam.header
        with pysam.AlignmentFile(output_bam_path, "wb", header=header) as output_bam:
            for chr in range(36,55,1):
                segments = input_bam.fetch(f"OX4570{chr}.1")
                for aligned_segment in segments:
                    meth_tag = aligned_segment.get_tag('XM')
                    sequence = list(aligned_segment.get_forward_sequence())
                    for i, char in enumerate(meth_tag):
                        if char != '.':
                            #Replace character in the same position in the sequence with 'N'
                            sequence[i] = 'N'
                    # Update the read's sequence
                    aligned_segment.query_sequence = ''.join(sequence)
                    meth_tag = str("test")

                    # Write the modified read to the new BAM file
                    output_bam.write(aligned_segment)

def process_multiple_samples(sample_file):
    with open(sample_file, 'r') as f:
        for line in f:
            input_bam, output_bam = line.strip().split()
            process_and_write_bam(input_bam, output_bam)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 script.py {input_bam_path},{output_bam_path")
        sys.exit(1)

    input_bam_path = sys.argv[1]
    output_bam_path = sys.argv[2]


    process_and_write_bam(input_bam_path, output_bam_path)
