import sys

import re
import pysam
import pandas as pd
import numpy as np

from collections import defaultdict

class ReadAssigner:

    def __init__(self, input_file: str, cell_types: list, output_file: str, save_meta: bool = False):
        
        self.save_meta = save_meta

        self.isoform_counts = defaultdict(int)
        self.cell_types = cell_types
        self.cell_ids = []

        if input_file.endswith(".bam"):
            self.read_bam_file(bam_file = input_file)
        elif input_file.endswith(".fastq"):
            self.read_fastq_file(fastq_file = input_file)
        elif input_file.endswith(".fasta"):
            self.read_fasta_file(fasta_file = input_file)
        else:
            print("unrecognized file type at INPUT:", input_file)
        
        # Possible speedup
        print("Sorting reads... ", end="")
        self.cell_ids.sort()
        print("Done!")

        self.write_assignement(output_file = output_file)

    
    def read_fasta_file(self, fasta_file: str) -> None:
        print("Reading fasta-file...", end = "")
        fastqdata = open(fasta_file, "r")
        for line in fastqdata:
            if line[0] == '>':
                isoform_id = re.findall(r"ENSMUST[0-9.]+", line)
                if isoform_id == []:
                    continue
                self.cell_ids.append((isoform_id[0], line[1:].strip("\n")))
                self.isoform_counts[isoform_id[0]] += 1
        fastqdata.close()
        print("Done!")

    def read_fastq_file(self, fastq_file: str) -> None:
        print("Reading fastq-file...", end = "")
        fastqdata = open(fastq_file, "r")
        for line in fastqdata:
            if line[0] == "@":
                isoform_id = re.findall(r"ENSMUST[0-9.]+", line)[0]
                self.cell_ids.append((isoform_id, line[1:].strip("\n")))
                self.isoform_counts[isoform_id] += 1
        fastqdata.close()
        print("Done!")


    def read_bam_file(self, bam_file: str) -> None:
        print("Reading bam-file...", end = "")
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for read in samfile.fetch():
            isoform_id = re.findall(r"ENSMUST[0-9.]+", read.query_name)[0]
            self.cell_ids.append((isoform_id, read.query_name))
            self.isoform_counts[isoform_id] += 1
        samfile.close()
        print("Done!")

    def create_isoform_distribution(self) -> pd.DataFrame:

        distributions = np.random.random(size = (len(self.isoform_counts), len(self.cell_types)))
        distributions = distributions / distributions.sum(axis = 1).reshape(-1, 1)
        isoform_distributions = pd.DataFrame(distributions, index = self.isoform_counts.keys(), columns = self.cell_types)

        return isoform_distributions


    def write_assignement(self, output_file: str) -> None:

        print("Creating Isoform distributions... ", end = "")
        read_groups = pd.DataFrame(columns = ["read_id", "cell_type"])
        cell_type_counts = pd.DataFrame(np.zeros(shape = (len(self.isoform_counts), len(self.cell_types))), index = self.isoform_counts.keys() ,columns = self.cell_types, dtype = int)
        isoform_distributions = self.create_isoform_distribution() 
        print("Done!")

        print("Assigning read groups... ", end = "")
        for i, cell_info in enumerate(self.cell_ids):
            isoform_id, cell_id = cell_info
            cell_type = str(np.random.choice(self.cell_types, p = isoform_distributions.loc[isoform_id]))
            read_groups.loc[i] = [cell_id, cell_type]
            cell_type_counts.loc[isoform_id, cell_type] += 1
        print("Done!")

        print("Saving results... ", end = "")
        isoform_distributions.to_csv("distributions." + output_file, sep = '\t')
        cell_type_counts.to_csv("counts." + output_file, sep = '\t')
        read_groups.to_csv(output_file, header = False, columns = None, index = None, sep = '\t')
        print("Done!")

def main():

    filename, groups, outputfile = sys.argv[1:4]
    
    ReadAssigner(
        input_file = filename,
        cell_types = groups.split(":"),
        output_file = outputfile,
        save_meta = False
    )


if __name__ == "__main__":
    main()