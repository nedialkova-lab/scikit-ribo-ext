#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The scikit-ribo is design for analysis of yeast data, which does not have UTR annotation.
In order to run with other genomes, the boundaries of exons in gtf should be confined to CDS regions (including stop codons).
This script merges the start and stop position of the exon with the start and the stop codon of the CDS.
"""

import os
import argparse
import subprocess
from gtfparse import read_gtf


def main(args):

    # Read the gtf file
    df = read_gtf(args.gtf_file)
    #print (df)

    file_name = str(args.gtf_file)[:-4]

    list_transcripts_ID = set(list(df['transcript_id']))

    dic_start_codons_pos = {}
    dic_stop_codons_pos = {}
    dic_strand = {}
    # Here we write a new file with the principal isoforms, filtering all the transcripts
    # having more than one annotated start/stop codon, or not having a start/stop codon. 
    outfile = open(file_name + ".txt", 'w')
    for ID in list_transcripts_ID:
        # Find all the feature associated with a single transript
        df_for_transcript_ID = df[df['transcript_id']==ID]
        #print (df_for_transcript_ID)
        # Find the position of the start codon, the stop codon, and the strand of the transcripts
        start_codon_feature = df_for_transcript_ID[df_for_transcript_ID['feature']=='start_codon']
        stop_codon_feature = df_for_transcript_ID[df_for_transcript_ID['feature']=='stop_codon']
        # Skip all the transcripts without a start codon or with more than one start codon
        if start_codon_feature.shape[0] == 0 or start_codon_feature.shape[0] > 1:
            continue
        # Skip all the transcripts without a stop codon or with more than one stop codon
        if stop_codon_feature.shape[0] == 0 or stop_codon_feature.shape[0] > 1:
            continue
        strand = str(start_codon_feature['strand'].values[0])
        outfile.write("%s\n" %(ID))
        # Positive strand
        if strand == '+':
            dic_start_codons_pos[ID] = int(start_codon_feature['start'].values[0])
            dic_stop_codons_pos[ID] = int(stop_codon_feature['end'].values[0])
            dic_strand[ID] = strand
        # Negative strand
        if strand == '-':
            dic_start_codons_pos[ID] = int(start_codon_feature['end'].values[0])
            dic_stop_codons_pos[ID] = int(stop_codon_feature['start'].values[0])
            dic_strand[ID] = strand
    outfile.close()

    # Create a new gtf file with the filtered principal isoforms
    subprocess.call("grep -f " + file_name + ".txt " + args.gtf_file + " > " + file_name + "_complete.gtf", shell=True)

    # Write the file with the correct positions for each features.
    # The input file is generated using the grep function above.
    # We need to cut all the regions outside the start and the stop codon to run scikit-ribo on multicellular organisms.
    # We generate a file with two columns: corrected start position, corrected stop positions
    infile = open(file_name + "_complete.gtf", 'r')
    outfile = open('correct_positions.txt', 'w')
    for line in infile:
        transcript_ID = line.split('	')[8].split('; ')[1].split(' ')[1][1:-1]
        # Positive strand (+)
        if dic_strand[transcript_ID] == '+':
            feature_start = int(line.split('	')[3])
            feature_end = int(line.split('	')[4])
            # Feature start < start codon
            if feature_start < dic_start_codons_pos[transcript_ID]:
                new_start = dic_start_codons_pos[transcript_ID]
            # Feature start > stop codon
            elif feature_start > dic_stop_codons_pos[transcript_ID]:
                new_start = dic_stop_codons_pos[transcript_ID]
            else:
                new_start = feature_start
            # Feature end < start codon
            if feature_end < dic_start_codons_pos[transcript_ID]:
                new_end = dic_start_codons_pos[transcript_ID]
            # Feature end > stop codon
            elif feature_end > dic_stop_codons_pos[transcript_ID]:
                new_end = dic_stop_codons_pos[transcript_ID]
            else:
                new_end = feature_end
            outfile.write("%s	%s\n" %(new_start,new_end))
        # Negative strand (-)
        if dic_strand[transcript_ID] == '-':
            feature_start = int(line.split('	')[4])
            feature_end = int(line.split('	')[3])
            # Feature start > start codon
            if feature_start > dic_start_codons_pos[transcript_ID]:
                new_start = dic_start_codons_pos[transcript_ID]
            # Feature start < stop codon
            elif feature_start < dic_stop_codons_pos[transcript_ID]:
                new_start = dic_stop_codons_pos[transcript_ID]
            else:
                new_start = feature_start
            # Feature end > start codon
            if feature_end > dic_start_codons_pos[transcript_ID]:
                new_end = dic_start_codons_pos[transcript_ID]
            # Feature end < stop codon
            elif feature_end < dic_stop_codons_pos[transcript_ID]:
                new_end = dic_stop_codons_pos[transcript_ID]
            else:
                new_end = feature_end
            outfile.write("%s	%s\n" %(new_end,new_start))
    infile.close()
    outfile.close()


    # Replace the corrected start positions using awk (OFS='\t' sets the output field separator)
    subprocess.call("awk -F'\t' 'FNR==NR{a[NR]=$1;next}{$4=a[FNR]}1' OFS='\t' correct_positions.txt " + file_name + "_complete.gtf" + " > " + file_name + "_corrected_temp_1.gtf", shell=True)
    # Replace the corrected stop positions using awk (OFS='\t' sets the output field separator)
    subprocess.call("awk -F'\t' 'FNR==NR{a[NR]=$2;next}{$5=a[FNR]}1' OFS='\t' correct_positions.txt " + file_name + "_corrected_temp_1.gtf" + " > " + file_name + "_corrected_temp_2.gtf", shell=True)
    # Remove all the features non-overlapping the coding regions
    subprocess.call("awk '$4!=$5' " + file_name + "_corrected_temp_2.gtf" + " > " + file_name + "_corrected_temp_3.gtf", shell=True)
    # Remove all the features associated with 'transcript' and 'UTR'
    subprocess.call("""awk 'BEGIN { OFS=FS="\t" } $3 !~ /^(transcript|UTR)/' """ + file_name + """_corrected_temp_3.gtf > """ + file_name + """_corrected_temp_4.gtf""", shell=True)
    # Remove all the attributes in the 9th column after gene_id and transcript_id to reduce the memory
    subprocess.call("""awk 'BEGIN{FS=OFS="\t"}  {sub(/gene_type.*/,"",$9)} 1' """ + file_name + """_corrected_temp_4.gtf > """ + file_name + """.corrected_annotation.gtf""", shell=True)

    # Remove the intermidiate files
    os.remove("correct_positions.txt")
    os.remove(file_name + "_corrected_temp_1.gtf")
    os.remove(file_name + "_corrected_temp_2.gtf")
    os.remove(file_name + "_corrected_temp_3.gtf")
    os.remove(file_name + "_corrected_temp_4.gtf")
    


# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Remove UTR regions from ORFS in GTF file, keeping only features overlapping coding regions.")
    # TODO Specify your real parameters here.
    parser.add_argument("gtf_file",
                        help="Source GTF File",
                        metavar="FILE")
    args = parser.parse_args()

    main(args)


