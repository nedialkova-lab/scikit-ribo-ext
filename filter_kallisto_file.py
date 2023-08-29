
__doc__= """


This script filter the output of Kallisto, taking into account only protein
#coding RNA.

Input files:    txtfile_1:           name of the txt file with the all the Ensembl transcript ID passing the filters
                                     (this file is the output of removes_3UTR_5UTR_from_GTF.py)
                txtfile_2:           name of the txt file with the association between Ensembl gene ID and transcript ID
                tsvfile:             name of the tsv file with the Kallisto TPM values
                

Options:        -h
                --help


By Sergio Forcelloni


"""


import csv
from sys import argv
from getopt import getopt


def main(argv):
    
    """
    main function
    """
    
    opts, args = getopt(argv[1:], "h", ["help"])

    print ('\n\nInput files:\n')
    for arg,i in zip(args,range(1, len(args)+1)):
        print ('Input file #%s is: %s\n' %(i,arg))

    if len(opts) > 0:
        if opts[0] == "-h" or opts[0] == "--help":
            print (__doc__)
            exit(0)       

    # Check if all the input files are provided by the user
    if len(args) < 3:
        print (__doc__)
        exit(0)

    # Reading input files
    txtfile_1 = args[0]
    txtfile_2 = args[1]
    tsvfile = args[2]

    # Dafault parameters
    prefix = "out"

    # Reading command line options
    for opt in opts:
        if opt[0] == "-p" or opt[0] == "--prefix":
            prefix = opt[1]
            
    list_ids = []
    infile = open(txtfile_1,'r')
    for line in infile:
        list_ids.append(line.strip())
    infile.close()


    # Definition of the dictionary that associates the ensembl gene ID to the corresponding transcript ID.
    dic_IDs = {}
    fi = open(txtfile_2,'r')
    for line in fi:
        splitted_line = line.strip().split()
        transcript_ID = splitted_line[1]
        gene_ID = splitted_line[0]
        try:
            dic_IDs[transcript_ID]
        except KeyError:
            dic_IDs[transcript_ID] = gene_ID


    read_tsv = csv.reader(open(tsvfile), delimiter="\t")
    outfile = open(tsvfile.split('.')[0]+"_protein_coding.tsv", 'w')
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['target_id','length','eff_length','est_counts','tpm'])
    for row in read_tsv:
        if 'target_id' in row: continue
        transcript_id = row[0]
        try:
            tpm = float(row[4])
        except ValueError: continue
        length = float(row[1])
        if tpm > 0 and transcript_id in list_ids:
            tsv_writer.writerow([dic_IDs[transcript_id],float(row[1]),float(row[2]),float(row[3]),float(row[4])])
    outfile.close()


if __name__=="__main__":
    main(argv)
