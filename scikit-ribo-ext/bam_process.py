#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo 
## ----------------------------------------
## a module for processing bam files
## ----------------------------------------
## author: Han Fang 
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 10/28/2016
## ----------------------------------------

from __future__ import print_function, division
import os
import sys
import argparse
import pybedtools as pbt
import pysam
import pandas as pd
import numpy as np
import csv
import re
import errno
import gffutils
from datetime import datetime
import polars as pl
import pyarrow


class BamProcess(object):
    '''extracting alignment, prepare training/testing data
    '''
    def __init__(self, bam=None, mapq=20, minRL=10, maxRL=35, RelE=None, output=None, startCodons=None, cds=None,
                 posRanges=None, nts=None, geneSubset=None, directory=None, prefix=None):
        self.bam = bam
        self.mapq = mapq
        self.output = output
        self.minRL = minRL
        self.maxRL = maxRL
        self.startCodons = startCodons
        self.cds = cds
        self.RelE = RelE
        self.posRanges = posRanges
        self.ntsFn = nts
        self.geneSubset = geneSubset
        self.directory = directory
        self.Idxprefix = prefix

    def filterRegion(self):
        # create bedtool, filter and sort
        self.startCodons = pbt.BedTool(self.startCodons).sort()
        self.cds = pbt.BedTool(self.cds).bed6().sort()
        # find overlapping regions
        distinctStartCodons = self.startCodons.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
        distinctCDSs = self.cds.merge(c=    "1,2", o="count").filter(lambda x: int(x[4]) == 1)
        distinctStartCodons = distinctStartCodons.intersect(distinctCDSs, wa=True, sorted=True)
        # filter start codon
        startCodonHash = set([(i[0], i[1], i[2]) for i in distinctStartCodons])
        self.startCodons = self.startCodons.filter(lambda x: (x[0], x[1], x[2]) in startCodonHash)

    # TODO: add a function to see whether there are too many soft-clipped alignment
    def filterBam(self):
        # create a template/header from the input bam file
        inBam = pysam.AlignmentFile(self.bam, "rb")
        self.prefix = self.output + "/riboseq"
        outBam = pysam.AlignmentFile(self.prefix + ".bam", "wb", template=inBam)
        # read a bam file and extract info
        cigar_to_exclude = set([1,2,3,4,5]) #set(['I','D','S','H'])
        for read in inBam.fetch():
            cigars = set([c[0] for c in read.cigartuples])
            if read.mapping_quality > self.mapq and self.minRL <= read.query_length <= self.maxRL and \
            not cigars.intersection(cigar_to_exclude):
            #len(cigars.intersection(cigar_to_exclude)) == 0: # read.reference_id != 'chrM': # and 'N' not in edges:
                read.mapping_quality = read.query_length # change to bamtobed
                outBam.write(read)
        inBam.close()
        outBam.close()
        sys.stderr.write("[status]\tFinished filtering the bam file" + "\n")
        # save the bedtool to local
        self.bedtool = pbt.BedTool(self.prefix + ".bam")
        self.bedtool = self.bedtool.bam_to_bed(bed12=False)
        self.bedtool.saveas(self.prefix + '.bed')
        ## Modified to avoid an error while reading the bam file
        ## Add the following command lines to sort the bam file 
        cmd = "sort -k1,1 -k2,2n " + self.prefix + ".bed > " + self.prefix + ".bed2"
        os.system(cmd)
        cmd = "rm " + self.prefix + ".bed"
        os.system(cmd)
        cmd = "mv " + self.prefix + ".bed2 " + self.prefix + ".bed"
        os.system(cmd)
        
    def posIndex(self):
        # create a dict to store the position read-frame and index info
        self.posOffsets, self.negOffsets = [], []
        with open(self.posRanges, 'r') as f:
            next(f) # skip header
            for line in f:
                gene, chr, strand, ranges = line.rstrip("\n").split("\t")
                boxes = [(int(i[0]), int(i[1]), int(i[2])) for i in [j.split(",") for j in ranges.split("|") ]]
                if strand == "+":
                    self.posOffsets.extend([[chr, pos, (abs(pos-(box[0]-15)) + box[2])%3] for box in boxes for pos in range(box[0]-15, box[1]+12)])
                else:
                    boxes = boxes[::-1] # flip the order
                    self.negOffsets.extend([[chr, pos, (abs(pos-(box[1]+15)) + box[2])%3] for box in boxes for pos in range(box[1]+15, box[0]-12, -1)])
        # convert to dataframe
        self.posOffsets = pd.DataFrame(self.posOffsets, columns=["chrom", "pos", "offset"]).drop_duplicates(subset=["chrom", "pos"])
        self.posOffsets = pl.from_pandas(self.posOffsets)
        self.negOffsets = pd.DataFrame(self.negOffsets, columns=["chrom", "pos", "offset"]).drop_duplicates(subset=["chrom", "pos"])
        self.negOffsets = pl.from_pandas(self.negOffsets)

    def sortBam(self):
        self.bedtool = pbt.BedTool(self.prefix + ".bam")
        self.bedtool = self.bedtool.bam_to_bed(bed12=False)
        self.bedtool = self.bedtool.sort()

    def geneSubsets(self):
        # create dict of transcript ID and gene ID based on subset list provided
        time = str(datetime.now())
        sys.stderr.write("[status]\tLoading annotation and fetching gene IDs for subsetting: " + time + "\n")
        gtfDedup = self.directory + "/" + self.Idxprefix + ".dedup.gtf"
        genes = dict()
        filteredGenes = dict()
        db = gffutils.create_db(gtfDedup, ":memory:", force=True, verbose=False)
        # retrieve a list of gene names with type "CDS" from db
        for gene in db.features_of_type("CDS", order_by=("seqid","start")):
            genes[gene['transcript_id'][0]] = gene['gene_id'][0]
        with open(self.geneSubset, "r") as geneList:
            for line in geneList:
                line = line.strip()
                # different annotations means transcript names don't alway match fully after "."
                # find correct one and match to gene
                shortTransc = line.split(".")[0]
                try:
                    transcrMatch = [a for a in genes.keys() if re.search(shortTransc, a)][0]
                except IndexError:
                    sys.stderr.write("[warning]\tFollowing gene in subset list not in index: " + shortTransc + "\n")
                filteredGenes[transcrMatch] = genes[transcrMatch]

        return filteredGenes

    def makeTrainingData(self):
        # intersect with start codons
        self.bedtool = pbt.BedTool(self.prefix + '.bed')
        training = self.bedtool.intersect(self.startCodons, wa=True, wb=True, sorted=True, s=True)
        time = str(datetime.now())
        sys.stderr.write("[status]\tFinished intersecting the bedtool with start codons: " + time + "\n")
        # convert bedtool to df
        colNames = ['chrom', 'start', 'end', 'name', 'read_length', 'strand', 'sc_chrom', 'sc_start', 'sc_end', 'gene',
                    'sc_score', 'gene_strand']
        training = training.to_dataframe(names=colNames)
        # convert to polars dataframe
        training = pl.from_pandas(training)
        # a-site
        if not self.RelE:
            training['asite'] = np.where(training['gene_strand'] == '+', training['sc_start'] - training['start'] + 3,
                                         training['end'] - training['sc_end'] + 3 )
        else:
            training['asite'] = np.where(training['gene_strand'] == '+', training['end'] - training['sc_start'] - 3,
                                         training['sc_end'] - training['start'] - 3 )
        # phasing 5'
        tmpA = training.join(self.posOffsets, left_on=["chrom", "start"], right_on=["chrom", "pos"])
        tmpB = training.join(self.negOffsets, left_on=["chrom", "end"], right_on=["chrom", "pos"])
        training = pl.concat([tmpA, tmpB])
        training = training.rename({'offset':'5_offset'})
        # phasing 3'
        tmpA = training.join(self.posOffsets, left_on=["chrom", "end"], right_on=["chrom", "pos"])
        tmpB = training.join(self.negOffsets, left_on=["chrom", "start"], right_on=["chrom", "pos"])
        training = pl.concat([tmpA, tmpB])
        training = training.rename({'offset':'3_offset'})
        # filter a read by whether it has a-site that satisfies [9,18] or [1,8]
        if not self.RelE:
            training = training.filter((pl.col("asite") >= 9) & (pl.col("asite") <= 18))
            training = training.filter((pl.col("asite") >= pl.col("read_length")/2 -1))
        else:
            training = training.filter((pl.col("asite") >= 1) & (pl.col("asite") <= 5))
        # get nts
        self.nts = pl.read_csv(self.ntsFn, has_header = True, sep = "\t")
        training['pos_-1'] = np.where(training['gene_strand'] == '+', training['start']-1,  training['end'])
        training['pos_0'] = np.where(training['gene_strand'] == '+', training['start'], training['end']-1)
        training['pos_n-1'] = np.where(training['gene_strand'] == '+', training['end']-1, training['start'])
        training['pos_n'] = np.where(training['gene_strand'] == '+', training['end'], training['start']-1)
        training = training.with_column(pl.col(["pos_-1", "pos_0", "pos_n-1", "pos_n"]).cast(pl.Int64))
        # merge training with nts
        training = training.join(self.nts, left_on=["chrom", "pos_-1"], right_on=["chrom", "pos"])
        training = training.drop(["pos_-1"])
        training = training.rename({'nt': 'nt_-1'})
        training = training.join(self.nts, left_on=["chrom", "pos_0"], right_on=["chrom", "pos"])
        training = training.drop(["pos_0"])
        training = training.rename({'nt': 'nt_0'})
        training = training.join(self.nts, left_on=["chrom", "pos_n-1"], right_on=["chrom", "pos"])
        training = training.drop(["pos_n-1"])
        training = training.rename({'nt': 'nt_n-1'})
        training = training.join(self.nts, left_on=["chrom", "pos_n"], right_on=["chrom", "pos"])
        training = training.drop(["pos_n"])
        training = training.rename({'nt': 'nt_n'})
        # slice the dataframe to the variables needed for training data, removed "start_nt", "end_nt"
        training = training[["asite", "read_length", "5_offset", "3_offset", "gene_strand", "nt_-1", "nt_0", "nt_n-1", "nt_n"]]
        training = training.to_pandas()
        return training

    def makeCdsData(self):
        # create pandas df from bedtools intersect
        self.bedtool = pbt.BedTool(self.prefix + '.bed')
        cds = self.bedtool.intersect(self.cds, wa=True, wb=True, sorted=True, s=True)
        time = str(datetime.now())
        sys.stderr.write("[status]\tFinished intersecting the bedtool within cds: " + time + "\n")
        colNames = ['chrom', 'start', 'end', 'name', 'read_length', 'strand', 'gene_chrom', 'gene_start', 'gene_end',
                    'gene', 'gene_score', 'gene_strand']
        cds = cds.to_dataframe(names=colNames)
        # convert to polars dataframe
        cds = pl.from_pandas(cds)
        # subset cds data if specified
        if self.geneSubset:
            time = str(datetime.now())
            sys.stderr.write("[status]\tSubsetting CDS data by gene IDs: " + time + "\n")
            self.geneSubset = self.geneSubsets()
            genes = list(self.geneSubset.values())
            cds = cds.filter(pl.col("gene").is_in(genes))
        # phasing 5'
        tmpA = cds.join(self.posOffsets, left_on=["chrom", "start"], right_on=["chrom","pos"])
        tmpB = cds.join(self.negOffsets, left_on=["chrom", "end"], right_on=["chrom","pos"])
        cds = pl.concat([tmpA, tmpB])
        cds = cds.rename({'offset':'5_offset'})
        # phasing 3'
        tmpA = cds.join(self.posOffsets, left_on=["chrom", "end"], right_on=["chrom","pos"])
        tmpB = cds.join(self.negOffsets, left_on=["chrom", "start"], right_on=["chrom","pos"])
        cds = pl.concat([tmpA, tmpB])
        cds = cds.rename({'offset':'3_offset'})
        # get nts
        cds['pos_-1'] = np.where(cds['gene_strand'] == '+', cds['start']-1,  cds['end'])
        cds['pos_0'] = np.where(cds['gene_strand'] == '+', cds['start'], cds['end']-1)
        cds['pos_n-1'] = np.where(cds['gene_strand'] == '+', cds['end']-1, cds['start'])
        cds['pos_n'] = np.where(cds['gene_strand'] == '+', cds['end'], cds['start']-1)
        cds = cds.with_column(pl.col(["pos_-1", "pos_0", "pos_n-1", "pos_n"]).cast(pl.Int64))
        # merge cds with nt
        cds = cds.join(self.nts, left_on=["chrom", "pos_-1"], right_on=["chrom", "pos"])
        cds = cds.drop(["pos_-1"])
        cds = cds.rename({'nt': 'nt_-1'})
        cds = cds.join(self.nts, left_on=["chrom", "pos_0"], right_on=["chrom", "pos"])
        cds = cds.drop(["pos_0"])
        cds = cds.rename({'nt': 'nt_0'})
        cds = cds.join(self.nts, left_on=["chrom", "pos_n-1"], right_on=["chrom", "pos"])
        cds = cds.drop(["pos_n-1"])
        cds = cds.rename({'nt': 'nt_n-1'})
        cds = cds.join(self.nts, left_on=["chrom", "pos_n"], right_on=["chrom", "pos"])
        cds = cds.drop(["pos_n"])
        cds = cds.rename({'nt': 'nt_n'})
        # slice the dataframe to the variables needed for training data
        cds = cds[["read_length", "5_offset", "3_offset", "gene_strand", "chrom", "start", "end", "nt_-1", "nt_0", "nt_n-1", "nt_n", 'strand']]
        cds = cds.to_pandas()
        return cds


# ----------------------------------------
#  main
# ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input bam file")
    parser.add_argument("-p", help="prefix for BED/index files")
    parser.add_argument("-q", help="minimum mapq allowed, Default: 20", default=20, type=int)
    parser.add_argument("-l", help="shortest read length allowed, Default: 10", default=10, type=int)
    parser.add_argument("-u", help="longest read length allowed, Default: 35", default=35, type=int)
    parser.add_argument("-e", help="whether the sample involved RelE, Default: F", default='F', type=str)
    parser.add_argument("-o", help="output path")
    # check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()
    # process the file if the input files exist
    if (args.i != None and args.p != None and args.o != None):
        # specify inputs
        inBam = args.i
        pre = args.o + "/" + args.p
        start, cds, posIdx, nt = pre + ".start.bed", pre + ".cds.bed", pre + ".pos_ranges.txt", pre + '.nt_table.txt'
        mapq = args.q
        minRL, maxRL = args.l, args.u
        RelE = False if args.e == 'F' else True
        output = args.o
        time = str(datetime.now())
        sys.stderr.write("[status]\tStart the module at " + time + "\n")
        sys.stderr.write("[status]\tProcessing the input bam file: " + inBam + "\n")
        sys.stderr.write("[status]\tOutput path: " + output + "\n")
        sys.stderr.write("[status]\tReading the start codon BED file: " + start + "\n")
        sys.stderr.write("[status]\tReading the open reading frame codon BED file: " + cds + "\n")
        sys.stderr.write("[status]\tReading the position-phase file: " + posIdx + "\n")
        # filter alignment
        sys.stderr.write("[execute]\tKeep reads with length [" + str(minRL) + ","+ str(maxRL) + "] and mapq >= " +
                         str(mapq) + "\n")
        aln = BamProcess(inBam, mapq, output, minRL, maxRL, start, cds, posIdx, RelE, nt)
        sys.stderr.write("[execute]\tFilter out overlapping regions" + "\n")
        aln.filterRegion()
        sys.stderr.write("[execute]\tImport the position ranges" + "\n")
        aln.posIndex()
        sys.stderr.write("[execute]\tFilter out un-reliable alignment from bam files" + "\n")
        aln.filterBam()
        time = str(datetime.now())
        sys.stderr.write("[execute]\tCreate training dataframe at " + time + "\n")
        aln.makeTrainingData()
        time = str(datetime.now())
        sys.stderr.write("[execute]\tCreate CDS dataframe at " + time + "\n")
        aln.makeCdsData()
        time = str(datetime.now())
        sys.stderr.write("[status]\tBam processing finished at " + time + "\n")
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()

