#!/usr/bin/env python3

#################################################################################################
# 
#   Author: Carlos Arevalo (caeareva)
#
#   File: covidFasta19.py
#   Executable: python fastaFinder.py > sars2Seq.fa 
#               python fastaFinder.py > sarsG13Seq.fa
#               python fastaFinder.py > coronavirusUrbaniSeq.fa
#               python fastaFinder.py > sarsRs4231Seq.fa
#               python fastaFinder.py > coronavirusBtRsSeq.fa
#
#   Required module: pybedtools.py   
#   Pupose: Obtain fasta sequences for open reading frames (ORF) in any genome and output file
#   Condition: Program requires to change bed and fasta genome files internally each time before running it
#   
#################################################################################################

import pybedtools
from pybedtools import BedTool
import sys

class getFasta:
    '''
    Class takes a bedfile input, extract ORF coordinates, and the fasta sequences.
    The program outputs a file containing the raw fastas sequences for each gene found in the genome
    '''
    
    def __init__(self):
        '''
        Initializes program and creates empty lists 
        ''' 
        self.outJunctions = [] # stores intron junctions
        self.orfFasta = {}
        self.filteredFrames = {}
        self.orfElements = {}

    def fromBedtoFasta(self):
        '''
        Read over bed file and create a bed file containing ORF coordinates. 
        Then,uses the ORF coordinates in new bed file and a reference genome to obtain 
        the fasta for each ORF. The original bed file has to be in a matrix array form:

        Example bedfile:
                            NC_004718.3   265 13413 13149 ORF# +1
                            NC_004718.3 13597 21485  7889 ORF# +3
                            NC_004718.3 21490 25259  3770 ORF# +3
        intronCoord bed file:
                            NC_004718.3 13597   21485
                            NC_004718.3 21490   25259
                            NC_004718.3 28120   29388
        newFasta file:
                            
                            NC_004718.3:2992-3295   
                            ATGATGCTGGTGAAGAAAACTTTTCATCACGTATGTATTGTTCCTTTTACCCTCCAGATGAGGAAGAAGAGGACGATGCAGAGT
                            GTGAGGAAGAAGAAATTGATGAAACCTGTGAACATGAGTACGGTACAGAGGATGATTATCAAGGTCTCCCTCTGGAATTTGGTG
                            CCTCAGCTGAAACAGTTCGAGTTGAGGAAGAAGAAGAGGAAGACTGGCTGGATGATACTACTGAGCAATCAGAGATTGAGCCAG
                            AACCAGAACCTACACCTGAAGAACCAGTTAATCAGTTTACTGGTTATTTAA.........
        ''' 

        self.bedFile = "sars2.bed" 
        self.referenceGenome = "SARSCov2.fa" 
        orfPositions = []
        orfFasta = []
        filteredFasta = []
        with open(self.bedFile, "r") as fh:
            # open file and go through each line, i.e. rows
            for lines in fh: 
                # split line into columns and append chromosome, start and end.
                columns = lines.rstrip().split()
                # defines each column
                ID, start, end, length, orf, frame = columns[:6]
                # map the start and end of each chromosome in each sample from bed file
                start, end, length = map(int, [start, end, length])
                # new bed file containing intron coordinates is stored in intronCoord
                orfCoordinate = ('\t'.join(columns[:3]))
                pass
                
                # newBed stores the intron coordinates
                newBed = pybedtools.BedTool(orfCoordinate, from_string=True)
                # stores the fastas for each intron coordinate
                seq = pybedtools.BedTool(self.referenceGenome)
                # mathches intron coordinates with fastas and stores them  
                newBed = newBed.sequence(fi=seq) #, tab=True) 
                # object stores introns fastas
                outFasta = open(newBed.seqfn).read()
                print(outFasta)
                
                element = outFasta.rstrip().split()
                header, sequence = element
                for header in element:
                    orfPositions.append(header)
                    for sequence in element:
                        orfFasta.append(sequence)
                        pass
            
                # create objects to store new lists for introns headers and sequences 
                newHeader = orfFasta[::2] # dictionary keys 
                newSeq = orfFasta[1::2] # dictionary values 
                # create a dictionary for introns coordinates and sequences 
                for key in newHeader:
                    for value in newSeq:
                        self.orfFasta[key] = value
                        newSeq.remove(value)
                        break

#################################################################################################   
# Main function
#################################################################################################  

def main ():
    '''
    Calls functions and return output file 
    '''
    myORF = getFasta()
    fasta = myORF.fromBedtoFasta() 
    print(fasta)

if __name__ == "__main__":
    main()
  



