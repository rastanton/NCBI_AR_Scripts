from Bio import SeqIO
import sys
import glob

##Written by Rich Stanton (njr5@cdc.gov)
##Removed the Saute-assembled contigs from an assembly downloaded from NCBI
##Usage: python Saute_Remover.py My_Fasta

def Saute_Remover(input_genome):
    Genome = list(SeqIO.parse(input_genome, 'fasta'))
    Out = open(input_genome, 'w')
    for contig in Genome:
        Description = contig.description
        if ('guided' in Description) == False:
            SeqIO.write(contig, Out, 'fasta')
    Out.close()

def Length_Adder(input_genome):
    Genome = list(SeqIO.parse(input_genome, 'fasta'))
    Out = open(input_genome, 'w')
    for contig in Genome:
        Length = len(str(contig.seq))
        Description = contig.description
        if ('guided' in Description) == False:
            contig.id = contig.id + '__' + str(Length)
            SeqIO.write(contig, Out, 'fasta')
    Out.close()

Length_Adder(sys.argv[1])
