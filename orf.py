import random

from dna import *

def restOfORF(DNA):
    """Takes a sequence starting with an ATG and finds first stop
    codon. Returns ORF. If no in-frame stop codon, return whole
    sequence."""
    for i in range(0,len(DNA),3):
        if DNA[i:(i+3)] in stopList:
            seq=DNA[:i]
            return seq
    return DNA

def oneFrame(DNA):
    """Begining at the start of DNA, searches that one frame for all
    ORFs. Returns their seqs as list."""
    seqL = []
    for i in range(0,len(DNA),3):
        if DNA[i:i+3] == "ATG":
            seq = restOfORF(DNA[i:])
            seqL = seqL + [seq]
    return seqL

def oneFrameV2(DNA):
    seqL = []
    for i in range(0,len(DNA),3):
        if DNA[i:i+3] == "ATG":
            seq = restOfORF(DNA[i:])
            seqL = seqL + [seq]
            break
    return seqL

def longestORF(DNA):
    """Finds the longest distance between a Start codon and the next
    in frame Stop. Returns this along with the corresponding DNA."""
    maxLn=0
    maxSeq=""
    for orf in oneFrame(DNA)+oneFrame(DNA[1:])+oneFrame(DNA[2:]):
        if len(orf)>maxLn:
            maxSeq=orf
            maxLn=len(orf)
    return(maxSeq)

def longestORFV2(DNA):
    """Finds the longest distance between a Start codon and the next
    in frame Stop. Returns this along with the corresponding DNA."""
    maxLn=0
    maxSeq=""
    for orf in oneFrameV2(DNA)+oneFrameV2(DNA[1:])+oneFrameV2(DNA[2:]):
        if len(orf)>maxLn:
            maxSeq=orf
            maxLn=len(orf)
    return(maxSeq)

def longestORFBothStrands(DNA):
    dnaLongestORF = longestORFV2(DNA)
    revDNA = reverseComplement(DNA)
    revDnaLongestORF = longestORFV2(revDNA)
    if (len(dnaLongestORF) > len(revDnaLongestORF)):
        return dnaLongestORF
    else:
        return revDnaLongestORF

def collapse(L):
   output = ""      # This is our initial output string
   for s in L:      # for each string in the list...
      output = output + s    #... construct a new output string
   return output    # ... and return the final output string


def longestORFNoncoding(DNA, numReps):
    listDNA = list(DNA)
    lengthNonCoding = 0
    for i in range(0, numReps):
        random.shuffle(listDNA)
        lengthNonCoding = collapse(longestORFBothStrands(collapse(listDNA)))
    return lengthNonCoding

def findORFs(DNA):
    seqL = []
    maxLn=0
    maxSeq=""
    for orf in range(3):
        orf_seq = DNA[orf:]
        orfs_frames=oneFrameV2(orf_seq)
        seqL.extend(orfs_frames)
    return seqL

def findORFsBothStrands(DNA):
    revDNA = reverseComplement(DNA)
    return oneFrame(DNA) + oneFrame(revDNA)

def getCoordinates(gene, DNA):
    pos = DNA.find(gene)
    if pos > 0:
        return pos, pos + len(gene)
    else:
        rev = reverseComplement(gene)
        posrev = DNA.find(rev)
        return posrev, posrev + len(gene)


# Compose a function that uses two functions: (1) findORFsBothStrands(DNA)
# and (2) getCoordinates(gene, DNA) to find genes in a segment of genomic
# DNA. This function should return the following: (1) a list of genes;
# (2) the start and (3) end coordinates for each gene; and (4) a sequence
# of the encoded protein. Write the function such that genes are listed
# order of their start coordinates (i.e., moving forward along the genomic
# DNA sequence)


def geneFinder(DNA, minLen):
    seq = findORFsBothStrands(DNA)
    finalList = []
    finalOutputList = []
    for orf in seq:
        ord = getCoordinates(orf,DNA)
        if ord is not None:
            start_pos, end_pos = ord
            lenorf = end_pos - start_pos + 1
            if lenorf > minLen:
                proseq = codingStrandToAA(orf)
                finalOutputList.append([start_pos,end_pos, proseq])

    finalOutputList.sort()
    return finalOutputList

def printGenes(geneL):
    print("Printing the Results: ")
    for index, gene in enumerate(geneL, start=1):
        start_pos, end_pos, proSequence = gene
        print(f"Gene {index}:")
        print(f" Start Position: {start_pos}")
        print(f" End Position: {end_pos}")
        print(f" Protein Sequence: {proSequence}")


# X73525 = loadSeq("X73525.fa")
# printGenes(geneFinder(X73525, 100))

X73525 = loadSeq("X73525.fa")
minLen = len(longestORFNoncoding(X73525, 1500))
print(minLen)
printGenes(geneFinder(X73525, minLen))



#print(len(longestORFNoncoding(X73525,10)))

#print(longestORFNoncoding(X73525,2))

##############


