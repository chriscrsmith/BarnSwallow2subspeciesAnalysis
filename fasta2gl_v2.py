
# v2, this script, is different because it deals with (1) missing data (N's), and (2) different header format

# this will take in 1 a sample file with populations, and 2 a fasta file, and output a file with 0,1,2 genotypes

import sys

popDict = {}
with open(sys.argv[1], "r") as infile:
    for line in infile:
        sample, pop = line.strip().split()[:]
        popDict[sample] = pop
outline = []
for s in popDict:
    outline.append(s)
print '\t'.join(outline) # header




with open(sys.argv[2], "r") as infile:
    newline = infile.readline().strip().split()
    locusID = newline[1] 
    sample = newline[2]
    seqs = []
    correspondingSamples = [sample]
    for line in infile:
        newline = line.strip().split()
        if len(newline) == 1: # sequence
            seqs.append(newline[0])
        else: # header
            if newline[1] == locusID:
                sample = newline[2]
                correspondingSamples.append(sample)
            else: # found new locus, process sequences
                # first, go through the seqs and see what nucleotides occur at each position
                nucs = []
                for pos in range(len(seqs[0])):
                    newNucs = []
                    for s in range(len(seqs)):
                        nuc = seqs[s][pos]
                        if nuc not in newNucs and nuc != "N" :
                            newNucs.append(seqs[s][pos])
                    nucs.append(newNucs)
                
                # check for snps, and output
                for pos in range(len(seqs[0])):
                    if len(nucs[pos]) > 1:
                        outline = ["NA"] * len(popDict)
                        outInd = 0
                        for s in popDict:
                            if s in correspondingSamples:
                                ind = correspondingSamples.index(s)
                                theirAlleles = [ seqs[ind][pos]  , seqs[ind+1][pos]  ]
                                countNs = theirAlleles.count("N")
                                if countNs == 0:
                                    geno = 0
                                    for a in theirAlleles:
                                        if a == nucs[pos][1]: # choosing random allele to be "alt"
                                            geno += 1
                                elif countNs == 2:
                                    geno = "NA"
                                else:
                                    print " \n\n Houston we have a problem \n\n "
                                    1/0                                
                            outline[outInd] = geno
                            outInd += 1
                        print '\t'.join(map(str, outline))

                # reset these variables
                seqs = []
                locusID = newline[1]
                sample = newline[2]
                correspondingSamples = [sample]




                    
                    
                    
            
            
                
