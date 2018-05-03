

# input for this one is a fasta file
# VERY important to analyze the fasta containing Ns, and not the ms file, because the ms file contains "0s" where there
# should be Ns which affects the folding of the SFS

import sys

with open(sys.argv[1], 'r') as infile:
    firstHeaderInfo = infile.readline().strip().split()
    locusID = firstHeaderInfo[1]
    seqs = []
    for line in infile:
        if line[0] != ">": # sequence
            seqs.append(line.strip())
        else: # header
            newline = line.strip().split()

            if newline[1] != locusID: # end of locus, process seqs
                locusLength = len(seqs[0])
                for position in range(locusLength):
                    nucs = []
                    for seq in seqs:
                        if seq[position] != "N":
                            nucs.append(seq[position])
                    for nuc in set(nucs):
                        if len(set(nucs)) > 2:
                            pass # for now, although I ned to deal with these at some point
                        elif len(set(nucs)) == 2:
                            a1 = nucs.count(list(set(nucs))[0])
                            a2 = nucs.count(list(set(nucs))[1])
                            if a1 <= a2:
                                print a1#, nucs
                            else:
                                print a2#, nucs
                            #print len(nucs)
                            
                # restart variables
                locusID = newline[1]
                seqs = []
            
            

            
