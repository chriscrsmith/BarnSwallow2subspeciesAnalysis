


# this script is for parsing the blast output from blasting swallow contigs against flycatcher scaffolds

# blastn -query ../BarnSwallow/masked_high_coveragesoap_assembly1_k47_min1000bp_scafs.fasta -db CombinedReference/flycatcherIncludingMT_plusChickenW.fa -outfmt "6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue gaps" -perc_identity 75 > New/blastOut.txt


import sys

# first read in chromosomes
chromDict = {}
with open(sys.argv[1], "r") as infile:
    for line in infile:
        scaff, chrom = line.strip().split()
        chromDict[scaff] = chrom



# now parse blast output
with open(sys.argv[2], "r") as infile:

    # to get things rolling, first read a line and add it to the dictionary
    qseqid, sseqid, pident, qlen, slen, length, qstart, qend, sstart, send, evalue, gaps = infile.readline().strip().split()
    currentScaffold = qseqid
    hits = {}
    hits[chromDict[sseqid]] = [ [qseqid, sseqid, pident, qlen, slen, length, qstart, qend, sstart, send, evalue, gaps] ]

    # then go through the rest of the lines
    for line in infile:
        qseqid, sseqid, pident, qlen, slen, length, qstart, qend, sstart, send, evalue, gaps = line.strip().split()
        evalue = float(evalue)
        if sseqid not in chromDict: # checking that the subject sequence is accounted for
            print "\n\nHouston we have a problem.\n\nspecifically, there's a rogue subject sequence\n\n"
            1/0
        else:
            if qseqid == currentScaffold:
                if evalue < 1e-50:
                    if chromDict[sseqid] not in hits:
                        hits[chromDict[sseqid]] = []
                    hits[chromDict[sseqid]].append( [qseqid, sseqid, pident, qlen, slen, length, qstart, qend, sstart, send, evalue, gaps] )
            else:
                # process the hits
                if len(hits) > 0:
                    chromsAligned = {}
                    bestLength = 0
                    for chrom in hits: # go through and find the chromosome with the longest alignment length
                        totAlignmentLengthOnThisChrom = 0
                        for hit in hits[chrom]:
                            qlen, alignLen = float(hit[3]), float(hit[5])
                            totAlignmentLengthOnThisChrom += alignLen
                        chromsAligned[chrom] = totAlignmentLengthOnThisChrom
                        if totAlignmentLengthOnThisChrom > bestLength:
                            bestLength = totAlignmentLengthOnThisChrom
                    if bestLength / qlen > 0.5: # make sure at least half of the swallow scaffold aligns (liberal, I would say)
                        outputIt = True
                        for chrom in chromsAligned:
                            if chromsAligned[chrom] == bestLength and chrom not in map(str,range(29)):
                                outputIt = False # checking if the best length is associated with Z, W, or MT chromosomes (or MT genome)
                        if outputIt == True:
                            print currentScaffold
                        

                # reset these variables
                hits = {}
                currentScaffold = qseqid
                if evalue < 1e-50:
                    hits[chromDict[sseqid]] = [ [qseqid, sseqid, pident, qlen, slen, length, qstart, qend, sstart, send, evalue, gaps] ]






