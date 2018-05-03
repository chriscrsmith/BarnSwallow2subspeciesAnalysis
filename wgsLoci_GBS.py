



# big script for processing sam files and finding well covered continuous sequences to output

import sys, numpy
sys.stderr.write("\n\n" +  
                 "\tNOTE:\n" +
                 "If these sequences are used in an ABC analysis, you can't use (i) recombination, nor " +
                 "(2) haplotype diversity statistics. Because we can't phase the data. " +
                 "This is because we don't have enough coverage to confidently gentoype all snps, which leads to " +
                 "phasing issues when a sample is inferred to be a heterozygote but doesn't have both alleles in the reads." + 
                 "\n    also, currently, my script doesn't deal with multiple alternate alleles"
                 "\n\n")


###################
# preliminary steps 
###################

# parameters
goodScaffoldsPath = sys.argv[1] # list of scaffolds that aligned to flycatcher on non sex chromosomes
referencePath = sys.argv[2]
vcfPath = sys.argv[3]
samDirectory = sys.argv[4]
populationMapPath = sys.argv[5] 
paralogThresholdPath = sys.argv[6]
minReadDepthPerSample = int(sys.argv[7]) # minimum read depth ("bad" reads excluded) to say a sample is covered
minSamplesEachPop = int(sys.argv[8]) # minimum samples covered in each pop to output a genomic position
minimumLocusLength = int(sys.argv[9]) # the minimum locus length. Toggle this to see how it affects the total length of outputted sequence
maximumLocusLength = int(sys.argv[10]) # we're not including intra-locus recombination in simulations, so need to keep them short
minDistanceBetweenLoci = int(sys.argv[11]) # loci are simulated entirely independently, so need to insure the observed loci are independent
outputType = sys.argv[12] # what type of output. Options: "locusLengths", "fasta", "locfile", "missingData", "vcf"

# variable initialization
minSnpQUAL = 0 # minimum QUAL in .vcf to include a snp. I decided to role without this filter (that's why it's 0)
locusCounter = 0 # this is only used if outputting fasta. Included at the start, here, because the output code is inside a loop
outputLocusNum = 0 # this is only used if outputting locfile
alreadyOutputLocfileHeader = False 



# define function for getting snps from a certain scaffold from a vcf
def getSnps(vcfFile, sampleIndices, scaf, previousLine):
    data = {}
    # check on previous line, saved from when we looked at the vcf last
    if previousLine[0] == scaf:
        if float(previousLine[5]) >= minSnpQUAL:
            if len(previousLine[4]) > 1:
                genotypeData = previousLine[9:]
                alleleCollection = []
                for g in genotypeData:
                    for al in g.split(":")[0].split("/"):
                        if al != ".":
                            alleleCollection.append(al)
                AC0 = alleleCollection.count("0")
                AC1 = alleleCollection.count("1")
                AC2 = alleleCollection.count("2")
                AC3 = alleleCollection.count("3")
                AFs = [AC0,AC1,AC2,AC3]
                previousLine[2] = AFs
            data[previousLine[1]] = previousLine
        for line in vcfFile:
            newline = line.strip().split()
            if newline[0] == scaf:
                if float(newline[5]) >= minSnpQUAL: 

                    # here, checking for multiple alts, and then saving the alle frequencies (this gets complicated)
                    if len(newline[4]) > 1:
                        genotypeData = newline[9:]
                        alleleCollection = []
                        for g in genotypeData:
                            for al in g.split(":")[0].split("/"):
                                if al != ".":
                                    alleleCollection.append(al)
                        AC0 = alleleCollection.count("0")
                        AC1 = alleleCollection.count("1")
                        AC2 = alleleCollection.count("2")
                        AC3 = alleleCollection.count("3")
                        AFs = [AC0,AC1,AC2,AC3]
                        newline[2] = AFs # clever! this is storing extra AF information in the VCF info where there was a blank column 
                    data[newline[1]] = newline # read in the whole line, because looking at all samples, sometimes
            else:
                if newline[0] in goodScaffolds:
                    break
                else:
                    for line in vcfFile:
                        newline = line.strip().split()
                        if newline[0] in goodScaffolds:
                            break
        return data, newline # returning the last line read as the new "previousLine"
    else: # it's not the scaffold you want, then it's not in vcf ("previousLine" is really the current line)
        return data, previousLine

            
    


# define function for checking individual sam read
def checkRead(r):
    verdict = "bad" # default setting
    flags, qual, cigar = int(r[1]), int(r[4]), r[5]
    if qual == 60: # only considering top quality reads
        bits = map(int,list(bin(flags)[2:])[::-1])
        # flags = bits + (12-len(bits))*[0] # just adding on zeros as place holders    
        # 0. template has multiple segments (paired end reads), invariable, it's always 1 
        # 1. each segment properly aligned
        # 2. segment unmapped
        # 3. next segment unmapped
        # 4. seq reverse complemented
        # 5. next segment reverse complemented
        # 6. first segment
        # 7. second segment
        # 8. secondary alignment (?)
        # 9. not passing filters (?)
        # 10. PCR duplicate
        # 11. supplementary alignment
        # going to pretend the last several don't happen without a "segment unmapped" flag too
        if ("H" not in cigar) and ("S" not in cigar): # these are (hard, or soft) clipped reads. Other cigar options were not found in our reads
                #if ("D" not in cigar) and ("I" not in cigar): # indels
            verdict = "good!"
    return verdict



# define function for getting alignment data for a certain scaffold from the individual sams
def getSamData(sams, scaf, lastLinesRead):
    data = []
    newLinesRead = []
    for s in range(numSamples):
        if lastLinesRead[s][2] == scaf: # if the current scaffold is the correct one
            if checkRead(lastLinesRead[s]) == "good!": # checking that the read was mapped properly 
                data.append( [ lastLinesRead[s] ] ) 
            else:
                data.append( [] ) # in this case, we're on the correct scaffold, just the first read wasn't good
            for line in sams[samples[s]]:
                newline = line.strip().split()
                if newline[2] == scaf:
                    if checkRead(newline) == "good!":
                        data[s].append(newline)
                else:
                    break 
            if newline[2] in goodScaffolds: 
                newLinesRead.append(newline)
            else: # want to keep parsing until reach a line in good scaffolds list 
                for line in sams[samples[s]]:
                    newline = line.strip().split()
                    if newline[2] in goodScaffolds:
                        break
                newLinesRead.append(line.strip().split()) 
        else: # the current position in this samples sam file is not the right scaffold, i.e. sample doesn't have reads on the scaffold
            data.append(None) # placeholder for sample without data on this scaffold
            newLinesRead.append(lastLinesRead[s])
    return data, newLinesRead



# define function for genotyping multiple alternate alleles using the PL field info
def genotypeMultipleAlts(PLfieldInfo, r, a, calledGenotype):
    best = min(PLfieldInfo)
    if PLfieldInfo.count(best) > 1: # if the PL fields are ambiguous, just output N because that shits way too complicated to be worth it
        newGeno = "NN"
    else:
        bestInd = PLfieldInfo.index(best)
        assignGenotype = False
        if calledGenotype == "0/2":
            if bestInd == 3:
                calls = [0,0]
                assignGenotype = True
            else:
                newGeno = "NN"
        elif calledGenotype == "1/2":
            if bestInd == 4:
                calls = [1,1]
                assignGenotype = True
            else:
                newGeno = "NN"
        elif calledGenotype == "0/3":
            if bestInd == 6:
                calls = [0,0]
                assignGenotype = True
            else:
                newGeno = "NN"
        elif calledGenotype == "1/3":
            if bestInd == 7:
                calls = [1,1]
                assignGenotype = True
            else:
                newGeno = "NN"
        else:
            print "/n/nHouston we have a problem./n/n"
        if assignGenotype: # genotype not assigned yet
            newGeno = ""
            for g in calls:
                if g == 0:
                    newGeno += r # ref
                else:
                    newGeno += a.split(',')[g-1] # alt
    return newGeno








###################
# start doing stuff
###################
# read in the list of good scaffolds
goodScaffolds = {}
with open(goodScaffoldsPath, "r") as infile:
    for line in infile:
        if line[0] != "#":
            goodScaffolds[line.strip()] = 0



# read in the reference                          
reference = {}
scaffoldOrder = []
with open(referencePath, "r") as infile:
    for line in infile:
        newline = line.strip()
        if newline[0] == ">":
            scaffold = newline[1:]
            scaffoldOrder.append(scaffold)
        else:
            reference[scaffold] = newline



# build population map
popMap = {} # keys are samples
reversePopMap = {} # keys are pops, used separately, for speed
pops = []
samples = [] # this sample list will be in the order specified in the pop map file
with open(populationMapPath, "r") as infile:
    for line in infile:
        newline = line.strip().split()
        samples.append(newline[0])
        popMap[newline[0]] = newline[1]
        if newline[1] not in pops:
            pops.append(newline[1])
            reversePopMap[newline[1]] = []
        reversePopMap[newline[1]].append(newline[0])
numSamples = len(samples)
numPops = len(pops)


            
# open vcf
vcf = open(vcfPath, "r")
newline = "##"
while newline[0:2] == "##" :
    newline = vcf.readline()
bamList = newline.strip().split()[9:] # header
sampleVCFindices = {} # indices in the (large) bam list for the (smaller number of) samples you're analyzing
for position in range(len(bamList)):
    newID = bamList[position].split('.')[0].split("/")[-1]
    if newID in samples:
        sampleVCFindices[newID] = position+9 # saving positions of each sample in the vcf, +9 for snp info fields
lastVcfLineRead = vcf.readline().strip().split() # need a default value for this variable, the first snp line



# open sams
samDict = {} # trying to fill a directory with open bam files (leaving open to keep progress in reading each file) [yes, this is unusual, and works!]
lastSamLinesRead = [] # keeping track of the last line read, so I can leave each file open and read quickly
for s in range(numSamples):
    indFile = samples[s] + ".sam"
    samDict[samples[s]] = open(samDirectory+indFile, "r")
    for line in samDict[samples[s]]:
        if line[0] != "@": # simply getting the header lines out of the way real quick
            newline = line.strip().split()
            if newline[2] in goodScaffolds: # if scaffold 1, etc, aren't in the list of good scaffolds, keep looking
                break
    lastSamLinesRead.append(line.strip().split()) # this list is ordered by the order of the "samples" list



# read in the paralog thresholds
paralogThresholds = {}
with open(paralogThresholdPath, "r") as infile:   # these thresholds were obtained using simulations
    for line in infile:                           # specifically, these are the proportion of heterozygotes that are expected to be real 10% of the time.
        numSam, thresh = line.strip().split()[:]  # i.e. If the proportion of heterozygotes is greater than this, it's likely 
        paralogThresholds[numSam] = float(thresh) # not real (is a paralog). And if the proportion  of heterozygotes is less, then it may be real. 








        
##########
# the meat
##########
for scaffold in scaffoldOrder:
    if scaffold in goodScaffolds:
        # read in snp info
        snps, lastVcfLineRead = getSnps(vcf, sampleVCFindices, scaffold, lastVcfLineRead)

        # next go through one read at a time, (1) counting the reads per each position and (2) saving the genotype for each sample
        samData, lastSamLinesRead = getSamData(samDict, scaffold, lastSamLinesRead)
        scaffoldLength = len(reference[scaffold])
        alignmentMatrix = [] # these two matrixes will have rows as long as the scaffold, but indexed so that [0] is position 1 in the sam and vcf
        coverageMatrix = []
        for sample in range(numSamples):
            alignmentMatrix.append( ['NN']*scaffoldLength ) # initializing the row for this sample
            coverageMatrix.append( [0]*scaffoldLength )
            if samData[sample] != None: # some samples have no reads on a given scaffold
                for read in samData[sample]:
                    refPosition = int(read[3]) # 1-indexed
                    for p in range(len(read[9])): ## temp. ok. if the length of read is 150, then this will iterate from 0 to 149
                        if (refPosition+p) <= scaffoldLength: # if the position doesn't go beyond the length of the scaffold
                            coverageMatrix[sample][refPosition+p-1] += 1 # the minus one is for proper indexing
                            if alignmentMatrix[sample][refPosition+p-1] == "NN": # if we haven't assigned a genotype yet, then genotype it

                                # if in vcf, genotype accordingly
                                if str(refPosition+p) in snps:
                                    snpData = snps[str(refPosition+p)]
                                    ref, alt = snpData[3:5]
                                    alleleFrequencies = snpData[2]
                                    sampleIndex = sampleVCFindices[samples[sample]]
                                    sampleDatum = snpData[sampleIndex].split(":")
                                    genotypeCall = sampleDatum[0]
                                    # first look at the allele frequencies
                                    if alleleFrequencies != ".":
                                        if alleleFrequencies[3] > 0: # if 4 alleles, don't even analyze the locus
                                            genotypeCall = "N/N" # only a placeholder, because this whole locus gets thrown out below (likely paralog)
                                        elif alleleFrequencies[2] > 0: # if three alleles, throw out the least common one
                                            if alleleFrequencies[0] >= alleleFrequencies[2] and alleleFrequencies[1] >= alleleFrequencies[2]: 
                                                genotypeCall = genotypeCall.replace("2", "N")
                                            elif alleleFrequencies[0] > alleleFrequencies[1] and alleleFrequencies[2] > alleleFrequencies[1]:
                                                genotypeCall = genotypeCall.replace("1", "N")
                                            elif alleleFrequencies[1] > alleleFrequencies[0] and alleleFrequencies[2] > alleleFrequencies[0]:
                                                genotypeCall = genotypeCall.replace("0", "N")
                                            else: # here, the third allele is bigger, and the first and second are tied for second. Just call NA.
                                                genotypeCall = "N/N"
                                    if genotypeCall == "0/0":
                                        genotype = ref*2
                                    elif genotypeCall == "1/1":
                                        genotype = alt.split(",")[0] *2
                                    elif genotypeCall == "0/1":
                                        genotype = ref+alt.split(",")[0]
                                    elif genotypeCall == "./.": # sometimes we get to this point when mpileup discards a read aligned to masked nucleotides
                                        genotype = "NN" # leave as "NN", for now
                                    elif genotypeCall == "0/2":
                                        genotype = ref+alt.split(",")[1]
                                    elif genotypeCall == "1/2":
                                        genotype = alt.split(",")[0] + alt.split(",")[1]
                                    elif genotypeCall == "2/2":
                                        genotype = alt.split(",")[1]*2
                                    elif genotypeCall == "0/N" or genotypeCall == "N/0":
                                        genotype = ref+"N"
                                    elif genotypeCall == "1/N" or genotypeCall == "N/1":
                                        genotype = alt.split(",")[0]+"N"
                                    elif genotypeCall == "2/N" or genotypeCall == "N/2":
                                        genotype = alt.split(",")[1]+"N"
                                    elif genotypeCall == "N/N":
                                        genotype = "NN"
                                    else:
                                        print genotypeCall
                                        print "houston we have a problem"
                                        1/0
                                    alignmentMatrix[sample][refPosition+p-1] = genotype 

                                else:
                                    alignmentMatrix[sample][refPosition+p-1] = reference[scaffold][refPosition+p-1]*2 # diploid *2


        # now that we have the coverage and genotype saved for each sample at each position on the scaffold, 
        # go through and find regions with good coverage
        newRegion = True # default
        goodRegions = [] 
        nucs = ["A","T","C","G"]
        for position in range(scaffoldLength): # this is going to be a little confusing, but position, here, starts at 0
            popCounts = [0]*numPops # default 0, but we're going to count the number of samples covered from each pop
            for sample in range(numSamples):
                if coverageMatrix[sample][position] >= minReadDepthPerSample:
                    pop = popMap[samples[sample]]
                    popCounts[pops.index(pop)] += 1
            outputIt = True # default setting
            for p in range(numPops):
                if popCounts[p] < minSamplesEachPop:
                    outputIt = False
            if outputIt == True and reference[scaffold][position] in nucs: # the ref contains Ns and lower case bases which I think represent masking
                if newRegion == True:
                    startPos = int(position)
                    newRegion = False
            else:
                if newRegion == False: # this is the end of the region of continuous good coverage
                    endPos = int(position)-1 # minus one, so that the end position includes the last covered position
                    goodRegions.append([startPos, endPos])
                    newRegion = True



        # next we need to filter the regions to include ones (i) short and (ii) well-spaced enough to ignore both intra- and inter-locus recombination
        scanPosition = 0 # position in the reference for the sliding window
        regionIndex = 0 # index, or ID, of the final region(s)
        numRegions = len(goodRegions) # simply the number of regions to search through
        filteredRegions = [] # final list of filtered regions
        while scanPosition < scaffoldLength and regionIndex < numRegions: 
            region = goodRegions[regionIndex]
            startPos = region[0]
            endPos = region[1]
            if startPos >= scanPosition: # if the region is spaced far enough apart from the previous region
                regionLength = endPos-startPos+1
                if regionLength >= minimumLocusLength: # this region is long enough to output
                    if regionLength > maximumLocusLength: # region is too long
                        newEnd = startPos + maximumLocusLength -1 # clipping the end of the region
                        ########### later, could dig deeper to maximize coverage within the region ##########
                    else:
                        newEnd = endPos

                    # as a final filter, check heterozygosity in all 168 samples to see if it looks like a paralog
                    paralog = False
                    countSnps = 0
                    for position in range(startPos, newEnd+1 ):
                        if str(position+1) in snps:
                            countSnps += 1 # counting the number of snps in the broader dataset (168 samples)
                            countCovered = 0
                            countHeteros = 0
                            for sample in range(9, len(snps[str(position+1)])):
                                sampleDatum = snps[str(position+1)][sample].split(':')
                                dp = int(sampleDatum[2])
                                if dp >= minReadDepthPerSample:
                                    countCovered += 1
                                    geno = sampleDatum[0]
                                    if geno == "0/1":
                                        countHeteros += 1
                                    elif geno == "0/2" or geno == "1/2": # other types of heterozygotes
                                        countHeteros += 1
                                    elif geno == "0/3" or geno == "1/3" or geno == "2/3": # dude these are not real snps, man. 4-alleles = paralog
                                        paralog = True
                            if countCovered >= (minSamplesEachPop*numPops): # e.g. min 10, if min 5 samples per each of 2 pops
                                if (float(countHeteros) / countCovered) > paralogThresholds[str(countCovered)]:
                                    paralog = True 
                            else: # this seems to happen when the vcf skips reads that align partially to masked reference sequence
                                pass
                    if float(countSnps) / (newEnd-startPos) >= 0.25:
                        paralog = True # if 25% of the sequence is snps in the full dataset, 168 birds (NOT stringent, might look at this later)
                    if paralog == False: # add to final list if it looks good
                        filteredRegions.append( [ startPos, newEnd ] )
                    scanPosition = newEnd + minDistanceBetweenLoci + 1 # reset scan position variable. Indented even with "for position in range", currently
            regionIndex += 1 # indented even with "region, startPos, endPos"









        # output locus lengths
        if outputType == "locusLengths":
            for locus in filteredRegions:
                print '\t'.join( map(str, [scaffold, locus[0], locus[1], locus[1]-locus[0]+1]))









        # output sequences in fasta format
        if outputType == "fasta" or outputType == "missing":
            for locus in filteredRegions:
                locusCounter += 1
                startPos = locus[0]
                endPos = locus[1]
                locusLength = endPos-startPos+1
                for pop in pops:
                    for sample in reversePopMap[pop]:
                        sampleInd = samples.index(sample)
                        seq1 = ['N']*locusLength
                        seq2 = ['N']*locusLength
                        posCounter = 0 # position in the output sequence, always starting at 0, which is separate from position in the reference
                        allMissingData = True # default setting, checking if the sample has any coverage at this locus
                        for position in range(startPos, endPos+1):
                            cov = coverageMatrix[sampleInd][position]
                            if cov >= minReadDepthPerSample:                         
                                nucs = alignmentMatrix[sampleInd][position]
                                if nucs != "NN":
                                    seq1[posCounter] = nucs[0] # assigning a random nucleotide to each sequence, because we're not phasing
                                    seq2[posCounter] = nucs[1]
                                    allMissingData = False
                                else:
                                    if outputType == "fasta":
                                        seq1[posCounter] = reference[scaffold][position] # giving non-"N" placeholder to fasta for ms-conversion
                                        seq2[posCounter] = reference[scaffold][position]
                                    elif outputType == "missingData":
                                        pass # leaving as Ns
                            else:
                                if outputType == "fasta":
                                    seq1[posCounter] = reference[scaffold][position]
                                    seq2[posCounter] = reference[scaffold][position]
                                elif outputType == "missingData":
                                    pass # if outputting for missing data analysis, then leave the N's in.
                            posCounter += 1
                        if allMissingData == False: # if at least one base is covered for this sample
                            print ">" + ' '.join([ scaffold, "Locus"+str(locusCounter), sample, "Allele1" ])
                            print ''.join(seq1)
                            print ">" + ' '.join([ scaffold, "Locus"+str(locusCounter), sample, "Allele2" ])
                            print ''.join(seq2)









        # output locfile
        # NOTE: very very very stupid. If the locfile specifies as few as 2 samples in a population, e.g. n=2 in a test locfile, it gives "Seg fault"
        if outputType == "locfile":
            if alreadyOutputLocfileHeader == False: # don't want to output the header for every scaffold: "id    n    length    pop"
                print '\t'.join(["id", "n", "length", "pop"])
                alreadyOutputLocfileHeader = True
            for locus in filteredRegions:
                # being careful to stay consistent with the order of samples in the fasta
                outputLocusNum += 1 # loc file indexing starts at 1
                startPos = locus[0]
                endPos = locus[1]
                locusLength = endPos-startPos+1
                popCounter = 0 # keeping track of which population
                for pop in pops:
                    popCounter += 1
                    n = 0 # counting the number of samples covered in this pop
                    for sample in reversePopMap[pop]:
                        sampleInd = samples.index(sample)
                        for position in range(startPos, endPos+1):
                            cov = coverageMatrix[sampleInd][position]
                            if cov >= minReadDepthPerSample:
                                n += 1 # counting the number of samples with at least 1 base (sufficiently) covered
                                break
                    print '\t'.join([ "f"+str(outputLocusNum), str(n*2), str(locusLength), str(popCounter) ]) # n*2 for diploids









        # output vcf
        if outputType == "vcf":
            for locus in filteredRegions:
                startPos = locus[0]
                endPos = locus[1]
                for position in range(startPos, endPos+1):
                    if str(position+1) in snps:
                        print '\t'.join(snps[str(position+1)])








                    
#############
# close files
#############
vcf.close()
for s in samples:
    indFile = s + ".sam"
    samDict[indFile] = open(samDirectory+indFile, "r")



    

    









