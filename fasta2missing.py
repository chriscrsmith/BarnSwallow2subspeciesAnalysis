

# convert fasta to missing Data

import sys, os

sys.stderr.write("\nWARNING: data given to this script must contain Ns, unlike the fasta used for converting to ms format.\n\n")


outName = sys.argv[2] # all files with this prefix will be REMOVED from the Missing/ directory, and replaced with new files


# remove temp files that get appended to
bashCommand = "mkdir Missing"
sys.stderr.write(bashCommand)
os.system(bashCommand)
bashCommand = "rm Missing/" + outName + "*" + "\n"
sys.stderr.write(bashCommand)
os.system(bashCommand)

# open the primary missing data file
missingFile = open( ("Missing/"+outName+".missing"), "w")



def processSeqs(sequences, locusID, outFile):
    with open("Missing/TEMP"+outName+".fa", "w") as tempFile:
        for s in sequences:
            tempFile.write(s + "\n")
    bashCommand = "perl ~/Software/msABC20120315_chrisVersion/missing_ms_report.pl in=" + ("Missing/TEMP"+outName+".fa")+" > Missing/" + (outName+"_"+locusID+".txt") + " 2>> Missing/tempFile.stderr\n"
    os.system(bashCommand)
    # write to primary missing data file
    outFile.write("Missing/"+outName+"_"+locusID+".txt\n")
  


with open(sys.argv[1], "r") as infile:
    newline = infile.readline().strip()
    currentLocus = newline.split()[1]
    data = [newline]
    for line in infile:
        if line[0] != ">": # sequence
            data.append(line.strip())
        else: # header
            locusID = line.strip().split()[1]
            if locusID != currentLocus: # encountered new locus, process the seqs
                processSeqs(data, currentLocus, missingFile)
                currentLocus = locusID # reset variables
                data = [line.strip()]
            else:
                data.append(line.strip())
processSeqs(data, currentLocus, missingFile) # final locus       
missingFile.close()
