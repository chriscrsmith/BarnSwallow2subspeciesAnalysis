

# convert fasta to ms format

import sys, os

sys.stderr.write("\nWARNING: data given to this script must not contain Ns. If any sample has N the site will be completely skipped over.\n\n")




# remove temp files that get appended to
bashCommand = "rm tempFile.ms tempFile.stderr\n"
sys.stderr.write(bashCommand)
os.system(bashCommand)


def processSeqs(sequences):
    with open("tempFile.fa", "w") as tempFile:
        for s in sequences:
            tempFile.write(s + "\n")
    bashCommand = "perl ~/Software/msABC20120315_chrisVersion/fas2ms.pl fas=tempFile.fa 1>> tempFile.ms 2>> tempFile.stderr\n"
    os.system(bashCommand)

    


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
                processSeqs(data)
                currentLocus = locusID # reset variables
                data = [line.strip()]
            else:
                data.append(line.strip())
processSeqs(data) # final locus       
