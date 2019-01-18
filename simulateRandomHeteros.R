

# Assuming allele frequencies p=0.5 and q=0.5 to remain conservative. HWE expects 1/2 heterozygous.
  # example:
sample(c(1,2), 10, replace = T) # 1 for hetero, 2 for homozygous.



# I want to know what proportion of heterozygotes in my sample would indicate a likely paralog.
# Finding the proportion of heterozygotes above which there is <0.001 probability of being sampled by chance in this very basic model.
# Therefore, I would expect the false positive rate to be about 0.001, but it should filter out many of the obvious paralogs.
proportionYouWant = 0.001
numReps = 100000
propsOutput = c()
maxSamples = 1547
for (s in 10: maxSamples)
{
        numSamples = s
        collection = rep(NA, numReps)
        print(s)
        
        # simulate numReps times
        for (i in 1:numReps)
        {
                sim = sample(c(1,2), numSamples, replace = T)
                count = length(which(sim==1))
                collection[i] = count
        }
        
        sortedCollection = sort(collection)
        threshold = sortedCollection[numReps - (proportionYouWant* numReps)] / numSamples
        propsOutput = c(propsOutput, threshold)
}
numCovered = seq(10, maxSamples)
df = cbind(numCovered, propsOutput)
df

write.table(df, "paralogThresholds_0.001_1500samples.txt", sep = "\t", row.names = F, col.names = F, quote = F)






