
# example:
sample(c(1,2), 10, replace = T) # 1 for hetero, 2 for not. If diploid, then 1/2 chance of being hetero by chance.

# find a good threshold for each different number of samples covered
proportionYouWant = 0.001
numReps = 100000
propsOutput = c()
for (s in 10:1547)
{
	numSamples = s
	collection = rep(NA, numReps)
	
	# simulate numReps times
	for (i in 1:numReps)
	{
		sim = sample(c(1,2), numSamples, replace = T)
		#sim = sample(c(1,2,2,2,2,2,2,2,2,2), numSamples, replace = T)
		count = length(which(sim==1))
		collection[i] = count
	}
	
	# find the point where the proportion you want is expected to be hetero
	searching = T
	currentProp = 0.99
	while (searching == T)
	{
		prop = length(collection[collection > numSamples*currentProp]) / numReps
		if (prop >= proportionYouWant)
		{
			print(c(s, currentProp))
			propsOutput = c(propsOutput, currentProp)
			searching = F
		}
		currentProp = currentProp - 0.01
	}
}
numCovered = seq(10,1547)
df = cbind(numCovered, propsOutput)
df

write.table(df, "paralogThresholds_0.001_1500samples.txt", sep = "\t", row.names = F, col.names = F, quote = F)












