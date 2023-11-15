# creating file of pairwise Fst values between all of my populations for IBE v IBD analysis in R

import os
from statistics import mean

vcf = "WGS_subset.vcf"

# initialize existing populations and samples
for line in open(vcf, "r"):
	if line.startswith("#CHROM"):
		line = line.split()
		samples = line[9:]
		pops = []
		for x in samples:
			if x[:-1] not in pops:
				pops.append(x[:-1])


# Do Fst for each existing paur
for i in range(len(pops)):
    # Compare the current element with all other elements: i+1 allows for no repeat comparisons
    for j in range(i + 1, len(pops)):
    	pop1 = pops[i]
    	pop2 = pops[j]
    	print(f"comparing {pop1} and {pop2}")
    	filename_sample1 = "temp_keep1.txt"
    	filename_sample2 = "temp_keep2.txt"
    	out1 = open(filename_sample1, "w")
    	out2 = open(filename_sample2, "w")
    	for x in samples:
    		if x[:-1] == pop1:
    			out1.write(x + "\n")
    		if x[:-1] == pop2:
    			out2.write(x + "\n")
    	out1.close()
    	out2.close()

    	fst_outprefix = f"{pop1}_and_{pop2}_fst"
    	command=f"/Users/whitneybrannen/vcftools_0.1.13//bin/vcftools --vcf {vcf} --weir-fst-pop {filename_sample1} --weir-fst-pop {filename_sample2} --out Fst_outputs/{fst_outprefix}"
    	os.system(command)


pairwise_fst = open("all_pairwise_fst.txt", "w")
pairwise_fst.write("Pop1\tPop2\tMeanFst\n")

for file in os.listdir("Fst_outputs"):
	pop1 = file.split("_")[0]
	pop2 = file.split("_")[2]
	fst_all = []
	for line in open(f"Fst_outputs/{file}","r"):

		fst = line.split()[2]
		if fst != 'nan' and fst != 'WEIR_AND_COCKERHAM_FST':
			if float(fst) >= 0:
				fst_all.append(float(fst))
	mean_fst = mean(fst_all)
	pairwise_fst.write(f"{pop1}\t{pop2}\t{mean_fst}\n")

pairwise_fst.close()




#os.system('/Users/whitneybrannen/vcftools_0.1.13//bin/vcftools')