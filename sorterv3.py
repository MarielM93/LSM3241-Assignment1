from __future__ import division
import csv

# this has to be done beforehand
# cut -f 3 probeMappingRMA.csv > cutRMA.csv
# cut -f 3 probeMappingMAS.csv > cutMAS.csv

# merging all the DEGs captured from RMA and MAS5
def merger(dict):
	for object in dict:
		if object not in deg_merged:
			deg_merged.append(object)

with open('cutRMA.csv', 'rb') as f:
	reader = csv.reader(f)
	rma_raw = list(reader)

rma = []
for gene in rma_raw:
	if gene not in rma:
		rma.append(gene)

print "no. of RMA hits: %d" % len(rma) 
	
with open('cutMAS.csv', 'rb') as f:
	reader = csv.reader(f)
	mas_raw = list(reader)

mas = []
for gene in mas_raw:
	if gene not in mas:
		mas.append(gene)	
print "no. of MAS5 hits: %d" % len(mas)
	
deg_merged = []
merger(mas)
merger(rma)

print "Here's all the DEGs after merging MAS5 and RMA data:"
print "No. of DEGs: %d" % (len(deg_merged))
print "Merged DEG list:"
print deg_merged

# finding all the DEGs not captured by RMA
uncaptrma = []
for gene in mas:
	if gene not in rma:
		if gene not in uncaptrma:
			uncaptrma.append(gene)

print "Here's what RMA normalization missed out but MAS5 captured:"
print "No. of DEGs: %d" % (len(uncaptrma))
print "DEG list:"
print uncaptrma

# finding all the DEGs not captured by MAS5
uncaptmas = []
for gene in rma:
	if gene not in mas:
		if gene not in uncaptmas:
			uncaptmas.append(gene)
			
print "Here's what MAS5 normalization missed out but RMA captured:"
print "No. of DEGs: %d" % (len(uncaptmas))
print "DEG list:"
print uncaptmas

# save this into a file:
# python sorterv3.py > normcomp.txt