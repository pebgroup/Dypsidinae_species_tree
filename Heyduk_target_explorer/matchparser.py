#!/usr/local/bin/python3
import re

exons = {}
products = {}

with open("bwamap.intersect.CDS.txt", "r") as matches:
	for line in matches:
		line = line.strip()
		LINE = line.split("\t")
		annot = {}
		for e in LINE[14].split(";"):
			e = e.split("=")
			annot[e[0]] = e[1]
		# list gene for exons
		if LINE[3] in exons.keys():
			if annot['gene'] not in exons[LINE[3]]:
				exons[LINE[3]].append(annot['gene']) 
		else:
			exons[LINE[3]] = [annot['gene']]
		# build dictionary of gene products
		prod = re.sub(r" isoform X[0-9]","",annot['product'])
		prod = re.sub(r"%2C",",",prod)
		if annot['gene'] not in products.keys():
			products[annot['gene']] = [prod]
		else:
			if prod not in products[annot['gene']]:
				products[annot['gene']].append(prod)
			
			
#print(exons)

genes = {}

with open("exons_vs_NCBIgenes.csv","w") as outfile:
	for i,j in exons.items():
		#if len(j) > 1:
			#print(i+" : "+str(len(j)))
	# 		Elaeis_136_5 : 2
	# 		Elaeis_136_3 : 2
	# 		Elaeis_357_6 : 2
			print(i+";"+",".join(j), file=outfile)
			gene = i.split("_")[1]
			if gene in genes.keys():
				for k in j:
					if k not in genes[gene]:
						genes[gene].append(k)
			else:
				genes[gene] = j

with open("genes_vs_NCBIgenes2.csv","w") as outfile:
	for i,j in genes.items():
		for g in j:
			for p in products[g]:
				print("Elaeis_"+i+";"+str(len(j))+";"+g+";"+str(len(products[g]))+";"+p, file=outfile)

