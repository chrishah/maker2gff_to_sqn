#! /usr/bin/python

import sys
## The script renames product names according to a tsv file
## usage: script.py new_names.tsv original.gff > new.gff

newnames=sys.argv[1]


fh = open(newnames)

names={}

for line in fh:
	temp = line.strip().split("\t")
	names[temp[0]] = {'old':temp[1], 'new':temp[2]}

fh.close()

#for gene in sorted(names):
#	print gene,names[gene]

fh=open(sys.argv[2])

for line in fh:
	line = line.strip()
	if len(line.strip().split("\t")) == 9:
		if line.strip().split("\t")[2] == 'gene':
			found=0
			temp = line.strip().split("\t")
			features = temp[-1].split(";")
			for feat in features:
				if feat.startswith("ID="):
					gene_ID = feat.replace("ID=","")
					break

#			print gene_ID
			
			if gene_ID in names: #searching for the gene by name
				found=1
#				print "FOUND"
				for i in reversed(range(len(features))):
					if features[i].startswith("description="):
#						print "Replacing %s with %s" %(features[i].replace("description=",""), names[gene_ID]['new'])
						if not names[gene_ID]['new'] == '-':
							features[i] = "description=%s" %names[gene_ID]['new']
						else:
							del(features[i])
						break

				feats = ";".join(features)
				temp[-1] = ";".join(features)
				line = "\t".join(temp)

		if line.strip().split("\t")[2] == 'CDS' and found:
			temp = line.strip().split("\t")
                        features = temp[-1].split(";")
			for i in reversed(range(len(features))):
				if features[i].startswith("product="):
					if not names[gene_ID]['new'] == '-':
						features[i] = "product=%s" %names[gene_ID]['new']
					else:
						del(features[i])
					break
			feats = ";".join(features)
                        temp[-1] = ";".join(features)
                        line = "\t".join(temp)

	print line
