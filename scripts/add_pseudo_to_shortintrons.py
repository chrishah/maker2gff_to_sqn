#! /usr/bin/python

import sys
## The script reads a discrepancy report from table2asn_GFF, identifies the genes that had a 'SHORT INTRON' issue and adds pseudo=true and a note to the respective gene to indicate that it's problematic to the gene feature in the gff file
##usage: script.py discrepancy_report.dr original.gff > new.gff

dr=sys.argv[1]


short_intron_gene_ids = []
short_intron_gene_numbers = []
on=0
count=0

fh = open(dr)

for line in fh:
	if on:
		if len(line.strip()) == 0:
			break
		else:
#			print line.strip()
			count = int(line.strip().split("\t")[3].split("_")[1])
			short_intron_gene_numbers.append(count)
			ID = line.strip().split("\t")[3]
			short_intron_gene_ids.append(ID)
	

        if 'SHORT' in line:
                count+=1
        if count == 2:
                on=1

#print len(short_intron_gene_ids),short_intron_gene_ids
#print len(short_intron_gene_numbers),short_intron_gene_numbers

fh.close()

fh=open(sys.argv[2])
genecount=1

for line in fh:
	line = line.strip()
	if len(line.strip().split("\t")) == 9:
		if line.strip().split("\t")[2] == 'gene':
			found = 0
			gene_ID = line.strip().split("\t")[-1].split("=")[-1]
			
			if gene_ID in short_intron_gene_ids: #searching for the gene by name
				found = 1
			else:
				if genecount in short_intron_gene_numbers: #searching for the gene by number
					found = 1
			
			if found:
				line += ";pseudo=true;Note=predicted gene contains short introns (<10 bp) - possibly nonfunctional due to frameshift or problems with the sequence, assembly and/or prediction" 
				

			genecount += 1
	print line
