#!/usr/bin/python

import sys

#The script transfers annotations from a blast2go annotation table to a maker gff file.
#usage: script.py blast2go_table.txt maker2.gff > maker2_annotated.gff

table=sys.argv[1]
gff=sys.argv[2]

#print table
#print gff

dictionary = {}

fh=open(table)
for line in fh.readlines()[1:]:
#	print line.strip()
	l = line.strip().split("\t")
	ID=l[0]
	dictionary[ID] = {'product':'', 'GO':'', 'ec_number':'', 'IPR':'', 'Dbxref':''}
	dbxref_list = []
	if '--NA--' in l[2] or 'LOC' in l[2] or 'hypothetical protein' in l[2]: #'--NA--' seems to be used when there was no significant blast hit; Some descriptions indicate similarity to uncharacterized proteins with an ID starting with LOC.... (uncharacterized proteins from a variety of organisms in the database) or to some hypothetical proteins from specific organisms. These cases are ignored by the script so the gene in question won't receive a functional annotation in the gff file 
		pass
	else:
		dictionary[ID]['product'] = l[2]
	if len(l) >= 9:
		GOs = l[8].split("; ")
		for i in range(len(GOs)):
			GOs[i] = ":".join(GOs[i].split(":")[1:])
		
		dictionary[ID]['GO'] = ",".join(GOs)
	if len(l) >= 11:
		dictionary[ID]['ec_number'] = ":".join(l[10].split(":")[1:])
	if len(l) >= 13:
		if not l[12] == 'no IPS match':
			IPRs = l[12].split("; ")
			final_IPRs = []
#			print IPRs
			for i in reversed(range(len(IPRs))):
				if IPRs[i].startswith('IPR'):
					final_IPRs.append("InterPro:%s" %IPRs[i].split(" ")[0])
				else:
					if not 'COILS' in IPRs[i]:
						db=IPRs[i].split(" ")[1].replace("(","").replace(")","")
						final_IPRs.append("%s:%s" %(db,IPRs[i].split(" ")[0]))


#			print final_IPRs
			dictionary[ID]['IPR'] = ",".join(final_IPRs)

#	print "PRELIMINARY: %s\t%s" %(ID, dictionary[ID])
	if dictionary[ID]['GO']:
		dbxref_list.append(dictionary[ID]['GO'])

	if dictionary[ID]['IPR']:
                dbxref_list.append(dictionary[ID]['IPR'])
	
#	print dbxref_list
	if dbxref_list:
		dictionary[ID]['Dbxref'] = ",".join(dbxref_list)
	
#	print "%s\t%s "%(ID,dictionary[ID])
#	1,3,9,11,13

for line in open(gff):
	line=line.strip()
	if "\tgene\t" in line:
		ID = line.split("\t")[8].split(";")[0].split("=")[1]
		found=0
		for ids in dictionary.keys():
			if ID == "-".join(ids.split("-")[:-2]):
#				print "found\t%s" %dictionary[ids]
				if dictionary[ids]['product']:
					line += ";description=%s" %dictionary[ids]['product']
				if dictionary[ids]['Dbxref']:
					line += ";Dbxref=%s" %dictionary[ids]['Dbxref']
				print line
				found=1
				break
		if not found:
			print line
	elif "\tCDS\t" in line:
#		print line
		mRNA_ID = line.split("\t")[8].split(";")[0].split("=")[1].split(":")[0]
#		print mRNA_ID,dictionary[mRNA_ID]

		for qual in ['product', 'Dbxref', 'ec_number']:
			if dictionary[mRNA_ID][qual]:
				line += ";%s=%s" %(qual, dictionary[mRNA_ID][qual])
		print line
	else:
		print line
