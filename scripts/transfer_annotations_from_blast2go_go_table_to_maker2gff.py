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

#Weed out problematic product descriptions

#These cases are ignored by the script so the gene in question won't receive a functional annotation in the gff file 
#'--NA--' seems to be used when there was no significant blast hit;
#Some descriptions indicate similarity to uncharacterized proteins with an ID starting with LOC.... (uncharacterized proteins from a variety of organisms in the database)
#hypothetical proteins from specific organisms.
#descriptions containing 'unnamed'
#Low quality protein causes a fatal error
#anything containing uncharacterized or uncharacterised I use 'ncharacteri' as the common pattern
#description containing unkown - 'nknown'
	exceptions = ['--NA--', 'LOC', 'hypothetical protein', 'Hypothetical protein', 'unnamed', 'LOW QUALITY PROTEIN', 'ow quality protein', 'ncharacteri','nknown']
	if l[2].isdigit(): #if the description contains only of numbers
		l[2] = '' #reset the description to nothing
	else:
		for case in exceptions: #if one of the exceptions is true
			if case in l[2]:
				l[2] = '' #reset the description to nothing
				break
	
	if l[2]: #if there is something left
		if 'partial' in l[2]:
			l[2] = l[2].replace("partial", "") #according to the annotation guidelnies the word partial is not allowed and will thus be removed, in the descriptions produced by Blast2GO 'partial' just indicates that it matched to a partial gene as far as I understand

		#check for any weird characters at the beginning of the description
		ok = 0
		while not ok: #the description should come out without any strange things at the start of this while loop
			if l[2].startswith(" "):
				l[2] = l[2].lstrip() # lstrip removes any spaces at the beginning of the string
			elif l[2].startswith('-'): # sometimes the description starts with a hyphen which causes problems - these are removed
				l[2] = l[2][:-1]
			elif l[2].startswith('['):
				l[2] = l[2][:-1] # sometimes the description starts with a '[' which causes problems - these are removed
			elif l[2].startswith(','):
				l[2] = l[2][:-1]
			else:
				ok = 1

		#check for any weird characters at the end of the description
		ok = 0
		while not ok: #the description should come out without any strange things at the end of this while loop
			if l[2].endswith(" "):
				l[2] = l[2].rstrip() # rstrip removes any spaces at the end of the string
			elif l[2].endswith('-'): # sometimes the description ends with a hyphen which causes problems - these are removed
				l[2] = l[2][:-1]
			elif l[2].endswith('['):
				l[2] = l[2][:-1] # sometimes the description ends with a '[' which causes problems - these are removed
			elif l[2].endswith(','):
				l[2] = l[2][:-1]
			else:
				ok = 1

		if l[2]: #if there's something left after the above steps
			dictionary[ID]['product'] = l[2] 

	if len(l) >= 9:
		GOs = l[8].split("; ")
		for i in range(len(GOs)):
			GOs[i] = ":".join(GOs[i].split(":")[1:])
		
		dictionary[ID]['GO'] = ",".join(GOs)
	if len(l) >= 11:
#		print l
		if l[10].startswith("EC"):
			ECs = l[10].split("; ")
			final_ECs = []
			for i in reversed(range(len(ECs))):
				if len(ECs[i].split(":")[1].split(".")) == 4: #The NCBI submission only accepts EC numbers that are in the format x.x.x.x (so 4 numbers separated by '.') -> x.x.x sometimes occurs in Blast2GO, but these will be excluded
					final_ECs.append(ECs[i].split(":")[1])
			dictionary[ID]['ec_number'] = ",".join(final_ECs)
	if len(l) >= 13:
		if not l[12] == 'no IPS match':
			IPRs = l[12].split("; ")
			final_IPRs = []
#			print IPRs
			for i in reversed(range(len(IPRs))):
				if IPRs[i].startswith('IPR'):
					final_IPRs.append("InterPro:%s" %IPRs[i].split(" ")[0])
				else:
					if 'PIRSF' in IPRs[i]:
						final_IPRs.append("PIR:%s" %IPRs[i].split(" ")[0])
					else:
						pass # other databases like GENE3D,HAMAP,PANTHER,PRINTS,PRODOM,PROSITE_PATTERNS,PROSITE_PROFILES,SMART,SUPERFAMILY seem not to be supported by NCBI and throw an error when trying to submit the annotation 
#					if not 'COILS' in IPRs[i]:
#						db=IPRs[i].split(" ")[1].replace("(","").replace(")","")
#						final_IPRs.append("%s:%s" %(db,IPRs[i].split(" ")[0]))


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
		pseudo = 0
		ID = line.split("\t")[8].split(";")[0].split("=")[1]
		if 'pseudo=true' in line.strip():
			pseudo = 1
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
		if not pseudo:
			mRNA_ID = line.split("\t")[8].split(";")[0].split("=")[1].split(":")[0]
#			print mRNA_ID,dictionary[mRNA_ID]

			for qual in ['product', 'Dbxref']: #add product name and Evidence codes from GO and IPR to the CDS if any exist
				if dictionary[mRNA_ID][qual]:
					line += ";%s=%s" %(qual, dictionary[mRNA_ID][qual])
			if dictionary[mRNA_ID]['product'] and dictionary[mRNA_ID]['ec_number']: #only add the EC number if the protein has a valid product name
				line += ";ec_number=%s" %dictionary[mRNA_ID]['ec_number']

		print line
	else:
		print line
