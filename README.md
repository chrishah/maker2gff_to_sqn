# maker2gff_to_sqn

For a recent paper we wanted to submit a functionally annotated draft genome to NCBI. NCBI has a script called [table2asn_GFF](ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn_GFF/) to convert gff3 files to sqn files that can then be submitted to Genbank. NCBI provide some instructions on how to use the script [here](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/). Nevertheless it was quite a process to get the original GFF file from MAKER2 converted, so this repo details what I had to do.



First, I downloaded the program to do the conversion - I got version 1.0.0 - I saw that by now they have upgraded to version 1.29.17 or whatever so no warranty that everything will work with newer versions in the same way. I put the v.1.0.0 in a directory called `backup/`.
```bash
wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn_GFF/linux64.table2asn_GFF.gz
gunzip linux64.table2asn_GFF.gz 
./linux64.table2asn_GFF -version

```

The original files
 - `maker2.nofasta.gff.gz`
 - `Brachionus_plicatilis_scaffold_min500_real.fasta.gz`
were too big to upload them to Github so I split them up into smaller pieces - see details in the `data/README.md` file.

First, I need to restor those:
```bash

for file in $(ls -1 data/x000000* | perl -ne 'chomp; @a=split(/\./); shift(@a); $suffix=join(".", @a); print "$suffix\n";' | sort -n | uniq)
do 
	echo $file
	for part in $(ls -1 data/x000000* | grep "$file")
	do 
		cat $part
	done > data/$file
done
```


Then I run NCBI's program on the GFF file that was produced by MAKER2. First produce a gff file that does not contain the fasta bits. The file `template.sbt` was produced during the preparations for the submission via the submission portal.
```bash

#gunzip the files for now
gunzip ./data/Brachionus_plicatilis_scaffold_min500_real.fasta.gz
gunzip ./data/maker2.nofasta.gff.gz

#Run conversion
./linux64.table2asn_GFF -M n -J -c w -euk -t ./data/template.sbt -gaps-min 10 -l paired-ends -i ./data/Brachionus_plicatilis_scaffold_min500_real.fasta -f ./data/maker2.nofasta.gff -o ./B_plicatilis.pass1.sqn -Z ./B_plicatilis.pass1.dr -locus-tag-prefix BpHYR1 -n Brachionus_plicatilis -taxid 10195 -V b &> output_table2asn_B_plicatilis.pass1.txt

```

The program produces a bunch of files including a 'discrepancy report' (dr) that details any problems. The file that is needed for the submission is called \*.sqn. We have quite a few errors in the discrepancy report but the NCBI told me that most of them can be ignored since we have a eukaryotic genome here. The one thing that needs to be fixed in our case are 2414 genes that appear to have introns that are shorter than 10db which they consider problematic. NCBI initially suggested to drop all these, but many of them seem to have reasonalbe hits to known proteins when translated so we decided to keep them but flag them up as potentially problematic, by adding a 'pseudo=true' (according to the guidelines: annotate with pseudo=true any genes that are 'broken' but are not thought to be pseudogenes. These are genes that do not encode the expected translation, for example because of internal stop codons. These are often caused by problems with the sequence and/or assembly.) and a note with some more details about the problem.

In the analyses that we did for the paper we identified a few candidate genes that we discussed a bit further and I wanted to check whether any of those was affected by this issue. A list of the candidate genes is in the file `data/candidates.txt`.

```bash
#extract list of genes that are affected by the short introns issue from the discrepancy report
cat B_plicatilis.pass1.dr | perl -ne 'chomp; if ($on){if ($_ =~ /scaffold/){@a=split("\t"); print "$a[-1]\n"}}if ($_ =~ /^$/){$on=0}; if ($_ =~ /SHORT/){$on=1}' | sort -n | uniq | cut -d "_" -f 2 > shorts.list.txt

#The table2asn program assigns new unique names to the genes and only those are reported in the discrepancy report (essentially it numbers them as it encounters them in the gff file), so I need create map of original names with new names of the genes
cat data/maker2.nofasta.gff | grep -P "\tgene\t" | perl -ne 'chomp; @a=split("\t"); @b=split("=", $a[-1]); $count = sprintf("%06d", $.); print "$b[-1]\t$count\n"' > map.txt

#find the original gene id for these genes
grep -f shorts.list.txt map.txt

#Do everything in one step
grep -f <(cat B_plicatilis.pass1.dr | perl -ne 'chomp; if ($on){if ($_ =~ /scaffold/){@a=split("\t"); print "$a[-1]\n"}}if ($_ =~ /^$/){$on=0}; if ($_ =~ /SHORT/){$on=1}' | sort -n | uniq | cut -d "_" -f 2) <(cat data/maker2.nofasta.gff | grep -P "\tgene\t" | perl -ne 'chomp; @a=split("\t"); @b=split("=", $a[-1]); $count = sprintf("%06d", $.); print "$b[-1]\t$count\n"')

#check which of those are 'interesting'
grep -f <(cat data/candidates.txt | sed 's/-mRNA-.*//') <(grep -f shorts.list.txt map.txt)

#Only one of the genes is affected by this -> snap_masked-scaffold103554_len21087_cov72-processed-gene-0.19   004562

#inspect and manually fix this if possible - safe as data/maker2.nofasta.curated.gff
cp data/maker2.nofasta.gff data/maker2.nofasta.curated.gff

#in the case of this gene it was possible to merge the two exons up and downstream of the problematic intron - doint this still gives a valid amino acid sequence (so no frameshifts or stop codons were introduced) which, when blasted against genbank actually received a better hit than the original protein

#Here's what I did in detail to gff file:
#The lines (showing the two relevant exons that are too close together):
#scaffold103554_len21087_cov72   maker   exon    2432    3365    .....
#scaffold103554_len21087_cov72   maker   exon    2286    2425    .....

#were replaced with - so making it one single exon:
#scaffold103554_len21087_cov72   maker   exon    2286    3365

#and the lines (showing the two relevant CDS):
#scaffold103554_len21087_cov72   maker   CDS     2432    3365    ....
#scaffold103554_len21087_cov72   maker   CDS     2286    2425    .      

#were replaced with:
#scaffold103554_len21087_cov72   maker   CDS     2286    3365    ....
```

Now, I reran the conversion with the curated gene model.
```bash
#rerun conversion to sqn
./linux64.table2asn_GFF -M n -J -c w -euk -t ./data/template.sbt -gaps-min 10 -l paired-ends -i ./data/Brachionus_plicatilis_scaffold_min500_real.fasta -f ./data/maker2.nofasta.curated.gff -o ./B_plicatilis.pass2.sqn -Z ./B_plicatilis.pass2.dr -locus-tag-prefix BpHYR1 -n Brachionus_plicatilis -taxid 10195 -V b &> output_table2asn_B_plicatilis.pass2.txt
```

Now, the dr flagged up only 2413 genes with the SHORT INTRON issue. To deal with those as described above I have written a script, that will add the pseudo=true flag and a note.
```bash
## The remaining genes with the 'SHORT INTRON' issue will be flagged with a pseudo=true flag and a note that the respective genes are problematic
./scripts/add_pseudo_to_shortintrons.py B_plicatilis.pass2.dr data/maker2.nofasta.curated.gff > fixed_shorts.gff
```

Now, I transfer the functional annotations, GO terms, evidence codes, and enzyme codes to the gff file from the annotation table that we have produced via Blast2GO. Wrote a script to do that.
```bash
## transfer functional annotatino from blast2GO table to the gff file
./scripts/transfer_annotations_from_blast2go_go_table_to_maker2gff.py <(zcat data/blast2go_go_table_20180407_1046.txt.gz) fixed_shorts.gff > annotated.gff
```

Rename
```bash
cat annotated.gff | perl -ne 'chomp; if ($_ =~ /\tgene\t/){$count++; @a=split("\t"); @b=split(";", $a[-1]); $b[0] =~ s/ID=//; $old=$b[0]; $new="BpHYR1_".sprintf("%06d",$count); $_ =~ s/$old/$new/g; $_.=";Note=original geneID (Maker2) $old"; print "$_\n"}else{if ($old){$_ =~ s/$old/$new/g;} print "$_\n"}' > annotated.renamed.gff
```

As a last step we run NCBI's conversion program again.
```bash
#rerun conversion to sqn
./linux64.table2asn_GFF -M n -J -c w -euk -t ./data/template.sbt -gaps-min 10 -l paired-ends -i ./data/Brachionus_plicatilis_scaffold_min500_real.fasta -f annotated.renamed.gff -o ./B_plicatilis.final.sqn -Z ./B_plicatilis.final.dr -locus-tag-prefix BpHYR1 -n Brachionus_plicatilis -taxid 10195 -V b &> output_table2asn_B_plicatilis.final.txt
```

I move the final files to a backup directory to keep them and remove the rest before committing to Github.
```bash
rm *.pass1.*
rm *.pass2.*
rm map.txt shorts.list.txt
rm fixed_shorts.gff
rm annotated.gff
rm data/maker2.nofasta.curated.gff
rm data/maker2.nofasta.gff
rm data/Brachionus_plicatilis_scaffold_min500_real.fasta

mkdir backup
mv annotated.renamed.gff *B_plicatilis.final.* backup/

cd backup
gzip *.*
#split up any file that is larger than 5 MB into chunks so that each is ~1MB in size
for file in $(ls -l --block-size=MB | perl -ne 'chomp; @a=split(" "); $a[4] =~ s/MB//; if ($a[4] >= 5){print "$a[-1],$a[4]\n"}')
do
	input=$(echo "$file" | cut -d "," -f 1)
	prefix=$(echo -e "$input" | sed 's/.gz$//')
	gunzip $input
	input=$prefix
	number=$(echo $file | cut -d "," -f 2)

	split --number=$number --numeric-suffixes=0 --suffix-length=8 --additional-suffix=.$prefix $input

	#compress all files
	for f in $(ls -1 x0* | grep "$prefix"); do gzip -v $f; done

	#delete original file
	rm $input

done

```
