
To be able to put the larger files on Github I have split them up into smaller chunks.

For the GFF file:
```bash

#First I removed the FASTA part from the original maker2 gff file (the original file is not part of the repo)
zcat data/maker2.all.gff.gz | grep -e "^##" -e "^scaffold" | grep -v "^##FASTA" > data/maker2.nofasta.gff

#split up the gff file into smaller files with 100.000 lines each
input=maker2.nofasta.gff.gz
prefix=maker2.nofasta.gff
lines=100000

split --lines=$lines --numeric-suffixes=0 --suffix-length=8 --additional-suffix=.$prefix <(zcat $input)

#compress all files
for file in $(ls -1 x0* | grep "$prefix"); do gzip -v $file; done

#delete original file
rm $input
```

For the FASTA file:
```bash
#split up the FASTA file into smaller files with 1.000 lines each
input=Brachionus_plicatilis_scaffold_min500_real.fasta.gz
prefix=Brachionus_plicatilis_scaffold_min500_real.fasta
lines=1000

split --lines=$lines --numeric-suffixes=0 --suffix-length=8 --additional-suffix=.$prefix <(zcat $input)

#compress all files
for file in $(ls -1 x0* | grep "$prefix"); do gzip -v $file; done

#delete original file
rm $input
```

