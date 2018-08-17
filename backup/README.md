The original files can be restored with this command:
```bash
for suffix in $(ls -1 x00* | perl -ne 'chomp; @a=split(/\./); shift(@a); print join(".", @a)."\n"' |sort -n | uniq)
do 
	for file in $(ls -1 x000000* | grep "$suffix")
	do 
		cat $file
	done > $suffix
done
```

Get rid of the small files via
```bash
rm x0*
```


If you want to split larger files up again (in this case anything larger than 5MB, will be split up into as many smaller files as necessary so that each is about 1MB in size when compressed. The loop expects that the original files are gzipped.
```bash

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
