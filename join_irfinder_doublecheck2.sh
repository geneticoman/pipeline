

run=4 ## 1 2 3

project=GSE141633
prefix=SRR
irtype=dir ## nondir for unstranded rnaseq lib; dir for stranded rnaseq lib; irfinder automaticaly detects lib type; if rnaseq is stranded, it will also output dir output file, in addition to the default nondir output file by default; so for unstranded lib, you will only see nondir output file, but for stranded lib, you will see both nondir and dir output file, and you should use the dir output file for further analysis; note the automate detection of lib type may not be accurate, so sometimes you will see some samples have both nondir and dir output file, some only have nondir output file ... in this case, use nondir output file which is available for *all samples* for further analysis! Otherwise, if dir output file is available for *all samples*, use dir output file.

## ------------------------------------------------

projPath=/gpfs/data/pitroda-lab/GEO/GSE141633/IRFinder_out

cd $projPath

## ------------------------------------------------

## create symbolic links to irfinder output folders 
if [ $run -eq 1 ]; then 

	echo -e "run = $run "
	
	mkdir irfinder 
	cd irfinder

	for file in `ls -d ${projPath}/${prefix}*/`; do echo $file; ln -s $file;done

	ls -al ${prefix}* | wc -l
	## 91 samples 

fi 

## ------------------------------------------------

cd $projPath
echo ${projPath}/${prefix}*/IRFinder-IR-${irtype}.txt
samples=`ls -al ${projPath}/${prefix}*/IRFinder-IR-${irtype}.txt | sed 's/\//\t/g' | cut -f 8 | sort `
echo $samples | awk '{print NF}' ## 91 samples 

## ------------------------------------------------

## join all IRFinder output files together ...
out="${project}.rnaseq.irfinder.tsv"

## ------------------------------------------------

## concat irfinder output files into one matrix
if [ $run -eq 2 ]; then 

	echo -e "run = $run "

	# cp -p ../CRC_multisample/irfinder.header.0 .

	cat /gpfs/data/pitroda-lab/ReferenceData/IRFinder/irfinder.header.0 > $out 

	for sample in $samples; do 
		in=""
		in="${projPath}/$sample/IRFinder-IR-${irtype}.txt"
		echo $sample 
		echo $in

		awk -F"\t" -v s=$sample '{print s"\t"$0}' $in >> $out

		## produce one filtered IR bed file per sample...
		awk -F"\t" '$4 ~/clean/ &&  $9 + $19 >= 6 && $21 ~/-|LowCov|NonUniformIntronCover/' $in > $in.clean_keepNonUniformIntronCoverOrLowCovDP6.txt 
	done 

	echo -e "GeneSymbol\tEnsemblID\tStaticWarning" > 1.txt 
	awk -F"\t" 'NR>1 {print $5}' $out | sed 's/\//\t/g' >> 1.txt 
	paste $out 1.txt > 2.txt 

	mv 2.txt $out 

	rm 1.txt 

fi 

## ------------------------------------------------

## filter irfinder using different filters ....
if [ $run -eq 3 ]; then 

	echo -e "run = $run "

	## select clean IRs 
	awk -F"\t" 'NR==1 || ( NR>1 && $25=="clean" && $22 == "-" )' $out  > ${project}.rnaseq.irfinder.clean.tsv

	## select clean or non uniform cov IRs
	awk -F"\t" 'NR==1 || ( NR>1 && $25=="clean" && ( $22 == "-" || $22 == "NonUniformIntronCover"  ) ) ' $out  > ${project}.rnaseq.irfinder.clean_keepNonUniformIntronCover.tsv

	## add intron depth>= 1 filter
	awk -F"\t" 'NR==1 || ( NR>1 && $25=="clean" && $22 == "-" && $10 >=1 )' $out  > ${project}.rnaseq.irfinder.clean.cov1.tsv
	awk -F"\t" 'NR==1 || ( NR>1 && $25=="clean" && ( $22 == "-" || $22 == "NonUniformIntronCover"  ) && $10 >=1 ) ' $out  > ${project}.rnaseq.irfinder.clean_keepNonUniformIntronCover.cov1.tsv

	## add intron depth>= 2 filter
	awk -F"\t" 'NR==1 || ( NR>1 && $25=="clean" && $22 == "-" && $10 >=2 )' $out  > ${project}.rnaseq.irfinder.clean.cov2.tsv
	awk -F"\t" 'NR==1 || ( NR>1 && $25=="clean" && ( $22 == "-" || $22 == "NonUniformIntronCover"  ) && $10 >=2 ) ' $out  > ${project}.rnaseq.irfinder.clean_keepNonUniformIntronCover.cov2.tsv

	## add intron depth>= 3 filter
	awk -F"\t" 'NR==1 || ( NR>1 && $25=="clean" && $22 == "-" && $10 >=3 )' $out  > ${project}.rnaseq.irfinder.clean.cov3.tsv
	awk -F"\t" 'NR==1 || ( NR>1 && $25=="clean" && ( $22 == "-" || $22 == "NonUniformIntronCover"  ) && $10 >=3 ) ' $out  > ${project}.rnaseq.irfinder.clean_keepNonUniformIntronCover.cov3.tsv

fi 

## ------------------------------------------------

## manually testing one file to see if my filter setting is right ....
if [ 0 -eq 1 ]; then 

	sample=SRR8618319
	file=${sample}.IRFinder-IR-${irtype}.txt
	ln -s ../${project}_samples/${sample}/intron_retention/irfinder/${sample}/IRFinder-IR-${irtype}.txt $file 
	awk -F"\t" '$4 ~/clean/ &&  $9 + $19 >= 6 && $21 ~/-|LowCov|NonUniformIntronCover/' $file > $file.clean_keepNonUniformIntronCoverOrLowCovDP6.txt
	wc -l $file.clean_keepNonUniformIntronCoverOrLowCovDP6.txt
	## 82086

	## to make sure when I filter from the joint matrix I get the same result as per-sample filter ....
	awk -F"\t" -v s=$sample '$1==s && $5 ~/clean/ && $10 + $20 >= 6 && $22 ~/-|LowCov|NonUniformIntronCover/' ${project}.rnaseq.irfinder.tsv | wc -l
	## 82086

	## yay! it is consistent! Proceed to the next step
fi 

## ------------------------------------------------

## filter irfinder with the final filter  ....
if [ $run -eq 4 ]; then 

	echo -e "run = $run "

	awk -F"\t" 'NR==1 || (NR>1 && $5 ~/clean/ && $10 + $20 >= 6 && $22 ~/-|LowCov|NonUniformIntronCover/)' ${project}.rnaseq.irfinder.tsv > ${project}.irfinder.clean_keepNonUniformIntronCoverOrLowCovDP6.tsv

	awk -F"\t" 'NR==1 || (NR>1 && $5 ~/clean/ && $10 + $20 >= 10 && $22 ~/-/)' ${project}.rnaseq.irfinder.tsv > ${project}.irfinder.clean_DP10.tsv

	## Sean wants ERAP2 and IFNAR1....
	#awk -F"\t" 'NR == 1 || $23 == "ERAP2" || $23 == "IFNAR1"' ${project}.irfinder.clean_keepNonUniformIntronCoverOrLowCovDP6.tsv > ${project}.irfinder.clean_keepNonUniformIntronCoverOrLowCovDP6.ERAP2_IFNAR1.tsv

fi 

# ## ------------------------------------------------
