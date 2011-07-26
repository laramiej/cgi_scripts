#!/bin/tcsh
set j1=1153
set j2=3719
set chr=chr1
set beg=209934007
set end=209937079
set evDnbs=/Users/jlaramie/Documents/test_data/GS10851-1100-36-ASM/GS00392-DNA_G01/ASM/SV/evidenceJunctionDnbsBeta-GS10851-1100-36-ASM.tsvecho browser position ${chr}:${beg}-${end} 
echo browser hide allforeach j ( $j1 $j2 ) 
	echo track name='"'Junction $j'"' description='"'mate pairs in evidence of junction $j'"' visibility=squish itemRgb=on 
	awk -v j=$j -v j2=$j2 -v chr=$chr -v beg=$beg -v end=$end -v OFS="\t" -v dq='"' \		'$1!=j{next}\ $6=="L"&&$7=="-"{color="0,255,0"}\ $6=="R"&&$7=="+"{color="100,100,0"}\ 
		$6=="L"&&$7=="+"{color="255,0,0"}\ $6=="R"&&$7=="-"{color="0,0,255"}\ 
		{print $8,$9,$15,"matePair" NR,1,$7,$9,$15,color,2,"42,42",0 "," $15-42-$9}' \	$evDnbsend