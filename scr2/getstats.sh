#!/bin/bash

source ./config

ind=$1

genome=$(find "$(realpath $basedir/rep)" -name "*.fa")
 
ref=$(echo $genome | rev | cut -d"/" -f1 | rev | cut -d"." -f1)

bam_file="$mapdir/$ind.$ref.bam"

cd "$gvcfdir"

filt_gvcf_file="$gvcfdir/$ind.$ref.bcftools.filt.gvcf.gz"

filt_gvcf_plink_file="$gvcfdir/$ind.$ref.bcftools.filt.plink.gvcf.gz"

## replace roman chromosome names with arabic 
bcftools annotate --rename-chrs "$basedir/rep/chrmap.txt" -o $filt_gvcf_plink_file -Oz $filt_gvcf_file

rm $filt_gvcf_file $filt_gvcf_file".csi"

bcftools index --threads 2 $filt_gvcf_plink_file

## get stats missing 
## plink --vcf $filt_gvcf_plink_file --missing --chr-set 16 no-xy --allow-no-sex -out $ind 
## rm $ind.lmiss $ind.nosex $ind.log
## mv $ind.imiss "$missdir"

## get stats alt/het/homo/missing
bcftools view "$filt_gvcf_plink_file" | awk -v ind="$ind" '{
    if ($0 ~ /1\/1/) ALThomo++
    else if ($0 ~ /0\/1/) ALThet++
    else if ($0 ~ /0\/0/) REFhomo++
    else if ($0 ~ /\.\/\./) missings++
} END { printf "%s\t%d\t%d\t%d\t%d\n", ind, ALThomo, ALThet, REFhomo, missings }' > "$gvcfstatdir/$ind.geno.txt"
