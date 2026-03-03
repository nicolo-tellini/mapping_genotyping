### Set in variable so that it corrrsponds to your gVCF ###
in=EU-list-clean-norareline-noazerb.gvcf.gz

# Extract biallelic SNPs in core genomic regions.
# -COUNT(GT~"1")>0: retain sites where ALT allele 1 is observed in at least one genotype.
# - -m2 -M2: keep  biallelic variants.
# - -R core.bed: restrict analysis to core regions (no subtelomeric regions).
bcftools view -i 'COUNT(GT~"1")>0' -m2 -M2 -R core.bed $in -Oz -o varonly.core.vcf.gz

# Filter on missingness.
# --geno 0.05: remove SNPs with >5% missing genotypes.
# --mind 0.1: remove strains with >10% missing sites.
plink2 --vcf varonly.core.vcf.gz --mind 0.1 --geno 0.05 --make-bed --recode vcf --out varonly.core.geno005mind01

# Compute minor allele frequency (MAF)
bcftools +fill-tags varonly.core.vcf.gz -- -t MAF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%MAF\n' > maf.txt

# Convert VCF to PLINK binary format.
# --set-all-var-ids @:# assigns variant IDs as CHR:POS.
plink2 --vcf varonly.core.vcf.gz --allow-extra-chr --make-bed --out step0 --set-all-var-ids @:#


# Retain variants with MAF ≥ 0.01 (≈ MAC ≥ 26 in 1293 diploid strains = 2586 chromosomes).
plink2 --bfile step1 --maf 0.01 --make-bed --out step2

# Compute pairwise LD (r²) within 10 kb windows.
# --ld-window 99999 removes SNP-count limit within the kb window.
plink --bfile step2 --r2 --ld-window 99999 --ld-window-kb 10 --ld-window-r2 0 --out ld

# Identify approximately independent SNPs using LD-based pruning.
# Window size: 50 SNPs, step: 5 SNPs, r² threshold: 0.2.
plink2 --bfile step2 --indep-pairwise 50 5 0.2 --out pruned

# Create LD-pruned dataset.
plink2 --bfile step2 --extract pruned.prune.in --make-bed --out step_pca

# Principal Component Analysis (first 10 PCs).
plink2 --bfile step_pca --pca 10 --out pca

# Export LD-pruned SNPs to bgzipped VCF for phylogenetic analyses.
plink2 --bfile step_pca --recode vcf bgz --out iqtree_snps
bcftools index iqtree_snps.vcf.gz

# Export non–LD-pruned SNP set (after missingness filtering only).
# Used for window-based PCA to retain local LD structure and rare/private variation.
plink2 --bfile step1 --export vcf bgz --out step1vcf

# Convert VCF to PHYLIP format for IQ-TREE.
vcf2phy -i iqtree_snps.vcf.gz --output-prefix iqtree_snps

# Maximum-likelihood phylogeny with ascertainment bias correction for SNP-only data.
# -m MFP+ASC: model selection with Lewis correction.
# -B 1000: ultrafast bootstrap.
# -alrt 1000: SH-aLRT branch test.
iqtree2 -s iqtree_snps.min4.phy -m MFP+ASC -B 1000 -alrt 1000 -T 44 --prefix iqtree

# Alternative alignment (e.g., variable sites only).
iqtree2 -s iqtree.varsites.phy -m MFP+ASC -B 1000 -alrt 1000 -T 44 \
--prefix iqtree2 > iqtree2.log 2> iqtree2.err &

# Window-based PCA (winpca) on non–LD-pruned dataset.
mkdir -p winpca
cp step1vcf.vcf.gz winpca
cd winpca
gzip -d step1vcf.vcf.gz

cp ../../rep/chrs.txt .

# chrs.txt format: CHR START END LABEL
# Sliding window PCA: 20 kb window, 500 bp increment.
while read -r line
do
  chr=$(echo $line | cut -d" " -f1)
  start=$(echo $line | cut -d" " -f2)
  end=$(echo $line | cut -d" " -f3)
  label=$(echo $line | cut -d" " -f4)
  echo $chr:$start-$end $label
  /home/ntellini/tools/winpca/winpca pca step1vcf.vcf \
    $chr:$start-$end $label --window_size 20000 --increment 500
done < chrs.txt
conda deactivate

# pixy requires invariant sites; regenerate core-restricted dataset without ALT filtering.
bcftools view -R core.bed $in -Oz -o pixy.core.gvcf.gz
bcftools index pixy.core.gvcf.gz

# Define populations.
# Each strain is treated as its own population (enables all-vs-all dXY).
# For FST, biologically meaningful multi-sample groups must be defined.
bcftools query -l pixy.core.gvcf.gz | awk '{print $1"\t"$1}' > populations.txt

pixy --vcf pixy.core.gvcf.gz \
     --populations populations.txt \
     --stats pi watterson_theta tajima_d dxy \
     --window_size 2000 \
     --n_cores 10 \
     --output_folder pixy_out \
     --output_prefix pixy_2kb
