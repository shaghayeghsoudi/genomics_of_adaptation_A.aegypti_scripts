#/share/lanzarolab/users/shaghayegh/data/A.aegypti/good_ind_Q30.DP5.caling_1_maf0.03/bed.files
awk '{print $1, $4+1, $5, $9}' Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.good.pos.GENES.gff3 | sed 's/ /\t/g' > genes_good_bed
awk -F__ '{print $1, $2-1, $2, $0}' chrom_pos.g | sed 's/ /\t/g' > snps.bed
intersectBed -a snps.bed -b genes_good.bed -wb

awk '
{
    str = $1$2$3$4; 
}
FNR == NR {
    arr[str] = $NF;
}
FNR != NR {
    gene_name = arr[str] ? arr[str] : "NA";
    print $0, gene_name;
}' intersections_genes_SNPs_bedtools.txt snps.bed > snps_with_assigned_genes



