###Project_U01_SB_Cytokine###
#hogwash
#grouped by gene
#sync


library(readxl)
library(ape)
library(hogwash)

p="IDSA_severe"

# FUNCTION ----
# return vector of non-variant sites
index_of_non_var_sites <- function(mat){
  index = as.numeric(c(which(rowSums(mat) == 0), which(rowSums(mat) == ncol(mat))))
  return(index)

}

#GENO ----
#import panaroo results
#panaroo_tab <- read.table('../../../2021_02_19_merged_panaroo_data/data/gene_presence_absence.Rtab',
#                   header = TRUE,
#                   row.names = 1)

#import sample data for phenotype information and which genomes to keep
sample_data <- read.csv('R21_Serum_Stool_Merged_Filtered_wGenomeID.csv',
                        header=TRUE,
                        row.names=1)

passed_genomes <- sample_data[sample_data$QC==0,c(p,"genome_id")]
print(passed_genomes)
passed_genomes <- passed_genomes[!is.na(passed_genomes[,p]),]
print(passed_genomes)


#subset geno to only genomes that passed QC and are in the sample set
#panaroo_tab_QCd <- panaroo_tab[,colnames(panaroo_tab) %in% passed_genomes$genome_id]

# REMOVE VARS
#This will drop any of the panaroo groups that are either consistently present or consistently absent across all the strains since these are not informative
#panaroo_tab_QCd <- panaroo_tab_QCd[-index_of_non_var_sites(panaroo_tab_QCd),]

# INDELS
load('INDEL_parsed.RData')
indel_parsed <- parsed
#print(indel_parsed)
rm(parsed)
indel_mat_code <- indel_parsed$bin$mat
#print(indel_mat_code)
indel_mat_annots <- indel_parsed$bin$annots
#View(indel_mat_annots)

colnames(indel_mat_code)<-gsub('_$', '', colnames(indel_mat_code))
#print(colnames)

# SUBSET ON THOSE IN PANAROO MAT
indel_mat_code_subset <- indel_mat_code[,colnames(indel_mat_code) %in% passed_genomes$genome_id]
#print(indel_mat_code_subset)
View(indel_mat_code_subset)
#print(indel_mat_code_subset)
#indel_mat_code_no_duplicated <- indel_mat_code_subset[!duplicated(indel_mat_code_subset),]
  
#View(indel_mat_code_no_duplicated)

# REMOVE VARS FROM INDEL MAT
to_remove <- index_of_non_var_sites(indel_mat_code_subset)
indel_mat_code_subset <- indel_mat_code_subset[-to_remove,]
indel_mat_annots_subset <- indel_mat_annots[-to_remove,]
#View(indel_mat_code_no_duplicated)
#indel_mat_code_no_duplicated <- indel_mat_code_subset[!duplicated(indel_mat_code_subset),]
#View(indel_mat_code_no_duplicated)

# REMOVE LOW IMPACT VARIANTS FROM INDEL MAT & ANNOTS
indel_mat_annots_subset <- indel_mat_annots_subset[!grepl('LOW', row.names(indel_mat_code_subset)),]
indel_mat_code_subset <- indel_mat_code_subset[!grepl('LOW', row.names(indel_mat_code_subset)),]
indel_mat_code_subset <- indel_mat_code_subset[!duplicated(indel_mat_code_subset),]
#View(indel_mat_code_no_duplicated)
View(indel_mat_code_subset)

# SNPS
load('SNP_parsed.RData')
snp_parsed <- parsed
rm(parsed)
snp_mat_code <- snp_parsed$bin$mat
snp_mat_annots <- snp_parsed$bin$annots
#print(snp_mat_code)
#print(snp_mat_annots)

colnames(snp_mat_code)<-gsub('_$', '', colnames(snp_mat_code))
#print(colnames(snp_mat_code))

# SUBSET TO THOSE IN PANAROO MAT
snp_mat_code_subset <- snp_mat_code[,colnames(snp_mat_code) %in% passed_genomes$genome_id]
#print(snp_mat_code_subset)

# REMOVE VARS FROM SNP MAT & ANNOTS
to_remove = index_of_non_var_sites(snp_mat_code_subset)
snp_mat_code_subset = snp_mat_code_subset[-to_remove,]
snp_mat_annots_subset = snp_mat_annots[-to_remove,]

# REMOVE LOW IMPACT VARIANTS FROM SNP MAT & ANNOTS
snp_mat_annots_subset = snp_mat_annots_subset[!grepl('LOW', row.names(snp_mat_code_subset)),]
snp_mat_code_subset = snp_mat_code_subset[!grepl('LOW', row.names(snp_mat_code_subset)),]
View(snp_mat_code_subset)
snp_mat_code_subset <- snp_mat_code_subset[!duplicated(snp_mat_code_subset),]
View(snp_mat_code_subset)

# COMBINE MATRICES:
# FIRST, CHECK TO MAKE SURE COLUMN NAMES ALIGN
check1 <- sum(colnames(indel_mat_code_subset) != colnames(snp_mat_code_subset))
#print(check1)

if (check1 == 0){
  geno = rbind(snp_mat_code_subset, indel_mat_code_subset)
}
print(geno)
View(geno)

# GROUP BY GENE
gene = c(as.character(snp_mat_annots_subset$locus_tag), as.character(indel_mat_annots_subset$locus_tag)) # same order as geno mat!
vars = c(as.character(row.names(snp_mat_code_subset)), as.character(row.names(indel_mat_code_subset)))
print(gene)
View(vars)
print(vars)

gene_key = cbind(SNP = vars, GENE = gene)
print(gene_key)
View(gene_key)


#PHENO ----
#subset to just phenotype information
pheno<-passed_genomes[,1,drop=FALSE]
#print(pheno)
View(pheno)
rownames(pheno)<-passed_genomes$genome_id
View(rownames(pheno))

pheno <- pheno[!(rownames(pheno) %in% setdiff(rownames(pheno), colnames(geno))), , drop=FALSE]
#print(pheno)
View(pheno)

#check to make sure all samples in pheno are in geno and vice versa

# colnames(geno) %in% rownames(pheno) #all the samples are in there
# colnames(geno) == rownames(pheno) #they are not in the same order
reorder_ids <- match(colnames(geno),rownames(pheno))
pheno<-pheno[reorder_ids,,drop=FALSE]
View(reorder_ids)
View(pheno)
# colnames(geno) == rownames(pheno)

#now they're in the same order

# TREE ----
tree = read.tree('2021_02_02_cdiff_630_genome_aln_w_alt_allele_unmapped_R21_Serum_Stool_subset.treefile')

tree$tip.label <-gsub('_$', '', tree$tip.label)

tree = drop.tip(tree, setdiff(tree$tip.label, colnames(geno)))
#print(tree)

tip_order_changes <- function(tr){
  tip_log <- tr$edge[, 2] <= ape::Ntip(tr)
  ordered_tips <- tr$edge[tip_log, 2]
  sum(tr$tip.label != tr$tip.label[ordered_tips]) > 0
}
print(tip_order_changes)

tip_order_changes(tree) # if false, don't need to rearrange tree

# RUN HOGWASH

# hogwash(pheno = as.matrix(pheno[tree$tip.label,, drop = FALSE]),
#         geno = as.matrix(t(geno))[tree$tip.label,],
#         tree = tree,
#         file_name = 'grouped_post-ar_Project_U01_SB_Cytokine',
#         dir = '../results',
#         group_genotype_key = gene_key,
#         grouping_method = "post-ar",
#         test='synchronous')

hogwash(pheno = as.matrix(pheno[tree$tip.label,, drop = FALSE]),
        geno = as.matrix(t(geno))[tree$tip.label,],
        tree = tree,
        file_name = 'Project_U01_SB_Cytokine_hogwash_out',
        dir = './')

system.time({hogwash()})


