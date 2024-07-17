library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


#CpG site IDs. It was done here by adding the values manually. However one CAN IMPORT THE FILE WITH IMPORTANT CPG SITES AS WELL AND start FROM THERE INSTEAD.
cpg_sites <- c("cg11995374", "cg06592160", "cg06317989", "cg00640147", "cg22467810", "cg14499563", "cg15439294", "cg16956804", "cg01889461", "cg05057976")


#Annotation data, finding exact location on chromosome
cpg_annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

filtered_annotation <- cpg_annotation[rownames(cpg_annotation) %in% cpg_sites, ]
filtered_annotation #only for cpg sites


#grange creation
chromosome <- filtered_annotation$chr
start_position <- as.numeric(filtered_annotation$pos)
end_position <- as.numeric(filtered_annotation$pos)
strand <- filtered_annotation$strand
genomic_region <- GRanges(seqnames = chromosome, ranges = IRanges(start = start_position, end = end_position), strand = strand)

#Mapping to gene
gene_annotations <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)# Retrieve gene annotations for hg38

# Find overlapping genes
overlapping_genes <- findOverlaps(genomic_region, gene_annotations)

# Print the associated gene(s)
if (length(overlapping_genes) > 0) {
    print(gene_annotations[subjectHits(overlapping_genes)])
} else {
    print("No genes found in the specified genomic region.")
}

gene_ids <- gene_annotations[subjectHits(overlapping_genes)]$gene_id
gene_ids
gene_names <- select(org.Hs.eg.db, keys = gene_ids, columns = "SYMBOL")
print(gene_names)