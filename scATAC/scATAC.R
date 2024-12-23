data <- read.table("ATAC/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz", header = FALSE, sep = "\t")


data_GRanges = GRanges(
  seqnames = Rle(data$V1),
  ranges = IRanges(start = data$V2, end = data$V3),
  Cell_ID = data$V4
)

flag =NULL

for (i in 1:length(gene_location))
{
  if (as.character(strand(gene_location[i])) == "+") {
    tss <- start(ranges(gene_location[i]))
  } else if (as.character(strand(gene_location[i])) == "-") {
    tss <- end(ranges(gene_location[i]))
  }
  gene_location$TSS[i]=tss
  if (tss %in% start(ranges(all_TSS)))
  {
    flag=c(flag,i)
  }
}


overlaps = findOverlaps(all_TSS,data_GRanges,maxgap = 200, minoverlap = 0)
no_promoter= data_GRanges[-overlaps@to]

RNA_data = read.csv('RNA/RNA_Count.csv')
rownames(RNA_data)<-RNA_data$X


cell_IDs = colnames(RNA_data[-1])

all_cell_promoters=GRanges()

find_no_promoter<-function(cell_ID)
{
    gene_names = RNA_data$X[RNA_data[cell_ID]>0]
    gene_TSS = gene_location[gene_location$symbol%in%gene_names]
    gene_TSS_GRanges=GRanges(
      seqnames = seqnames(gene_TSS),
      ranges = IRanges(start = gene_TSS$TSS),
    )
    cell_GRanges = data_GRanges[data_GRanges$Cell_ID== gsub("\\.", "-", cell_ID)]
    cell_overlaps = findOverlaps(gene_TSS_GRanges,cell_GRanges)
    cell_promoter= cell_GRanges[cell_overlaps@to]
    all_cell_promoters=c(all_cell_promoters,cell_promoter)
}

cl <- makeCluster(14)
registerDoParallel(cl)
foreach(i = cell_IDs) %dopar% find_no_promoter(i) 
stopCluster(cl)

for (cell_ID in cell_IDs){
  gene_names = RNA_data$X[RNA_data[cell_ID]>0]
  gene_TSS = gene_location[gene_location$symbol%in%gene_names]
  gene_TSS_GRanges=GRanges(
    seqnames = seqnames(gene_TSS),
    ranges = IRanges(start = gene_TSS$TSS),
  )
  cell_GRanges = data_GRanges[data_GRanges$Cell_ID== gsub("\\.", "-", cell_ID)]
  cell_overlaps = findOverlaps(gene_TSS_GRanges,cell_GRanges)
  cell_promoter= cell_GRanges[cell_overlaps@to]
  all_cell_promoters=c(all_cell_promoters,cell_promoter)
}
saveRDS(no_promoter,file = "no_promoter.rds")
saveRDS(all_cell_promoters,file = 'ATAC/all_promoters.rds')

library(GenomicRanges)
library(Biostrings)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38) 
genome <- BSgenome.Hsapiens.UCSC.hg38

sequences <- DNAStringSet()

# 保存序列名称

# 遍历提取每个范围的序列，并将其保存到FASTA格式

# library(foreach)
# library(doParallel)
# library(GenomicRanges)
# cl <- makeCluster(10)
# registerDoParallel(cl)
# for (i in seq_along(no_promoter)) 
# foreach(i = 1:200) %dopar% 
  
set.seed(123) 
random_integers <- sample(1:length(no_promoter), 30000, replace = FALSE)
for (i in random_integers)
{
  gr <- no_promoter[i]
  seq <- getSeq(genome, gr)
  seq_name <- paste0(as.character(seqnames(gr)), ":", 
                     start(ranges(gr)), "-", end(ranges(gr)), "|", 
                     gr$Cell_ID)
  names(seq)<-seq_name
  #writeXStringSet(seq, "promoter_sequences.fasta",append = TRUE)
  sequences<-append(sequences,seq)
  print(length(sequences))
}


writeXStringSet(sequences, "no_promoter_sequences.fasta")
