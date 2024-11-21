setwd("D:/tss_atac")
library(Seurat)
library(SeuratObject)

#--------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("Biostrings")
install.packages("GenomicRanges")

library(Biostrings)
#下载hg38.fa文件
hg38 <- "C:/Users/1/Desktop/tss_atac/hg38.fa"
hg38_sequences <- readDNAStringSet(hg38)

#--------------------------------------------------------------------------------

#数据输入
rna <- load("C:/Users/1/Desktop/pbmc.rna.rda")

atac_data <-"C:/Users/1/Desktop/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
atac <- read.table(atac_data,sep="\t")

#提取id  行、列名
RNA_Gene = row.names(pbmc.rna)
RNA_Cell_ID = colnames(pbmc.rna)

ATAC_cell_ID  = unique(atac$V4)
# 提取 RNA 数据和 ATAC 数据中共同的细胞 ID
common_cells <- intersect(RNA_Cell_ID, ATAC_cell_ID)
# 从 ATAC 数据中选择与 RNA 数据中的细胞 ID 匹配的行
use_ATAC = atac[atac$V4 %in% RNA_Cell_ID,]
#复制ATAC数据
co_use_ATAC <- use_ATAC
#给ATAC 数据添加一列数据，这列数据为行号
co_use_ATAC$V6 <- seq_len(nrow(co_use_ATAC))
# 数据框名为 co_use_ATAC，要保存为 BED 文件
write.table(co_use_ATAC, "co_use_ATAC.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 将数据框保存为 CSV 文件
write.csv(co_use_ATAC, file = "co_use_ATAC.csv", row.names = FALSE)


# 找到高表达基因
pbmc <- FindVariableFeatures(pbmc.rna, selection.method = "vst",nfeatures = 2000 )

top_2000_genes <- head((pbmc@assays$RNA@var.features), 2000)
top_2000_genes

# 提取这些基因的表达数据
expr_data_top_2000 <- pbmc@assays$RNA@counts[top_2000_genes, ]
expr_data_top_2000

#-----------------------------------------------------------------------------------------------
#人类tss数据输入
tss_data <- read.table(gzfile("C:/Users/1/Desktop/refFlat.txt.gz"))

# 找出前2000基因tss
top_2000_genes_tss <- tss_data[tss_data$V1 %in% top_2000_genes, ]
top_2000_genes_tss

num_rows <- nrow(top_2000_genes_tss)
num_rows

#-----------------------------------------------------------------------------------------------

gene_infom1 <- list()

# 遍历每一行数据
for (i in 1:nrow(top_2000_genes_tss)) {
  line1 <- top_2000_genes_tss[i, 1]  

  if (!(line1 %in% gene_infom1[["V1"]])) {
    gene_infom1[[line1]] <- top_2000_genes_tss[i, ]
  }
}  

genes_2000_tss_df1 <- do.call(rbind, gene_infom1)
nrow(genes_2000_tss_df1)

# 保存
write.csv(genes_2000_tss_df1, "gene_first_lines1.csv", row.names = FALSE)


----------------------------------------------------------------------------------------------
# 将所有负链改为正链

for (i in 1:nrow(genes_2000_tss_df1)) {
  if (genes_2000_tss_df1[i, "V4"] == "-") {
    genes_2000_tss_df1[i, "V4"] <- "+"
    
    temp_start <- genes_2000_tss_df1[i, "V6"]
    genes_2000_tss_df1[i, "V6"] <- genes_2000_tss_df1[i, "V5"]
    genes_2000_tss_df1[i, "V5"] <- temp_start
  }
}

head(genes_2000_tss_df1)
# 保存
write.csv(genes_2000_tss_df1, "genes_2000_tss_df1.csv",  row.names = FALSE)
#----------------------------------------------------------------------------------------

#将数据保存为bed格式文件
#复制数据
co_genes_2000_tss_df <- genes_2000_tss_df1
#完成第一和第三列换列
temp_col <- co_genes_2000_tss_df[, 1]
co_genes_2000_tss_df[, 1] <- co_genes_2000_tss_df[, 3]
co_genes_2000_tss_df[, 3] <- temp_col

#完成第二和第六列换列
temp_col <- co_genes_2000_tss_df[, 2]
co_genes_2000_tss_df[, 2] <- co_genes_2000_tss_df[, 6]
co_genes_2000_tss_df[, 6] <- temp_col

#完成第三和第五列换列
temp_col <- co_genes_2000_tss_df[, 3]
co_genes_2000_tss_df[, 3] <- co_genes_2000_tss_df[, 5]
co_genes_2000_tss_df[, 5] <- temp_col
co_genes_2000_tss_df
nrow(co_genes_2000_tss_df)
#  将数据框导出为BED格式的文件
write.table(co_genes_2000_tss_df, file = "co_genes_1613_11_col_tss.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

genes_name_1613 = row.names(co_genes_2000_tss_df)
co_genes_1613_tss_df <- co_genes_2000_tss_df
#添加行号，这个行号被存储在V12列中
co_genes_1613_tss_df$V12 <- seq_len(nrow(co_genes_1613_tss_df))
write.table(co_genes_1613_tss_df, file = "co_genes_1613_12_col_tss.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#------------------------------------------------------------------------------------
library(dplyr)
top_2000_genes_tss_atac_data_list <- list()

# 遍历 use_ATAC 中的每一行
for (i in 1:nrow(use_ATAC)) {
  use_ATAC_fragment <- use_ATAC[i, ] # 获取当前行片段信息
  tss_atac_matches <- data.frame()  
  # 遍历 genes_2000_tss_df 中的每一行
  for (j in 1:nrow(copy_genes_2000_tss_df)) {
    tss_gene <- copy_genes_2000_tss_df[j, ] # 获取当前行基因的信息

      if (use_ATAC_fragment$V2 <= tss_gene$V6 && use_ATAC_fragment$V3 >= tss_gene$V5) 
      {
        tss_atac_matches <- bind_rows(tss_atac_matches, use_ATAC_fragment)  # 将匹配的片段合并到 matches 数据框中
        break
      }
    } 
  #添加
  top_2000_genes_tss_atac_data_list <- c(top_2000_genes_tss_atac_data_list, list(tss_atac_matches))
}


top_2000_genes_tss_atac_data <- bind_rows(top_2000_genes_tss_atac_data_list)

promoter_line_number <- read.table("C:/Users/1/Documents/promoter_line_number.txt",header = FALSE, sep = "\t")
###----------------------------------------------atac的启动子-----------------------------------------------
promoter_row_atac <- co_use_ATAC[co_use_ATAC$V6 %in% promoter_line_number$V1, ]
write.table(promoter_row_atac, 
            file = "promoter_row_atac.bed",
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)


###----------------------------------------------atac的非启动子-----------------------------------------------
unpromoter_row_atac <- co_use_ATAC[!co_use_ATAC$V6 %in% promoter_line_number$V1, ]
write.table(unpromoter_row_atac, 
            file = "unpromoter_row_atac.bed",
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

getwd()


###----------------------------------------------atac的启动子序列-----------------------------------------------
#(base) [jyli@FAT-SERVER bin]$ bedtools getfasta -fi "/data/jyli/Seurat and Archr data/hg38.fa" -bed "/data/jyli/Seurat and Archr data/promoter_unpromoter_tss/promoter_row_atac.bed" -fo "/data/jyli/Seurat and Archr data/promoter_unpromoter_tss/promoter_sequence.fasta"

###----------------------------------------------atac的非启动子序列-----------------------------------------------
#(base) [jyli@FAT-SERVER bin]$ bedtools getfasta -fi "/data/jyli/Seurat and Archr data/hg38.fa" -bed "/data/jyli/Seurat and Archr data/promoter_unpromoter_tss/unpromoter_row_atac.bed" -fo "/data/jyli/Seurat and Archr data/promoter_unpromoter_tss/unpromoter_sequence.fasta"


###通过细胞barcode，关联pbmc.rna和co_use_ATAC
library(dplyr)

#提取pbmc.rna中的元数据
pbmc.rna_metadata_df <- pbmc.rna@meta.data
head(pbmc.rna_metadata_df,5)

rna_counts_matrix <- pbmc.rna@assays[["RNA"]]@counts
rna_counts_matrix

head(rna_counts_matrix,5)

rna_t_counts_matrix <- t(rna_counts_matrix)
rna_t_counts_matrix_df <- as.data.frame(rna_t_counts_matrix)
rna_t_counts_matrix_df
head(rna_t_counts_matrix_df,1)
#添加cell_id列
rna_t_counts_matrix_cr_df<-rna_t_counts_matrix_df
rna_t_counts_matrix_cr_df$Cell_ID <- rownames(rna_t_counts_matrix_df)
rna_t_counts_matrix_cr_df$Cell_ID
#读取启动子的atac数据
promoter_co_use_ATAC <- read.csv("C:/Users/1/Documents/process_data/promoter_co_use_ATAC.csv", header = TRUE, quote = "")
promoter_co_use_ATAC$X.V4.
# 移除启动子的atac数据中细胞id中的双引号
promoter_co_use_ATAC$X.V4. <- gsub("\"", "", promoter_co_use_ATAC$X.V4.)

combined_rna_promoter_ATAC_data <- inner_join(rna_t_counts_matrix_cr_df, promoter_co_use_ATAC, by = c("Cell_ID" = "X.V4."))

###---------------------------------------------------fa格式文件中序列的最长和最短长度----------------------------------------------------------
library(Biostrings)

# 读取文件
promoter_sequence_fasta_file <- "C:/Users/1/Desktop/tss_atac/promoter_sequence.fasta"
unpromoter_sequence_fasta_file <-"C:/Users/1/Desktop/tss_atac/unpromoter_sequence.fasta"

p_sequences <- readDNAStringSet(promoter_sequence_fasta_file)
unp_sequences <- readDNAStringSet(unpromoter_sequence_fasta_file)

# 计算序列长度
p_sequence_lengths <- width(p_sequences)
unp_sequence_lengths <- width(unp_sequences)
print(p_sequence_lengths)
print(unp_sequence_lengths)

df_p_sequence_lengths <- data.frame(col1 = p_sequence_lengths)
df_unp_sequence_lengths <- data.frame(col1 = unp_sequence_lengths)

# 计算向量中小于500的元素个数
p_count <- sum(unp_sequence_lengths > 100)
print(p_count)

# 打印最长和最短序列长度
cat("Maximum p_sequence length:", max(p_sequence_lengths), "\n")
#Maximum p_sequence length: 4830
cat("Minimum p_sequence length:", min(p_sequence_lengths), "\n")
#Minimum p_sequence length: 15

cat("Maximum unp_sequence length:", max(unp_sequence_lengths), "\n")
#Maximum unp_sequence length: 5041
cat("Minimum unp_sequence length:", min(unp_sequence_lengths), "\n")
#Minimum unp_sequence length: 10

#---------------------------------------------------------------------------------------------------------------------
#把文件中长度小于100且大于300的序列删除
library(Biostrings)
promoter_sequence_delete <- readDNAStringSet("D:/tss_atac/promoter_sequence.fasta")
unpromoter_sequence_delete <- readDNAStringSet("D:/tss_atac/unpromoter_sequence.fasta")
# 删除长度小于100且大于300的序列
promoter_sequence_delete <- promoter_sequence_delete[width(promoter_sequence_delete) >= 100 & width(promoter_sequence_delete) <= 300]
unpromoter_sequence_delete <- unpromoter_sequence_delete[width(unpromoter_sequence_delete) >= 100 & width(unpromoter_sequence_delete) <= 300]

writeXStringSet(promoter_sequence_delete, "promoter_sequence_delete.fa")
writeXStringSet(unpromoter_sequence_delete, "unpromoter_sequence_delete.fa")
#---------------------------------------------------------------------------------------------------------------------
#启动子序列的分片操作
promoter_sequence_delete <- readDNAStringSet("D:/tss_atac/promoter_sequence_delete.fa")
# 固定长度为100的分片操作
promoter_sequence_fragment_delete <- DNAStringSet()
for (i in 1:length(promoter_sequence_delete)) {
  seq <- promoter_sequence_delete[i]
  seq_length <- nchar(seq)
  num_fragments <- ceiling(seq_length / 100)
  for (j in 1:num_fragments) {
    start <- (j - 1) * 100 + 1
    end <- min(j * 100, seq_length)
    fragment <- subseq(seq, start, end)
    promoter_sequence_fragment_delete <- c(promoter_sequence_fragment_delete, fragment)
  }
}

writeXStringSet(promoter_sequence_fragment_delete, "promoter_sequence_fragment_delete.fa")
#----------------------------------------------------------------------------------------------------------------------
#非启动子序列的分片操作
unpromoter_sequence_unrepeat_delete <- readDNAStringSet("D:/tss_atac/unpromoter_sequence_unrepeat_delete.fa")
# 固定长度为100
unpromoter_sequence_fragment_unrepeat_delete <- DNAStringSet()
for (i in 1:length(unpromoter_sequence_unrepeat_delete)) {
  seq <- unpromoter_sequence_unrepeat_delete[i]
  seq_length <- nchar(seq)
  num_fragments <- ceiling(seq_length / 100)
  for (j in 1:num_fragments) {
    start <- (j - 1) * 100 + 1
    end <- min(j * 100, seq_length)
    fragment <- subseq(seq, start, end)
    unpromoter_sequence_fragment_unrepeat_delete <- c(unpromoter_sequence_fragment_unrepeat_delete, fragment)
  }
}

# 将结果写入文件
writeXStringSet(unpromoter_sequence_fragment_unrepeat_delete, "unpromoter_sequence_fragment_unrepeat_delete.fa")
#-----------------------------------------------------------------------------------------------------------------------
#将分片之后的启动子的序列转化为one-hot编码
library(Biostrings)
# 读取文件
promoter_sequence_fragment_delete <- "D:/tss_atac/promoter_sequence_fragment_delete.fa" 
sequences <- readDNAStringSet(promoter_sequence_fragment_delete)

# one_hot_encode函数
one_hot_encode <- function(sequence) {

  seq_char <- as.character(sequence)
  len <- nchar(seq_char)
  alphabet <- c("A", "C", "G", "T")
  one_hot <- matrix(0, nrow = len, ncol = length(alphabet), dimnames = list(NULL, alphabet))

  for (i in 1:len) {
    nucleotide <- substr(seq_char, i, i)
    if (nucleotide %in% alphabet) {
      col_index <- match(nucleotide, alphabet)
      one_hot[i, col_index] <- 1
    }
  }
  
  return(one_hot)
}
# 调用函数
promoter_one_hot_sequences <- lapply(sequences, one_hot_encode)  
print(promoter_one_hot_sequences[[1]])

promoter_100_one_hot <- lapply(promoter_one_hot_sequences, function(one_hot_matrix) {
  current_length <- nrow(one_hot_matrix)
  if (current_length < 100) {
    rows_to_add <- 100 - current_length
    padding_matrix <- matrix(0, ncol = 4, nrow = rows_to_add)
    one_hot_matrix <- rbind(one_hot_matrix, padding_matrix)
  }
  return(one_hot_matrix)
})
# 打印补充后的矩阵
print(promoter_100_one_hot[[3]])

promoter_100_one_hot_padded <- lapply(promoter_100_one_hot, t)
print(promoter_100_one_hot_padded[[3]])

save_one_hot_to_csv <- function(one_hot_sequences, file_path) {
  lines <- lapply(one_hot_sequences, function(matrix) {
    apply(matrix, 1, function(row) paste(row, collapse = ","))
  })
  lines <- unlist(lapply(lines, function(seq_lines) paste(seq_lines, collapse = "\n")))
  write(lines, file = file_path)
}
# 将one-hot编码保存到CSV文件
save_one_hot_to_csv(promoter_100_one_hot_padded, "promoter_100_one_hot_padded.csv")
#-----------------------------------------------------------------------------------------------------------------------
#将分片之后的非启动子的序列转化为one-hot编码
library(Biostrings)

# 读取文件
unpromoter_sequence_fragment_delete <- "D:/tss_atac/unpromoter_sequence_fragment_unrepeat_delete.fa"
sequences <- readDNAStringSet(unpromoter_sequence_fragment_delete)

# one_hot_encode函数
one_hot_encode <- function(sequence) {
  seq_char <- as.character(sequence)
  len <- nchar(seq_char)
  alphabet <- c("A", "C", "G", "T")
  one_hot <- matrix(0, nrow = len, ncol = length(alphabet), dimnames = list(NULL, alphabet))
  for (i in 1:len) {
    nucleotide <- substr(seq_char, i, i)
    if (nucleotide %in% alphabet) {
      col_index <- match(nucleotide, alphabet)
      one_hot[i, col_index] <- 1
    }
  }
  
  return(one_hot)
}

# 调用函数
unpromoter_one_hot_sequences <- lapply(sequences, one_hot_encode)  
print(unpromoter_one_hot_sequences[[3]])

unpromoter_100_one_hot <- lapply(unpromoter_one_hot_sequences, function(one_hot_matrix) {
  current_length <- nrow(one_hot_matrix)
  if (current_length < 100) {
    rows_to_add <- 100 - current_length
    padding_matrix <- matrix(0, ncol = 4, nrow = rows_to_add)
    one_hot_matrix <- rbind(one_hot_matrix, padding_matrix)
  }
  return(one_hot_matrix)
})
# 打印补充后的矩阵
print(unpromoter_100_one_hot[[3]])
unpromoter_100_one_hot_padded <- lapply(unpromoter_100_one_hot, t)

print(unpromoter_100_one_hot_padded[[3]])

save_one_hot_to_csv <- function(one_hot_sequences, file_path) {
  lines <- lapply(one_hot_sequences, function(matrix) {
    apply(matrix, 1, function(row) paste(row, collapse = ","))
  })
  lines <- unlist(lapply(lines, function(seq_lines) paste(seq_lines, collapse = "\n")))
  write(lines, file = file_path)
}

# 将one-hot编码保存到CSV文件
save_one_hot_to_csv(unpromoter_100_one_hot_padded, "unpromoter_100_one_hot_padded.csv")

