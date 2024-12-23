test_val = read.csv('test_val_results.csv')
test_conf_matrix <- table(Predicted = test_val$test_predictions, Actual = test_val$test_labels)
val_conf_matrix<- table(Predicted = test_val$val_predictions, Actual = test_val$val_labels) 

library(ggplot2)
library(dplyr)

test_metrics <- get_metrics(test_conf_matrix) %>% mutate(Set = "Test")
val_metrics <- get_metrics(val_conf_matrix) %>% mutate(Set = "Validation")


get_metrics <- function(conf_matrix) {
  TN <- conf_matrix["0", "0"]
  FP <- conf_matrix["0", "1"]
  FN <- conf_matrix["1", "0"]
  TP <- conf_matrix["1", "1"]
  data.frame(Metric = c("TP", "TN", "FP", "FN"), Count = c(TP, TN, FP, FN))
}

# 获取测试和验证集的统计数据
test_metrics <- get_metrics(test_conf_matrix) %>% mutate(Set = "Test")
val_metrics <- get_metrics(val_conf_matrix) %>% mutate(Set = "Validation")
# 合并数据
metrics_data <- rbind(test_metrics, val_metrics)

# 绘制柱状图

ggplot(metrics_data, aes(x = Metric, y = Count, fill = Set)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(0.9), size = 4, family = "Times New Roman") +
  scale_fill_manual(values = c("Test" = "orange", "Validation" = "#66C2A5")) + 
  theme_minimal(base_family = "Times New Roman") + 
  theme(
    axis.text = element_text(color = "black", family = "Times New Roman", size = 12),  # 坐标刻度颜色和字体
    axis.title = element_blank(),  # 去掉X/Y轴标题
    axis.line = element_line(color = "black"),  # 添加x和y轴线，黑色
    legend.title = element_blank(),  # 去掉图例标题
    legend.text = element_text(family = "Times New Roman", size = 12),  # 图例字体设置
    panel.grid.major = element_blank(), # 去掉大网格线
    panel.grid.minor = element_blank()  # 去掉小网格线
  )

# 

test_index_TP = test_val[test_val$test_predictions==1 & test_val$test_labels==1,'test_index']
val_index_TP = test_val[test_val$val_predictions==1 & test_val$val_labels==1,'val_index']

cell_position = read.csv('Cell_position.csv')

test_TP_cell = cell_position[test_index_TP+1,'Cell_ID']
val_TP_cell =cell_position[val_index_TP+1,'Cell_ID']



#_______________绘制验证集和独立集细胞ID的的VN________________

a <-list('test_TP_cell'= test_TP_cell,
         'val_TP_cell'= val_TP_cell )

color <- brewer.pal(3, "Set3")

library(VennDiagram)
library(RColorBrewer)

# 数据
a <- list('test_TP_cell' = test_TP_cell,
          'val_TP_cell' = val_TP_cell)

# 配色
color <- brewer.pal(3, "Set3")

# 修正后的 Venn Diagram
venn.diagram(
  x = a,
  category.names = c("test_TP_cell", "val_TP_cell"),
  filename = 'venn.tiff',
  output = TRUE,
  imagetype = "tiff",
  resolution = 600,
  compression = "lzw",
  lwd = 2,
  lty = 2,
  fill = color[1:2],  # 限制为两个颜色，因为有两个集合
  col = c("red", "blue"),
  cex = 1.5,
  fontfamily = "serif",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20),  # 仅设置两个位置值
  cat.dist = c(0.05, 0.05),  # 同样仅为两个值
  rotation.degree = 0 # 修正 rotation 参数
)

test_val_intersect = intersect(val_TP_cell,test_TP_cell)
test_TP_cell =  setdiff(unique(test_TP_cell),test_val_intersect)
val_TP_cell =  setdiff(unique(val_TP_cell),test_val_intersect)
library(Seurat)
pbmc_cluster$custom_label <- "Unlabeled" 
pbmc_cluster$custom_label[val_TP_cell] <- "TP"
pbmc_cluster$custom_label[test_TP_cell] <- "TP"
pbmc_cluster$custom_label[test_val_intersect] <- "TP"

Idents(pbmc_cluster) <- "custom_label"
DimPlot(pbmc_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 


features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
             "CD8A")
features_location =gene_location[gene_location$symbol%in%features]
features_TSS <- GRanges(
  seqnames = seqnames(features_location),
  ranges = IRanges(
    start = ifelse(strand(features_location) == "+",
                   start(features_location),
                   end(features_location)),
    end = ifelse(strand(features_location) == "+",
                 start(features_location),
                 end(features_location))
  ),
  strand = strand(features_location),
  gene_id = features_location$gene_id,
  symbol = features_location$symbol
)

TP_infromation = cell_position[c(test_index_TP+1,val_index_TP+1),]

position_split <- strsplit(TP_infromation$position, "[:-]")
chromosomes <- sapply(position_split, `[`, 1)
start_positions <- as.numeric(sapply(position_split, `[`, 2))
end_positions <- as.numeric(sapply(position_split, `[`, 3))

TP_GR <- GRanges(
  seqnames = chromosomes,
  ranges = IRanges(start = start_positions, end = end_positions)
  
)
TP_GR$Cell_ID =TP_infromation$Cell_ID
TP_GR$index = c(test_index_TP,val_index_TP)




feature_overlaps = findOverlaps(features_TSS,TP_GR)
gene_cell_promoter <- data.frame(gene = character(), cell = character(), index=integer(),stringsAsFactors = FALSE)
for (i in c(1:length(feature_overlaps)))
{
  gene = features_TSS$symbol[feature_overlaps@from[i]]
  cell = TP_GR$Cell_ID[feature_overlaps@to[i]]
  index = TP_GR$index[feature_overlaps@to[i]]
  gene_cell_promoter<-rbind(gene_cell_promoter,data.frame(gene=gene,cell=cell,index=index))
}
gene_cell_promoter$cluster =  pbmc_cluster$seurat_clusters[gene_cell_promoter$cell]
gene_cell_promoter$mark_cluster = c(6,2,2)







