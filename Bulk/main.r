# -*- coding: utf-8 -*-
# =============================================================================
# 基因差异表达分析 (limma)
# 流程：读取表达矩阵 → 标准化 → 线性模型拟合 → 差异基因筛选 → 可视化
# =============================================================================

library(limma)
library(ggplot2)
library(here)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
library(pathview)

base_dir    <- file.path(here::here(), "复现", "基因数据分析")
inputed_dir <- file.path(base_dir, "原始数据")
output_dir  <- file.path(base_dir, "输出")

# 读取数据，第一列为基因名
raw_data    <- read.csv(file.path(inputed_dir, "expr_set.csv"))
gene_names  <- raw_data$name
expr_matrix <- raw_data[, sapply(raw_data, is.numeric)]
rownames(expr_matrix) <- gene_names

# 压缩数据范围，让分布更接近正态，当数据 ＞ 50 时启用
if (max(expr_matrix) > 50) {
  expr_matrix <- log2(expr_matrix + 1)
}

# 芯片/样本间分位数标准化（行名自动继承）
expr_normalized <- normalizeBetweenArrays(expr_matrix)

# 前 6 列为对照组 (CK)，后 6 列为处理组 (A)
group  <- factor(c(rep("CK", 6), rep("A", 6)), levels = c("CK", "A"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 拟合线性模型
fit <- lmFit(expr_normalized, design, method = "ls")

# 构建对比：处理组 vs 对照组
contrast_matrix <- makeContrasts(A - CK, levels = design)
fit_contrast    <- contrasts.fit(fit, contrast_matrix)

# 经验贝叶斯方差收缩（借用全局基因信息稳定小样本方差估计）
fit_contrast <- eBayes(fit_contrast)

# 提取全部基因的差异分析结果，FDR 校正；行名已自动对应基因名
deg_results <- topTable(fit_contrast, number = Inf, adjust.method = "fdr")

# 阈值：|log2FC| >= 1 且 P.Value < 0.05
deg_results$change <- ifelse(
  deg_results$P.Value < 0.05 & abs(deg_results$logFC) >= 1,
  ifelse(deg_results$logFC > 0, "Up", "Down"),
  "Stable"
)

message("差异基因统计：")
print(table(deg_results$change))

# 将差异分析数据导出
write.csv(deg_results, file = file.path(output_dir, "results.csv"))

# ========================= 火山图 =========================

fc_threshold   <- 1
pval_threshold <- 0.05
axis_limit     <- 6

# 确保绘图顺序：Stable 在底层，差异基因在上层
deg_results$change <- factor(deg_results$change, levels = c("Stable", "Down", "Up"))
deg_results <- deg_results[order(deg_results$change), ]

volcano_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(P.Value), color = change)) +
  geom_point(alpha = 1, size = 1.2) +
  scale_color_manual(values = c("Down" = "#66CCFF", "Stable" = "grey", "Up" = "#EE0000")) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold),
             linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = -log10(pval_threshold),
             linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_x_continuous(limits = c(-axis_limit, axis_limit),
                     breaks = seq(-axis_limit, axis_limit, 2)) +
  labs(x = "Log2(Fold Change)", y = "-Log10(P.Value)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title     = element_blank(),
    axis.text        = element_text(colour = "black", size = 11),
    axis.title       = element_text(colour = "black", size = 11)
  )

ggsave(file.path(output_dir, "火山图.png"), volcano_plot, width = 5.5, height = 5, dpi = 300)

# ========================= 热图 =========================

# 取 P.Value 最小的前 50 个基因
sig_genes <- head(rownames(deg_results[order(deg_results$P.Value), ]), 50)

# 从标准化矩阵中提取这些基因
heatmap_data <- expr_normalized[sig_genes, ]

# 列注释（样本分组标签）
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(heatmap_data)

ann_colors <- list(
  Group = c(CK = "#66CCFF", A = "#EE0000")
)

pheatmap(
  heatmap_data,
  scale            = "row",
  annotation_col   = annotation_col,
  annotation_colors = ann_colors,
  color            = colorRampPalette(c("#66CCFF", "white", "#EE0000"))(50),
  cluster_cols     = FALSE,
  cluster_rows     = TRUE,
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  fontsize         = 10,
  filename         = file.path(output_dir, "热图.png"),
  width            = 6,
  height           = 8
)

# ========================= 小提琴图 =========================

# 取差异基因（Up + Down）的表达数据
deg_genes <- rownames(deg_results)[deg_results$change != "Stable"]

# 如果差异基因太多，取 P.Value 最小的前 50 个
if (length(deg_genes) > 50) {
  deg_genes <- head(rownames(deg_results[order(deg_results$P.Value), ]), 50)
}

violin_data <- as.data.frame(t(expr_normalized[deg_genes, ]))
violin_data$Group <- group

# 宽表转长表
violin_long <- tidyr::pivot_longer(
  violin_data,
  cols      = -Group,
  names_to  = "Gene",
  values_to = "Expression"
)

violin_plot <- ggplot(violin_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("CK" = "#66CCFF", "A" = "#EE0000")) +
  labs(x = NULL, y = "Log2(Expression)", title = "DEGs Expression Distribution") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, size = 13),
    axis.text        = element_text(colour = "black", size = 11),
    axis.title       = element_text(colour = "black", size = 11)
  )

ggsave(file.path(output_dir, "小提琴图.png"), violin_plot, width = 4, height = 5, dpi = 300)
# ========================= 箱线图 =========================

box_data <- as.data.frame(expr_normalized)
box_long <- tidyr::pivot_longer(
  box_data,
  cols      = everything(),
  names_to  = "Sample",
  values_to = "Expression"
)

# 保持样本原始顺序，按样本名匹配分组
box_long$Sample <- factor(box_long$Sample, levels = colnames(expr_normalized))
group_map <- setNames(as.character(group), colnames(expr_normalized))
box_long$Group <- group_map[as.character(box_long$Sample)]

box_plot <- ggplot(box_long, aes(x = Sample, y = Expression, fill = Group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  scale_fill_manual(values = c("CK" = "#66CCFF", "A" = "#EE0000")) +
  labs(x = NULL, y = "Log2(Expression)", title = "Sample Expression Distribution") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, size = 13),
    axis.text.x      = element_text(colour = "black", size = 9, angle = 45, hjust = 1),
    axis.text.y      = element_text(colour = "black", size = 11),
    axis.title       = element_text(colour = "black", size = 11),
    legend.title     = element_blank()
  )

ggsave(file.path(output_dir, "箱线图.png"), box_plot, width = 7, height = 5, dpi = 300)

# ========================= GO / KEGG 富集分析 =========================

# 提取差异基因名
deg_up   <- rownames(deg_results)[deg_results$change == "Up"]
deg_down <- rownames(deg_results)[deg_results$change == "Down"]
deg_all  <- c(deg_up, deg_down)

# 基因 ID 转换：基因 symbol → Entrez ID
# 如果你的基因名本身就是 Entrez ID 可以跳过这步
gene_map <- bitr(deg_all,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

# -------------------- GO 富集 --------------------

go_result <- enrichGO(
  gene         = gene_map$ENTREZID,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "ALL",        # BP / CC / MF 全部跑
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

# GO 柱状图
go_data <- as.data.frame(go_result)

# 每个 ontology 取前 10 条
go_plot_data <- go_data %>%
  group_by(ONTOLOGY) %>%
  slice_min(p.adjust, n = 10) %>%
  ungroup() %>%
  arrange(ONTOLOGY, Count)

go_plot_data$Description <- factor(go_plot_data$Description, levels = go_plot_data$Description)

go_bar <- ggplot(go_plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#EE0000", high = "#66CCFF") +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Count", y = NULL, fill = "p.adjust") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(colour = "black", size = 9),
    axis.text.x = element_text(colour = "black", size = 10),
    strip.text  = element_text(size = 11)
  )
ggsave(file.path(output_dir, "GO富集柱状图.png"), go_bar, width = 10, height = 12, dpi = 300)

# GO 气泡图 
go_plot_data$GeneRatio_num <- sapply(go_plot_data$GeneRatio, function(x) {
  parts <- as.numeric(strsplit(x, "/")[[1]])
  parts[1] / parts[2]
})

go_dot <- ggplot(go_plot_data, aes(x = GeneRatio_num, y = Description, 
                                    size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#EE0000", high = "#66CCFF") +
  scale_size_continuous(range = c(2, 6)) +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
  labs(x = "GeneRatio", y = NULL, color = "p.adjust", size = "Count") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(colour = "black", size = 9),
    axis.text.x = element_text(colour = "black", size = 10),
    strip.text  = element_text(size = 11)
  )
ggsave(file.path(output_dir, "GO富集气泡图.png"), go_dot, width = 10, height = 12, dpi = 300)

# -------------------- KEGG 富集 --------------------

kegg_result <- enrichKEGG(
  gene         = gene_map$ENTREZID,
  organism     = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# KEGG 柱状图
kegg_bar_data <- as.data.frame(kegg_result)
kegg_bar_data <- kegg_bar_data[order(kegg_bar_data$Count), ]
kegg_bar_data$Description <- factor(kegg_bar_data$Description, levels = kegg_bar_data$Description)

kegg_bar <- ggplot(kegg_bar_data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#EE0000", high = "#66CCFF") +
  labs(x = "Count", y = NULL) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black", size = 10)
  )
ggsave(file.path(output_dir, "KEGG富集柱状图.png"), kegg_bar, width = 8, height = 6, dpi = 300)

# KEGG 气泡图
kegg_dot <- ggplot(kegg_bar_data, aes(x = Count/as.numeric(sub("/.*", "", BgRatio)), 
                                       y = Description, 
                                       size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#EE0000", high = "#66CCFF") +
  labs(x = "GeneRatio", y = NULL) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black", size = 10)
  )
ggsave(file.path(output_dir, "KEGG富集气泡图.png"), kegg_dot, width = 8, height = 6, dpi = 300)

write.csv(as.data.frame(kegg_result), file.path(output_dir, "KEGG富集结果.csv"), row.names = FALSE)

# ========================= PCA 主成分分析 =========================

# 对标准化后的表达矩阵做 PCA（样本在行，基因在列）
pca_result <- prcomp(t(expr_normalized), scale. = TRUE)

# 提取前两个主成分的方差解释比例
pca_summary <- summary(pca_result)
pc1_var <- round(pca_summary$importance[2, 1] * 100, 1)
pc2_var <- round(pca_summary$importance[2, 2] * 100, 1)

# 构建绘图数据
pca_data <- data.frame(
  PC1    = pca_result$x[, 1],
  PC2    = pca_result$x[, 2],
  Group  = group,
  Sample = colnames(expr_normalized)
)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3.5) +
  stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.6) +
  scale_color_manual(values = c("CK" = "#66CCFF", "A" = "#EE0000")) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    title = "PCA"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, size = 13),
    axis.text        = element_text(colour = "black", size = 11),
    axis.title       = element_text(colour = "black", size = 11),
    legend.title     = element_blank()
  )

ggsave(file.path(output_dir, "PCA图.png"), pca_plot, width = 6, height = 5, dpi = 300)

# ========================= GSEA 富集分析 =========================

# 1. 准备排序的基因列表 (使用全部基因，不设 P 值阈值)
gsea_data <- data.frame(
  SYMBOL = rownames(deg_results),
  logFC  = deg_results$logFC
)

# 转换 ID (SYMBOL -> ENTREZID)
gsea_map <- bitr(gsea_data$SYMBOL, 
                 fromType = "SYMBOL", 
                 toType   = "ENTREZID", 
                 OrgDb    = org.Hs.eg.db)

# 合并数据，处理多对一的情况（如果有重复，取 logFC 绝对值最大的那个）
gsea_data <- merge(gsea_data, gsea_map, by = "SYMBOL")
gsea_data <- gsea_data[order(abs(gsea_data$logFC), decreasing = TRUE), ]
gsea_data <- gsea_data[!duplicated(gsea_data$ENTREZID), ]

# 构建 named vector 并按 logFC 降序排列（GSEA 的标准输入格式）
gene_list <- gsea_data$logFC
names(gene_list) <- gsea_data$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# -------------------- GSEA GO 分析 --------------------
message("开始运行 GSEA GO 分析，全部基因参与，可能需要一点时间...")
gsea_go <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  keyType       = "ENTREZID",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  seed          = 42  # 固定随机种子，保证每次跑结果一样
)

if (nrow(as.data.frame(gsea_go)) > 0) {
  # 挑 NES (标准化富集分数) 绝对值最大的一条通路画图
  top_go_id <- gsea_go$ID[which.max(abs(gsea_go$NES))]
  top_go_title <- gsea_go$Description[which.max(abs(gsea_go$NES))]
  
  p_gsea_go <- gseaplot2(gsea_go, geneSetID = top_go_id, title = top_go_title, 
                         color = "#EE0000", pvalue_table = TRUE)
  ggsave(file.path(output_dir, "GSEA_GO_Top1.png"), p_gsea_go, width = 8, height = 6, dpi = 300)
  
  write.csv(as.data.frame(gsea_go), file.path(output_dir, "GSEA_GO_结果.csv"), row.names = FALSE)
  message("GSEA GO 完成！找到了 ", nrow(as.data.frame(gsea_go)), " 条显著变化的整体趋势。")
} else {
  message("GSEA GO 无显著结果。")
}

# -------------------- GSEA KEGG 分析 --------------------
message("开始运行 GSEA KEGG 分析...")
gsea_kegg <- gseKEGG(
  geneList      = gene_list,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  seed          = 42
)

if (nrow(as.data.frame(gsea_kegg)) > 0) {
  top_kegg_id <- gsea_kegg$ID[which.max(abs(gsea_kegg$NES))]
  top_kegg_title <- gsea_kegg$Description[which.max(abs(gsea_kegg$NES))]
  
  p_gsea_kegg <- gseaplot2(gsea_kegg, geneSetID = top_kegg_id, title = top_kegg_title, 
                           color = "#66CCFF", pvalue_table = TRUE)
  ggsave(file.path(output_dir, "GSEA_KEGG_Top1.png"), p_gsea_kegg, width = 8, height = 6, dpi = 300)
  
  write.csv(as.data.frame(gsea_kegg), file.path(output_dir, "GSEA_KEGG_结果.csv"), row.names = FALSE)
  message("GSEA KEGG 完成！找到了 ", nrow(as.data.frame(gsea_kegg)), " 条显著变化的整体趋势。")
} else {
  message("GSEA KEGG 无显著结果。")
}

# ========================= KEGG 通路图 =====================

pathview(
  gene.data  = gene_list,
  pathway.id = "hsa03010",
  species    = "hsa",
  limit      = list(gene = max(abs(gene_list)), cpds = 1),
  kegg.dir   = output_dir,
  out.suffix = "ribosome"
)

# pathview 默认输出到工作目录，手动把生成的图移到输出文件夹
file.rename(
  from = "hsa03010.ribosome.png",
  to   = file.path(output_dir, "KEGG_Ribosome通路图.png")
)

pathview(
  gene.data  = gene_list,
  pathway.id = "hsa04141",
  species    = "hsa",
  limit      = list(gene = max(abs(gene_list)), cpds = 1),
  kegg.dir   = output_dir,
  out.suffix = "ER_processing"
)

file.rename(
  from = "hsa04141.ER_processing.png",
  to   = file.path(output_dir, "KEGG_内质网蛋白加工通路图.png")
)