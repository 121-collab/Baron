# Baron et al., 2016 人胰腺 scRNA-seq 复现仓库
> 个人生信复现笔记 —— 复现 **Baron et al., 2016** 胰腺单细胞转录组 (scRNA-seq) 主要图谱  
> 目标：在公司服务器 + 本地 R 环境下，尽量接近原文 Fig1–4（人 / 小鼠）的分析流程和图像效果。
## 📄 原文与数据来源
- 论文：Baron et al., 2016 — *A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas*  
- 数据库：GEO  
  - 人胰腺：**GSE84133**（Human pancreas）  
  - 小鼠胰腺：**GSE84133**（Mouse pancreas，同一系列中的 mouse 样本）  
- 使用数据类型：
  - 作者提供的 **UMI 计数矩阵（UMI count matrix, csv.gz）**
  - 不再从 FASTQ 重新对齐，只在 UMI 矩阵基础上完成 Seurat 下游分析  
> 说明：本仓库定位为“下游分析复现”，不涵盖 Cell Ranger / 比对 / 定量环节。

---

## 🧰 运行环境（复现时使用）

- 服务器：Ubuntu 24.04 + 公司共享服务器  
- R 版本：R 4.3.x  
- 主要 R 包：
  - `Seurat (v5)`、`Matrix`、`data.table`、`ggplot2`、`patchwork`
  - `pheatmap`、`dplyr`、`tidyr` 等

---

## 📁 仓库结构（简要）

```text
Baron2016_pancreas_scRNA_repro/
├── README.md                       # 本说明文档
├── Baron2016_pancreas_pipeline.md  # 全流程笔记（Markdown 版本）
│
├── scripts/                        # 分步骤 R 脚本（已在流水线中整合）
│   ├── 01_build_seurat.R           # 读取 Baron UMI 矩阵，构建 Seurat 对象（人 + 小鼠）
│   ├── 02_fig1_markers.R           # Fig1：UMAP/tSNE + DotPlot/Top markers 热图
│   ├── 03_fig2_corr.R              # Fig2：Human/Mouse donor × celltype 相关性热图
│   └── 04_fig3_ductal_fig4_beta_PCA.R  # Fig3：ductal PCA/PC1；Fig4：beta 细胞 PCA/PC1
│
├── Baron_results/                  # 人胰腺相关结果（Seurat 对象 + 图）
│   ├── human_seurat_celltyped.rds          # 人胰腺 Seurat 对象（已注释 celltype）
│   ├── Human_DotPlot_markers_paper.pdf     # Fig1A/B 风格 DotPlot
│   ├── Human_TopMarkers_Heatmap_paper.pdf  # Fig1 热图（Top markers）
│   ├── Human_Fig2*_corr_heatmap.pdf        # Fig2B–C 风格相关性热图
│   ├── Human_Fig3*_Ductal_PCA_PC1_PC2.pdf  # Fig3 ductal PCA/PC1
│   └── Human_Fig4*_Beta_PCA_PC1_PC2*.pdf   # Fig4 beta PCA/PC1（含 moving average 等）
│
├── Baron_mouse_figs/               # 小鼠胰腺相关图像
│   ├── Mouse_DotPlot_markers_paper.pdf
│   ├── Mouse_TopMarkers_Heatmap_paper.pdf
│   ├── Mouse_Fig2*_corr_heatmap.pdf
│   ├── Mouse_Fig3*_Ductal_PCA_PC1_PC2.pdf
│   └── Mouse_Fig4*_Beta_PCA_PC1_PC2*.pdf
│
└── .gitignore                      # 忽略大体积中间文件（例如 *.rds, 临时日志等）
