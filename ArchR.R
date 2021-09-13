#ArchR 
#ssh fat02
#cd /share2/pub/zhenggw/zhenggw/HemeFragments/
conda activate R410
R
library(ArchR)
library(parallel)
library(Cairo)
# mclapply

#make arrow files
addArchRGenome("hg19") # hg38, mm9, mm10
addArchRThreads(threads = 16)

pathFragments <- "/share2/pub/zhenggw/zhenggw/HemeFragments/"

inputFiles <- list.files(pathFragments, pattern = ".gz",
        full.names = TRUE)
    names(inputFiles) <- gsub(".fragments.tsv.gz", "", list.files(pathFragments,
        pattern = ".gz"))
    inputFiles <- inputFiles[!grepl(".tbi", inputFiles)]
    inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later 信噪比，根据TSS富集分数进行计算 
  filterFrags = 1000, #ArchR默认会过滤TSS富集得分低于4或唯一比对数小于1000（也就是保留TSS富集得分大于4且唯一比对数大于1000的细胞）#细胞质控
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# > ArrowFiles
# [1] "scATAC_BMMC_R1.arrow"      "scATAC_CD34_BMMC_R1.arrow"
# [3] "scATAC_PBMC_R1.arrow"

#添加doublets信息，后续去除
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

#创建ArchRProject 允许多个Arrow文件整理到单个项目之中
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)



#一些操作
# paste0("Memory Size = ", round(object.size(projHeme1) / 10^6, 3), " MB")
# # "Memory Size = 37.477 MB"

# 我们还可以检查当前的ArchRProject中存放了哪些矩阵数据，这些数据一旦增加后就可以在下游分析中使用。
# getAvailableMatrices(projHeme1)
# # "GeneScoreMatrix" "TileMatrix"
# head(projHeme1$cellNames)
# head(projHeme1$Sample)
# quantile(projHeme1$TSSEnrichment)
#projHeme1$[TAB]

# 从ArchRProject中提取部分细胞
# 以前学习的R语言向量/矩阵/数据框提取数据的方法可以无缝的迁移到ArchRProject上，但区别在于ArchR并不会直接提取数据，而是构建一个新的ArchRProject对象。

# 最简单的方法是根据位置和细胞名
# # 根据位置
# projHeme1[1:100, ]
# # 根据细胞名
# projHeme1[projHeme1$cellNames[1:100], ]

# 复杂一些就是先根据逻辑判断提取细胞名，接着根据细胞名提取数据。例如挑选scATAC_BMMC_R1的细胞，或是TSSEnrichment大于8的细胞。

# # sample name
# idxSample <- BiocGenerics::which(projHeme1$Sample %in% "scATAC_BMMC_R1")
# cellsSample <- projHeme1$cellNames[idxSample]
# projHeme1[cellsSample, ]
# # TSS enrichment score
# idxPass <- which(projHeme1$TSSEnrichment >= 8)
# cellsPass <- projHeme1$cellNames[idxPass]
# projHeme1[cellsPass, ]

# 用getCellColData提取我们需要的两列，其中nFrages需要进行log10运算
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
    x = df[,1],
    y = df[,2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

p

# png("TSS-vs-Frags.png") 
# plot(p)
# dev.off()
# getwd()

library(ggplot2)
ggsave("TSS-vs-Frags.pdf")


# ArchR提供了小提琴图(violin plot)和山脊图(ridge plot)用来展示不同组之间的信息。这两种类型的图可以用一个函数plotGroups进行绘制。除了根据样本进行分组绘图外，还可以使用下游分析得到的分组信息（例如聚类）。
# 根据TSS富集得分为每个样本绘制山脊图。设置plotAs = "ridges"
p1 <- plotGroups(
    ArchRProj = projHeme1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p1

#绘制小提琴图 plotAs = "violin"
p2 <- plotGroups(
    ArchRProj = projHeme1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p2

p3 <- plotGroups(
    ArchRProj = projHeme1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p3

p4 <- plotGroups(
    ArchRProj = projHeme1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p4

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 4, height = 4)

# 绘制样本的TSS富集谱和Fragment大小分布
#Fragments大小分布
p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p1

#TSS富集谱
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
p2

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)

# #这两个图有问题，报错如下，不知道是不是hdf5的原因
# Error in .safelapply(seq_along(uniqGroups), function(x) { :
# Error Found Iteration 1 :
#         [1] "Error in H5Fopen(file) : HDF5. File accessibility. Unable to open file.\n"
#         <simpleError in H5Fopen(file): HDF5. File accessibility. Unable to open file.>
# Error Found Iteration 2 :
#         [1] "Error in H5Fopen(file) : HDF5. File accessibility. Unable to open file.\n"
#         <simpleError in H5Fopen(file): HDF5. File accessibility. Unable to open file.>


# --------------------------------------------------------
# 目前的问题
# Cairo
# parallel的mclapply容易出问题                                       !!!!!解决办法：library(parallel)
# HDF5. File accessibility. Unable to open file.                    !!!!!解决办法：chmod 777 *.arrow




#-------------------------------------------------------------------------------------------------------------

# 保存和加载ArchRProject
saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = FALSE)
# #load等于False不会改变当前环境的projHeme 如果想要将当前的项目复制到新的目录，可以设置load=TRUE

projHeme1 <- readRDS("./Save-ProjHeme1/Save-ArchR-Project.rds")

#------------------------------------------------------------------------------------------------
#从ArchRProject过滤doublets
projHeme2 <- filterDoublets(projHeme1)  #projHemeTmp <- filterDoublets(projHeme1, filterRatio = 1.5) 提高filterRatio会过滤更多的细胞


#------------------------------------------------------------------------------------------------
#ArchR降维分析 隐语义(Latent Semantic Indexing)迭代 （增加LSI的迭代次数，也可以用来处理批次效应）
projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

#使用Harmony矫正批次效应
projHeme2 <- addHarmony(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

#使用ArchR聚类
projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

head(projHeme2$Clusters)
table(projHeme2$Clusters)
#   C1  C10  C11  C12   C2   C3   C4   C5   C6   C7   C8   C9
# 1575  720 1221 1384 1079  307  388  432  353 1276  816  699

library(pheatmap)

saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme1", load = FALSE)
projHeme2 <- readRDS("./Save-ProjHeme1/Save-ArchR-Project.rds")

cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
cM

cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)


#使用scran聚类
projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI",
    method = "scran",
    name = "ScranClusters",
    k = 15
)
# table(projHeme2$Clusters)

#   C1  C10  C11  C12   C2   C3   C4   C5   C6   C7   C8   C9
# 1575  720 1221 1384 1079  307  388  432  353 1276  816  699

#UMAP和t-SNE嵌入
projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
#ArchR实现了一组适用于大多数情况的默认输入参数，要根据不同的细胞数、复杂度和质量进行调整。

# projHeme2@embeddings
# List of length 1
# names(1): UMAP


p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
#区分 样本 UMAP

p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
#区分 簇 UMAP

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)


# 还可以使用plotEmbedding()可视化之前用scran聚类的结果
# p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
# p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "ScranClusters", embedding = "UMAP")
# ggAlignPlots(p1, p2, type = "h")

 # t-Stocastic Neighbor Embedding (t-SNE)
projHeme2 <- addTSNE(
    ArchRProj = projHeme2, 
    reducedDims = "IterativeLSI", 
    name = "TSNE", 
    perplexity = 30
)

# 可以继续使用之前的colorBy和name参数
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

#同样tsne也可以使用scran聚类的cluster 以将scran的结果和Seurat::FindClusters()的结果进行比较
# p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
# p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "ScranClusters", embedding = "TSNE")
# ggAlignPlots(p1, p2, type = "h")
# plotPDF(p1,p2, name = "Plot-tSNE-Sample-ScranClusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)


#Harmony后降维
# 通过UMAP或t-SNE对结果进行嵌入可视化，对迭代LSI结果和Harmony校正结果进行比较，评估Harmony的作用。

# 保持和之前UMAP嵌入一样的参数，只修改reducedDims="Harmony"
projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

#tSNE
projHeme2 <- addTSNE(
    ArchRProj = projHeme2, 
    reducedDims = "Harmony", 
    name = "TSNEHarmony", 
    perplexity = 30
)
p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "Plot-TSNE2Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
#---------------------------------------------------------------------------------------------------------------------------------


# ArchR的基因得分和标记基因
#鉴定标记基因
markersGS <- getMarkerFeatures(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
#SummarizedExperiment对象，里面记录着每个标记特征的相关信息。这是ArchR的一个常见输出，是下游数据分析的关键输出格式。SummarizedExperiment对象和矩阵类似，行是感兴趣的特征（例如基因），列表示样本。一个SummarizedExperiment能够记录一个或多个assay矩阵

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6

markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
    "CD14", "CEBPB", "MPO", #Monocytes
    "IRF8", 
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells
  )

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

#在嵌入上可视化标记基因
markerGenes  <- c(
    "CD34",  #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "MME", #B-Cell Trajectory
    "CD14", "MPO", #Monocytes
    "CD3D", "CD8A"#TCells
  )

p <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

#选择绘图列表的其中一个基因进行展示
p$CD14

#如果是需要绘制所有基因，那么可以使用cowplot将不同的标记基因整合到一张图中

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = projHeme2, 
    addDOC = FALSE, width = 5, height = 5)

#使用MAGIC填充标记基因  因为scATAC-seq数据太过稀疏,基因得分图变化很大
projHeme2 <- addImputeWeights(projHeme2)

# 填充权重会在之后绘制UMAP嵌入图里的基因得分时传入到plotEmbedding()
markerGenes  <- c(
    "CD34",  #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "MME", #B-Cell Trajectory
    "CD14", "MPO", #Monocytes
    "CD3D", "CD8A"#TCells
  )

p <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2)
)
#和之前一样，可以只画其中一个基因

#也可以用cowplot绘制所有的标记基因
# Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

#使用ArchRBrowser绘制Track
# 除了在UMAP上绘制每个细胞的基因得分外，还可以在基因组浏览器上浏览这些标记基因的局部染色体开放状态。使用plotBrowserTrack()函数，它会构建一系列绘图对象，每个对象对应一个标记基因。函数会根据groupBy输入的组别信息在不同track上绘制每组的开放状态。

markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", #B-Cell Trajectory
    "CD14", #Monocytes
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells
  )

p <- plotBrowserTrack(
    ArchRProj = projHeme2, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000
)

# 通过选择列表中的给定基因来绘制最终结果

grid::grid.newpage()
grid::grid.draw(p$CD14)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = projHeme2, 
    addDOC = FALSE, width = 5, height = 5)

#使用scRNA-seq定义cluster类型
# ArchR接受未经修改的Seurat对象作为整合流程的输入。我们使用download.file下载数据

if(!file.exists("scRNA-Hematopoiesis-Granja-2019.rds")){
    download.file(
        url = "https://jeffgranja.s3.amazonaws.com/ArchR/TestData/scRNA-Hematopoiesis-Granja-2019.rds",
        destfile = "scRNA-Hematopoiesis-Granja-2019.rds"
    )
}

seRNA <- readRDS("scRNA-Hematopoiesis-Granja-2019.rds")
seRNA

# class: RangedSummarizedExperiment
# dim: 20287 35582
# metadata(0):
# assays(1): counts
# rownames(20287): FAM138A OR4F5 ... S100B PRMT2
# rowData names(3): gene_name gene_id exonLength
# colnames(35582): CD34_32_R5:AAACCTGAGTATCGAA-1
#   CD34_32_R5:AAACCTGAGTCGTTTG-1 ...
#   BMMC_10x_GREENLEAF_REP2:TTTGTTGCATGTGTCA-1
#   BMMC_10x_GREENLEAF_REP2:TTTGTTGCATTGAAAG-1
# colData names(10): Group nUMI_pre ... BioClassification Barcode

colnames(colData(seRNA))

# 使用table()，我们可以看到scRNA-seq细胞类型每一群的细胞数
table(colData(seRNA)$BioClassification)


#无约束整合 为后续更加精细的约束分析奠定基础
# 使用addGeneIntegrationMatrix()对scATAC-seq和scRNA-seq数据进行整合
projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "BioClassification",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)


#约束整合
cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments


# 首先，我们检查在无约束整合中用到的scRNA-seq数据里细胞类型标签。
unique(unique(projHeme2$predictedGroup_Un))

# 之后从scATAC-seq提取和scRNA-seq对应的聚类

#From scRNA
cTNK <- paste0(paste0(19:25), collapse="|")
cTNK
# [1] "19|20|21|22|23|24|25"

    # 其余的聚类就称之为"Non-T cell, Non-NK cell"
    cNonTNK <- paste0(c(paste0("0", 1:9), 10:13, 15:18), collapse="|")
    cNonTNK
    #[1] "01|02|03|04|05|06|07|08|09|10|11|12|13|15|16|17|18"

# 接着再用字符串模式在preClust找到对应的scATAC-seq列名，然后使用列名从混合矩阵提取对应的列。
# 对于T细胞和NK细胞，scATAC-seq聚类ID就是C7, C8, C9

#Assign scATAC to these categories
clustTNK <- rownames(cM)[grep(cTNK, preClust)]
clustTNK
#[1] "C8" "C9" "C7"
  
    #对于" Non-T cells and Non-NK cells", ID就是scATAC-seq聚类余下的部分
    clustNonTNK <- rownames(cM)[grep(cNonTNK, preClust)]
    clustNonTNK
    # [1] "C10" "C3"  "C11" "C4"  "C1"  "C12" "C5"  "C2"  "C6" 

#接着在scRNA-seq中做相同的操作 筛选出相同的细胞类型。首先，我们鉴定scRNA-seq数据中T细胞和NK细胞
#RNA get cells in these categories
rnaTNK <- colnames(seRNA)[grep(cTNK, colData(seRNA)$BioClassification)]
head(rnaTNK)

    # 然后，鉴定scRNA-seq数据中"Non-T cell Non-NK cell cells"
    rnaNonTNK <- colnames(seRNA)[grep(cNonTNK, colData(seRNA)$BioClassification)]
    head(rnaNonTNK)
    #[1] "CD34_32_R5:AAACCTGAGTATCGAA-1" "CD34_32_R5:AAACCTGAGTCGTTTG-1"
    #[3] "CD34_32_R5:AAACCTGGTTCCACAA-1" "CD34_32_R5:AAACGGGAGCTTCGCG-1"
    #[5] "CD34_32_R5:AAACGGGAGGGAGTAA-1" "CD34_32_R5:AAACGGGAGTTACGGG-1"


groupList <- SimpleList(
    TNK = SimpleList(
        ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustTNK],
        RNA = rnaTNK
    ),
    NonTNK = SimpleList(
        ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustNonTNK],
        RNA = rnaNonTNK
    )    
)


# 将该列表传递给addGeneIntegrationMatrix()函数的groupList参数
projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE, 
    groupList = groupList,
    groupRNA = "BioClassification",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)

#生成调色板
pal <- paletteDiscrete(values = colData(seRNA)$BioClassification)

#在scATAC-seq数据根据无约束整合得到的scRNA-seq细胞类型进行可视化
p1 <- plotEmbedding(
    projHeme2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal
)
p1

# do.call(cowplot::plot_grid, c(list(ncol = 3), p2c))


#约束整合得到scATAC-seq对应的scRNA-seq的细胞类型进行可视化

p2 <- plotEmbedding(
    projHeme2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co", 
    pal = pal
)
p2

plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme2", load = FALSE)

# 为每个scATAC-seq细胞增加拟scRNA-seq谱
#~5 minutes
projHeme3 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "BioClassification",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)

#GeneIntegrationMatrix已经被添加到Arrow文件中 填充权重值(impute weights)
projHeme3 <- addImputeWeights(projHeme3)

#生成一些UMAP图，里面的基因表达量值来自于GeneIntegrationMatrix
markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", #B-Cell Trajectory
    "CD14", #Monocytes
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells
  )

p1 <- plotEmbedding(
    ArchRProj = projHeme3, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    continuousSet = "horizonExtra",
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme3)
)

#相同UMAP图，但是使用GeneScoreMatrix里的基因得分值
p2 <- plotEmbedding(
    ArchRProj = projHeme3, 
    colorBy = "GeneScoreMatrix", 
    continuousSet = "horizonExtra",
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme3)
)


# 最后用cowplot将这些标记基因绘制在一起

p1c <- lapply(p1, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})

p2c <- lapply(p2, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})

do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))

do.call(cowplot::plot_grid, c(list(ncol = 3), p2c))


plotPDF(plotList = p1, 
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", 
    ArchRProj = projHeme3, 
    addDOC = FALSE, width = 5, height = 5)

#使用scRNA-seq信息标记scATAC-seq聚类 确定了scATAC-seq和scRNA-seq数据间的对应关系，我们就可以使用scRNA-seq数据中细胞类型对我们的scATAC-seq聚类进行定义
#首先，我们会在scATAC-seq和整合分析得到predictedGroup之间构建一个混合矩阵
cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup)
labelOld <- rownames(cM)
labelOld

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

# 对每一个scRNA-seq的聚类，重新定义标签，以便更好解释。

remapClust <- c(
    "01_HSC" = "Progenitor",
    "02_Early.Eryth" = "Erythroid",
    "03_Late.Eryth" = "Erythroid",
    "04_Early.Baso" = "Basophil",
    "05_CMP.LMPP" = "Progenitor",
    "06_CLP.1" = "CLP",
    "07_GMP" = "GMP",
    "08_GMP.Neut" = "GMP",
    "09_pDC" = "pDC",
    "10_cDC" = "cDC",
    "11_CD14.Mono.1" = "Mono",
    "12_CD14.Mono.2" = "Mono",
    "13_CD16.Mono" = "Mono",
    "15_CLP.2" = "CLP",
    "16_Pre.B" = "PreB",
    "17_B" = "B",
    "18_Plasma" = "Plasma",
    "19_CD8.N" = "CD8.N",
    "20_CD4.N1" = "CD4.N",
    "21_CD4.N2" = "CD4.N",
    "22_CD4.M" = "CD4.M",
    "23_CD8.EM" = "CD8.EM",
    "24_CD8.CM" = "CD8.CM",
    "25_NK" = "NK"
)
remapClust <- remapClust[names(remapClust) %in% labelNew]

# 接着使用mapLables()函数进行标签转换，将旧的标签映射到新的标签上
labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
labelNew2
# 合并labelsOld和labelsNew2，我们现在可以用mapLables()函数在cellColData里新建聚类标签。
projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew2, oldLabels = labelOld)

p1 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Clusters2")
p1
plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Save-ProjHeme3", load = FALSE)

