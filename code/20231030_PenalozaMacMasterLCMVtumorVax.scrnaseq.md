LCMV scRNA-Seq processing
================
Slim FOURATI
2024-02-09

Load required packages

``` r
suppressPackageStartupMessages(library(package = "knitr"))
suppressPackageStartupMessages(library(package = "hdf5r"))
suppressPackageStartupMessages(library(package = "Seurat"))
suppressPackageStartupMessages(library(package = "celldex"))
suppressPackageStartupMessages(library(package = "SingleR"))
suppressPackageStartupMessages(library(package = "BiocParallel"))
suppressPackageStartupMessages(library(package = "biomaRt"))
suppressPackageStartupMessages(library(package = "ggbeeswarm"))
suppressPackageStartupMessages(library(package = "ggpubr"))
suppressPackageStartupMessages(library(package = "ggrepel"))
suppressPackageStartupMessages(library(package = "xlsx"))
suppressPackageStartupMessages(library(package = "tidyverse"))
```

``` r
opts_chunk$set(echo = TRUE, fig.path = "../figure/")
options(readr.show_col_types   = FALSE,
        dplyr.summarise.inform = FALSE)
workDir <- dirname(getwd())
```

``` r
seqFiles <- list.files(path       = file.path(workDir, "input"), 
                      full.names = TRUE,
                      pattern = ".+LCMV.+h5$")
seuratObj <- NULL
for (seqFile in seqFiles) {
    seuratTemp <- Read10X_h5(filename = seqFile)
    sampleId <- gsub(pattern     = "\\..+",
                     replacement = "",
                     basename(seqFile))
    seuratTemp <- CreateSeuratObject(seuratTemp, 
                                     project = sampleId)
    seuratTemp <- RenameCells(seuratTemp,
                              new.names = paste0(sampleId,
                                                 "_",
                                                 colnames(seuratTemp)))
    seuratTemp <- DietSeurat(seuratTemp)

    if (is.null(seuratObj)) {
        seuratObj <- seuratTemp
    } else {
        seuratObj <- merge(x = seuratObj, y = seuratTemp)
    }
}
```

``` r
rm(seuratTemp)
```

# Quality control

Percentage of mitochondrial reads

``` r
ensembl <- useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl")
gene2chr <- getBM(attributes = c("mgi_symbol", "chromosome_name"), 
                  filters = "mgi_symbol", 
                  values = rownames(seuratObj$RNA), 
                  mart = ensembl)

mito.genes <- filter(gene2chr, chromosome_name %in% "MT") %>%
  .$mgi_symbol

percent.mito <- Matrix::colSums(seuratObj$RNA@counts[mito.genes, ])/
  Matrix::colSums(seuratObj$RNA@counts)

# AddMetaData adds columns to object@meta.data, and is a great place to
seuratObj <- AddMetaData(object   = seuratObj,
                         metadata = percent.mito,
                         col.name = "percent.mito")
```

``` r
ggplot(data    = seuratObj@meta.data,
       mapping = aes(x = orig.ident, y = percent.mito)) +
  geom_boxplot() +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percentage of reads that are mitochondrial") +
  theme_bw()
```

![](../figure/plot-mito-1.png)<!-- -->

Percent of ribosomal reads

``` r
# look at ribosomal genes
ribo.genes <- grep(pattern = "^Rps|^Rpl", 
                   rownames(x = seuratObj$RNA@counts), 
                   value   = TRUE)
percent.ribo <- Matrix::colSums(seuratObj$RNA@counts[ribo.genes, ])/Matrix::colSums(seuratObj$RNA@counts)

# AddMetaData adds columns to object@meta.data, and is a great place to
seuratObj <- AddMetaData(object   = seuratObj,
                         metadata = percent.ribo,
                         col.name = "percent.ribo")
```

``` r
ggplot(data    = seuratObj@meta.data,
       mapping = aes(x = orig.ident, y = percent.ribo)) +
  geom_boxplot() +
  theme_bw()
```

![](../figure/plot-ribo-1.png)<!-- -->

Number of cell detected

``` r
nbCellDF <- table(seuratObj@meta.data$orig.ident) %>%
  as.data.frame() %>%
  rename(orig.ident                  = Var1,
         `Estimated Number of Cells` = Freq)
meanReadsPerCellDF <- colSums(seuratObj$RNA@counts) %>%
  data.frame(eta = .) %>%
  rownames_to_column() %>%
  mutate(orig.ident = seuratObj@meta.data$orig.ident) %>%
  group_by(orig.ident) %>%
  summarize(`Mean Reads per Cell` = mean(eta))
medianGenesPerCell <- colSums(seuratObj$RNA@counts > 0) %>%
  data.frame(eta = .) %>%
  rownames_to_column() %>%
  mutate(orig.ident = seuratObj@meta.data$orig.ident) %>%
  group_by(orig.ident) %>%
  summarize(`Median Genes per Cell` = median(eta))

plotDF <- merge(x    = nbCellDF,
                y    = meanReadsPerCellDF,
                by   = "orig.ident") %>%
  merge(y  = medianGenesPerCell,
        by = "orig.ident") %>%
  pivot_longer(cols = -orig.ident)

ggplot(data = plotDF,
       mapping = aes(x = orig.ident, y = value)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = ~name, scale = "free", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

<img src="../figure/read-qc-ge-1.png" style="display: block; margin: auto auto auto 0;" />

``` r
plotDF %>%
  group_by(name) %>%
  summarize(median  = median(value),
            min     = min(value),
            max     = max(value)) %>%
  rename(metric = name) %>%
  kable()
```

| metric                    |   median |      min |      max |
|:--------------------------|---------:|---------:|---------:|
| Estimated Number of Cells |  6595.00 |  5243.00 |  7947.00 |
| Mean Reads per Cell       | 11294.36 | 10698.85 | 11889.88 |
| Median Genes per Cell     |  2110.00 |  2051.00 |  2169.00 |

# Cell annotation

``` r
DefaultAssay(seuratObj) <- "RNA"
seuratObj <- NormalizeData(seuratObj) %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()

seuratObj <- RunUMAP(seuratObj, dims = 1:10, n.components = 2L) 
```

Expression of canonical markers

``` r
feat1 <- FeaturePlot(object   = seuratObj,
                     features = c("Ptprc", "Tyrp1"),
                     raster   = FALSE)
feat1
```

<img src="../figure/feature-plot-cd45-1.png" style="display: block; margin: auto auto auto 0;" />
Ptprc=Cd45; cells clustering on UMAP1 \> 5 are tumors cells (CD45-).
Tyrp1 is a marker of melanocytes.

``` r
feat2 <- FeaturePlot(object   = seuratObj, 
                         features = c("Cd3g", "Cd3g", "Cd3e",
                                        "Cd4", "Cd8a", "Cd8b1"),
                     raster   = FALSE)
feat2
```

![](../figure/feature-plot-t-1.png)<!-- -->

``` r
feat3 <- FeaturePlot(object = seuratObj,
                     features = c("Cd14",
                                  "Fcgr3",
                                  "Cd33",
                                  "Itgax" #CD11c
                                  ),
  raster = FALSE)
feat3
```

![](../figure/feature-plot-immune-1.png)<!-- -->

``` r
feat4 <- FeaturePlot(object = seuratObj,
                     features = c("Cd19",
                                  "Ms4a1", # CD20
                                  "Cd79a",
                                  "Cd79b"),
  raster = FALSE)
feat4
```

![](../figure/feature-plot-immune-2.png)<!-- -->

``` r
feat5 <- FeaturePlot(object = seuratObj,
                     features = c("LCMV-chrL", "LCMV-chrS"),
  raster = FALSE,
  split.by = "orig.ident")
feat5
```

<img src="../figure/feature-plot-lcmv-1.png" style="display: block; margin: auto auto auto 0;" />

SingleR at the cell levels

``` r
scaledMat <- seuratObj$RNA@data
immgen <- celldex::ImmGenData()
predSubset <- SingleR(test    = scaledMat,
                      ref     = immgen,
                      labels  = immgen$label.main,
                      BPPARAM = MulticoreParam(workers = 7))
seuratObj <- AddMetaData(object   = seuratObj,
                         metadata = predSubset$labels,
                         col.name = "immgen.main")

umap1 <- DimPlot(seuratObj,
reduction  = 'umap',
group.by   = 'immgen.main',
raster     = FALSE) +
theme(legend.position = "bottom",
      legend.key.size = unit(x = 0.1, units = "in"),
      legend.text = element_text(size = 7))
umap1
```

<img src="../figure/singler-1.png" style="display: block; margin: auto auto auto 0;" />

SingleR at the cluster level

``` r
seuratObj <- FindNeighbors(seuratObj, reduction = "umap", dims = 1:2)
seuratObj <- FindClusters(seuratObj, resolution = 0.5)
```

``` r
umap2 <- DimPlot(seuratObj,
        reduction  = 'umap', 
        group.by   = 'seurat_clusters',
        raster     = FALSE,
        label = TRUE) + 
  theme(legend.position = "bottom")
umap2
```

<img src="../figure/umap-cluster-1.png" style="display: block; margin: auto auto auto 0;" />

``` r
scaledMat <- seuratObj$RNA@data
immgen <- celldex::ImmGenData()
predSubset <- SingleR(test    = scaledMat,
                      ref     = immgen,
                      clusters = seuratObj@meta.data$seurat_clusters,
                      labels  = immgen$label.main,
                      BPPARAM = MulticoreParam(workers = 7))
# correct cluster annotation
immgen.main <- predSubset$pruned.labels[match(seuratObj$seurat_clusters,
                        table = rownames(predSubset))]
immgen.main[seuratObj$seurat_clusters %in% c(12, 15, 27)] <- "Tumor"
immgen.main[seuratObj$seurat_clusters %in% 9] <- "Stroma"
immgen.main[seuratObj$seurat_clusters %in% 25] <- "DC"
immgen.main[seuratObj$seurat_clusters %in% c(3, 4, 8, 10, 13)] <- "Monocytes/Macrophages"
immgen.main[seuratObj$seurat_clusters %in% c(0, 11, 18, 21)] <- "T cells"
immgen.main[seuratObj$seurat_clusters %in% c(5, 6, 16, 22)] <- "NK cells"

if ("immgen.main" %in% colnames(seuratObj@meta.data)) {
  seuratObj$immgen.main <- NULL
}
seuratObj <- AddMetaData(object   = seuratObj,
                         metadata = immgen.main,
                         col.name = "immgen.main")
```

``` r
umap3 <- DimPlot(seuratObj,
        reduction  = 'umap', 
        group.by   = 'immgen.main',
        raster     = FALSE) + 
  theme(legend.position = "bottom")
umap3
```

![](../figure/umap-singler-cluster-1.png)<!-- -->

``` r
seuratObj$orig.ident <- gsub(pattern     = "_LCMV",
                             replacement = "",
                             seuratObj$orig.ident)
umap4 <- DimPlot(seuratObj,
                 reduction  = 'umap',
                 group.by   = 'immgen.main',
                 split.by   = 'orig.ident',
                 raster     = FALSE) +
  theme(legend.position = "bottom")
umap4 +
  geom_point(data = umap4$data[(Matrix::colSums(seuratObj$RNA@counts[c("LCMV-chrL", "LCMV-chrS"), ]) >= 1), ],
             mapping = aes(x = UMAP_1, y = UMAP_2), cex = 0.5) +
  labs(title = NULL)
```

![](../figure/umap-by-treat-1.png)<!-- -->

``` r
save(seuratObj, file = file.path(workDir, "output/lcmv.seuratObj.RData"))
```

# LCMV reads

``` r
vln2 <- VlnPlot(object = seuratObj,
          features = c("LCMV-chrL", "LCMV-chrS"),
        group.by= "immgen.main",
        split.by = "orig.ident")
vln2
```

![](../figure/lcmv-vln-1.png)<!-- --> LCMV reads are only detected in
the R3 samples (right) and mostly in tumor cells and then in some
monocytes and DCs.

# Type I ISGs expression

``` r
isgLS <- fgsea:::gmtPathways("/Users/iew5629/Desktop/Projects/Utils/MSigDB/h.all.v2023.1.Hs.symbols.gmt") %>%
  .[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]]

human <- useMart(biomart = "ensembl", 
                 dataset = "hsapiens_gene_ensembl",
                 host    = "https://dec2021.archive.ensembl.org/")
mouse <- useMart(biomart = "ensembl", 
                 dataset="mmusculus_gene_ensembl",
                 host    = "https://dec2021.archive.ensembl.org/")
human2mouse <- getLDS(mart = human, attributes = "hgnc_symbol", 
                      filters = "hgnc_symbol", values = isgLS, 
                      attributesL = "mgi_symbol", martL = mouse)
seuratObj <- AddModuleScore(seuratObj, features = list(isg = human2mouse$MGI.symbol),
                            name = "isg")
```

``` r
FeaturePlot(seuratObj, 
            features = "isg1",
            split.by = "orig.ident")
```

![](../figure/umap-is-1.png)<!-- -->

``` r
vln1 <- VlnPlot(object = seuratObj,
          features = "isg1",
        group.by= "orig.ident",
        split.by = "immgen.main")

vln1$data %>%
  ggplot(mapping = aes(x = ident, y = isg1)) +
  geom_beeswarm(size = 0.5, cex = 0.75) +
  geom_boxplot(outlier.color = "transparent", fill = "transparent", color = "red") +
  labs(x = NULL, y = "Type I ISGs") +
  stat_compare_means(cex = 2) +
  facet_wrap(facets = ~split, scale = "free_y", nrow = 2) +
  theme_bw()
```

![](../figure/boxplot-isg-1.png)<!-- -->

``` r
tumorObj <- subset(seuratObj, immgen.main == "Tumor")
tumorObj <- subset(tumorObj, seurat_clusters != "9") # low Tyrp1 (melanocyte marker)
Matrix::colSums(tumorObj$RNA@counts[c("LCMV-chrS", "LCMV-chrL"), ]) %>%
  data.frame(LCMV = .) %>%
  rownames_to_column(var = "cellbarcode") %>%
  merge(y = rownames_to_column(tumorObj@meta.data, var = "cellbarcode"),
        by = "cellbarcode") %>%
  mutate(`LCMV > 0` = (LCMV > 0),
         `LCMV > 0` = paste0(gsub(pattern = "_.+",
                                  replacement = "",
                                  orig.ident),
                             ".",
                             `LCMV > 0`)) %>%
  ggplot(mapping = aes(x = `LCMV > 0`, y = isg1)) +
  geom_beeswarm(cex = 0.75) + 
  geom_boxplot(outlier.colour = "transparent", color = "red", fill = "transparent") +
  labs(title = "Tumor cells", y = "Type I ISGs") +
  stat_compare_means(ref.group = "PBS.FALSE") +
  theme_bw()
```

<img src="../figure/cor-tumor-lcmv-isg-1.png" style="display: block; margin: auto auto auto 0;" />
Tumor cells with LCMV reads express higher levels of Type I ISG than
tumor cells without LCMV reads or tumor cells from the PBS.

``` r
freqDF <- seuratObj@meta.data %>%
  group_by(orig.ident, immgen.main) %>%
  summarize(n = n(), .groups = "drop")

seuratObj@meta.data %>%
  group_by(orig.ident) %>%
  summarize(tot = n(), .groups = "drop") %>%
  merge(x = freqDF, by = "orig.ident") %>%
  mutate(freq = n/tot) -> freqDF

ggplot(data = freqDF,
       mapping = aes(x = orig.ident, y = freq * 100)) +
  geom_bar(stat = "identity", mapping = aes(fill = immgen.main)) +
  scale_y_continuous() +
  labs(y = "Frequency (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](../figure/test-prop-1.png)<!-- -->

``` r
Idents(seuratObj) <- interaction(seuratObj$immgen.main,
                                 seuratObj$orig.ident,
                                 drop = TRUE)
fits <- NULL
for (GOI1 in grep(pattern = "R3$", levels(seuratObj), value = TRUE)) {
  GOI2 <- gsub(pattern = "R3$", replacement = "PBS", GOI1)
  de.markers <- FindMarkers(seuratObj, ident.1 = GOI1, ident.2 = GOI2) %>%
    mutate(contrast = paste0(GOI1, "-", GOI2))
  fits[[paste0(GOI1, "-", GOI2)]] <- de.markers
}

save(fits, file = file.path(workDir, "output/lcmv.fits.RData"))
```

``` r
tumorObj <- subset(seuratObj, immgen.main == "Tumor")

lcmv <- ifelse(test = Matrix::colSums(tumorObj$RNA@counts[c("LCMV-chrL", "LCMV-chrS"), ]) >= 1,
                                           yes = "LCMVpos",
                                           no  = "LCMVneg")
tumorObj <- AddMetaData(tumorObj, 
                        metadata = lcmv,
                        col.name = "LCMV")
Idents(tumorObj) <- tumorObj$LCMV
set.seed(seed = 1)
de.markers <- FindMarkers(tumorObj, 
                          test.use = "negbinom",
                          ident.1 = "LCMVpos",
                          ident.2 = "LCMVneg",
                          logfc.threshold = 0)
```

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

``` r
ggplot(data = de.markers,
       mapping = aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_point() +
  geom_point(data = rownames_to_column(de.markers, var = "gene_name") %>%filter(gene_name %in% human2mouse$MGI.symbol & p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25),
             color = "red") +
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0, 15)) +
  geom_hline(yintercept = -1 * log10(filter(de.markers, p_val_adj <= 0.05) %>% .$p_val %>% max()), linetype = 2) +
  geom_text_repel(data = rownames_to_column(de.markers, var = "gene_name") %>%filter(gene_name %in% human2mouse$MGI.symbol & p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25),
                  mapping = aes(label = gene_name), color = "red") +
  theme_bw()
```

    ## Warning: ggrepel: 13 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](../figure/volcano-isg-tumor-1.png)<!-- -->

``` r
tumorR3Obj <- subset(seuratObj, immgen.main == "Tumor" & orig.ident %in% "R3")

lcmv <- ifelse(test = Matrix::colSums(tumorR3Obj$RNA@counts[c("LCMV-chrL", "LCMV-chrS"), ]) >= 1,
                                           yes = "LCMVpos",
                                           no  = "LCMVneg")
tumorR3Obj <- AddMetaData(tumorR3Obj, 
                        metadata = lcmv,
                        col.name = "LCMV")
Idents(tumorR3Obj) <- tumorR3Obj$LCMV
set.seed(seed = 1)
de.markers <- FindMarkers(tumorR3Obj, 
                          test.use = "negbinom",
                          ident.1 = "LCMVpos",
                          ident.2 = "LCMVneg",
                          logfc.threshold = 0)
```

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in sqrt(1/i): NaNs produced

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

``` r
ggplot(data = de.markers,
       mapping = aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_point() +
  geom_point(data = rownames_to_column(de.markers, var = "gene_name") %>%filter(gene_name %in% human2mouse$MGI.symbol & p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25),
             color = "red") +
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0, 20)) +
  geom_hline(yintercept = -1 * log10(filter(de.markers, p_val_adj <= 0.05) %>% .$p_val %>% max()), linetype = 2) +
  geom_text_repel(data = rownames_to_column(de.markers, var = "gene_name") %>%filter(gene_name %in% human2mouse$MGI.symbol & p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25),
                  mapping = aes(label = gene_name), color = "red") +
  theme_bw()
```

    ## Warning: ggrepel: 37 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](../figure/volcano-isg-tumor-2.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.3.2 (2023-10-31)
    ## Platform: aarch64-apple-darwin23.0.0 (64-bit)
    ## Running under: macOS Sonoma 14.3
    ## 
    ## Matrix products: default
    ## BLAS:   /opt/homebrew/Cellar/openblas/0.3.26/lib/libopenblasp-r0.3.26.dylib 
    ## LAPACK: /opt/homebrew/Cellar/r/4.3.2/lib/R/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Chicago
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] lubridate_1.9.3             forcats_1.0.0              
    ##  [3] stringr_1.5.1               dplyr_1.1.4                
    ##  [5] purrr_1.0.2                 readr_2.1.5                
    ##  [7] tidyr_1.3.1                 tibble_3.2.1               
    ##  [9] tidyverse_2.0.0             xlsx_0.6.5                 
    ## [11] ggrepel_0.9.5               ggpubr_0.6.0               
    ## [13] ggbeeswarm_0.7.2            ggplot2_3.4.4              
    ## [15] biomaRt_2.58.2              BiocParallel_1.36.0        
    ## [17] SingleR_2.4.1               celldex_1.12.0             
    ## [19] SummarizedExperiment_1.32.0 Biobase_2.62.0             
    ## [21] GenomicRanges_1.54.1        GenomeInfoDb_1.38.5        
    ## [23] IRanges_2.36.0              S4Vectors_0.40.2           
    ## [25] BiocGenerics_0.48.1         MatrixGenerics_1.14.0      
    ## [27] matrixStats_1.2.0           SeuratObject_5.0.1         
    ## [29] Seurat_4.4.0                hdf5r_1.3.9                
    ## [31] knitr_1.45                 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RcppAnnoy_0.0.22              splines_4.3.2                
    ##   [3] later_1.3.2                   bitops_1.0-7                 
    ##   [5] filelock_1.0.3                polyclip_1.10-6              
    ##   [7] XML_3.99-0.16.1               lifecycle_1.0.4              
    ##   [9] rstatix_0.7.2                 globals_0.16.2               
    ##  [11] lattice_0.22-5                MASS_7.3-60.0.1              
    ##  [13] backports_1.4.1               magrittr_2.0.3               
    ##  [15] plotly_4.10.4                 rmarkdown_2.25               
    ##  [17] yaml_2.3.8                    httpuv_1.6.14                
    ##  [19] sctransform_0.4.1             spam_2.10-0                  
    ##  [21] sp_2.1-3                      spatstat.sparse_3.0-3        
    ##  [23] reticulate_1.35.0             cowplot_1.1.3                
    ##  [25] pbapply_1.7-2                 DBI_1.2.1                    
    ##  [27] RColorBrewer_1.1-3            abind_1.4-5                  
    ##  [29] zlibbioc_1.48.0               Rtsne_0.17                   
    ##  [31] RCurl_1.98-1.14               xlsxjars_0.6.1               
    ##  [33] rappdirs_0.3.3                GenomeInfoDbData_1.2.11      
    ##  [35] irlba_2.3.5.1                 listenv_0.9.1                
    ##  [37] spatstat.utils_3.0-4          goftest_1.2-3                
    ##  [39] spatstat.random_3.2-2         fitdistrplus_1.1-11          
    ##  [41] parallelly_1.36.0             DelayedMatrixStats_1.24.0    
    ##  [43] leiden_0.4.3.1                codetools_0.2-19             
    ##  [45] DelayedArray_0.28.0           xml2_1.3.6                   
    ##  [47] tidyselect_1.2.0              farver_2.1.1                 
    ##  [49] ScaledMatrix_1.10.0           BiocFileCache_2.10.1         
    ##  [51] spatstat.explore_3.2-6        jsonlite_1.8.8               
    ##  [53] ellipsis_0.3.2                progressr_0.14.0             
    ##  [55] ggridges_0.5.6                survival_3.5-7               
    ##  [57] progress_1.2.3                tools_4.3.2                  
    ##  [59] ica_1.0-3                     Rcpp_1.0.12                  
    ##  [61] glue_1.7.0                    gridExtra_2.3                
    ##  [63] SparseArray_1.2.3             xfun_0.41                    
    ##  [65] withr_3.0.0                   BiocManager_1.30.22          
    ##  [67] fastmap_1.1.1                 fansi_1.0.6                  
    ##  [69] rsvd_1.0.5                    digest_0.6.34                
    ##  [71] timechange_0.3.0              R6_2.5.1                     
    ##  [73] mime_0.12                     colorspace_2.1-0             
    ##  [75] scattermore_1.2               tensor_1.5                   
    ##  [77] spatstat.data_3.0-4           RSQLite_2.3.5                
    ##  [79] utf8_1.2.4                    generics_0.1.3               
    ##  [81] data.table_1.15.0             prettyunits_1.2.0            
    ##  [83] httr_1.4.7                    htmlwidgets_1.6.4            
    ##  [85] S4Arrays_1.2.0                uwot_0.1.16                  
    ##  [87] pkgconfig_2.0.3               rJava_1.0-11                 
    ##  [89] gtable_0.3.4                  blob_1.2.4                   
    ##  [91] lmtest_0.9-40                 XVector_0.42.0               
    ##  [93] htmltools_0.5.7               carData_3.0-5                
    ##  [95] fgsea_1.28.0                  dotCall64_1.1-1              
    ##  [97] scales_1.3.0                  png_0.1-8                    
    ##  [99] rstudioapi_0.15.0             tzdb_0.4.0                   
    ## [101] reshape2_1.4.4                nlme_3.1-164                 
    ## [103] curl_5.2.0                    zoo_1.8-12                   
    ## [105] cachem_1.0.8                  BiocVersion_3.18.1           
    ## [107] KernSmooth_2.23-22            vipor_0.4.7                  
    ## [109] parallel_4.3.2                miniUI_0.1.1.1               
    ## [111] AnnotationDbi_1.64.1          pillar_1.9.0                 
    ## [113] grid_4.3.2                    vctrs_0.6.5                  
    ## [115] RANN_2.6.1                    promises_1.2.1               
    ## [117] car_3.1-2                     BiocSingular_1.18.0          
    ## [119] beachmat_2.18.0               dbplyr_2.4.0                 
    ## [121] xtable_1.8-4                  cluster_2.1.6                
    ## [123] beeswarm_0.4.0                evaluate_0.23                
    ## [125] cli_3.6.2                     compiler_4.3.2               
    ## [127] rlang_1.1.3                   crayon_1.5.2                 
    ## [129] ggsignif_0.6.4                future.apply_1.11.1          
    ## [131] labeling_0.4.3                plyr_1.8.9                   
    ## [133] stringi_1.8.3                 viridisLite_0.4.2            
    ## [135] deldir_2.0-2                  munsell_0.5.0                
    ## [137] Biostrings_2.70.2             lazyeval_0.2.2               
    ## [139] spatstat.geom_3.2-8           Matrix_1.6-5                 
    ## [141] ExperimentHub_2.10.0          hms_1.1.3                    
    ## [143] patchwork_1.2.0               sparseMatrixStats_1.14.0     
    ## [145] bit64_4.0.5                   future_1.33.1                
    ## [147] KEGGREST_1.42.0               shiny_1.8.0                  
    ## [149] highr_0.10                    interactiveDisplayBase_1.40.0
    ## [151] AnnotationHub_3.10.0          ROCR_1.0-11                  
    ## [153] broom_1.0.5                   igraph_2.0.1.1               
    ## [155] memoise_2.0.1                 fastmatch_1.1-4              
    ## [157] bit_4.0.5
