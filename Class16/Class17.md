# Class17
Yuxuan Jiang A17324184

- [Downstream analysis](#downstream-analysis)
  - [PCA](#pca)
- [OPTIONAL: Differential-expression
  analysis](#optional-differential-expression-analysis)

## Downstream analysis

``` r
library(tximport)
library(rhdf5)
# setup the folder and filenames to read
folders <- dir(pattern="SRR21568*")
samples <- sub("_quant", "", folders)
files <- file.path( folders, "abundance.h5" )
names(files) <- samples
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
```

    1 2 3 4 

``` r
head(txi.kallisto$counts)
```

                    SRR2156848 SRR2156849 SRR2156850 SRR2156851
    ENST00000539570          0          0    0.00000          0
    ENST00000576455          0          0    2.62037          0
    ENST00000510508          0          0    0.00000          0
    ENST00000474471          0          1    1.00000          0
    ENST00000381700          0          0    0.00000          0
    ENST00000445946          0          0    0.00000          0

``` r
colSums(txi.kallisto$counts)
```

    SRR2156848 SRR2156849 SRR2156850 SRR2156851 
       2563611    2600800    2372309    2111474 

``` r
sum(rowSums(txi.kallisto$counts)>0)
```

    [1] 94561

``` r
to.keep <- rowSums(txi.kallisto$counts) > 0
kset.nonzero <- txi.kallisto$counts[to.keep,]
keep2 <- apply(kset.nonzero,1,sd)>0
x <- kset.nonzero[keep2,]
```

### PCA

``` r
pca <- prcomp(t(x), scale=TRUE)
summary(pca)
```

    Importance of components:
                                PC1      PC2      PC3   PC4
    Standard deviation     183.6379 177.3605 171.3020 1e+00
    Proportion of Variance   0.3568   0.3328   0.3104 1e-05
    Cumulative Proportion    0.3568   0.6895   1.0000 1e+00

``` r
plot(pca$x[,1],pca$x[,2],col=c("blue","blue","red","red"),xlab="PC1",ylab="PC2",pch=16)
```

![](Class17_files/figure-commonmark/unnamed-chunk-7-1.png)

#### ggplot PC1 vs PC2

``` r
library(ggrepel)
```

    Loading required package: ggplot2

``` r
colData <- data.frame(condition = factor(rep(c("control", "treatment"), each = 2)))
rownames(colData) <- colnames(txi.kallisto$counts)
y <- as.data.frame(pca$x)
y$Condition <- as.factor(colData$condition)

library(ggplot2)
ggplot(y)+
  aes(PC1,PC2,col=Condition)+
  labs(x="PC1",y="PC2",title="PC1 vs PC2 ")+
  geom_point()+
  geom_text_repel(label=rownames(y))+
  theme_bw()
```

![](Class17_files/figure-commonmark/unnamed-chunk-8-1.png)

#### ggplot PC1 vs PC3

``` r
library(ggplot2)
ggplot(y)+
  aes(PC1,PC3,col=Condition)+
  labs(x="PC1",y="PC3",title="PC1 vs PC3 ")+
  geom_point()+
  geom_text_repel(label=rownames(y))+
  theme_bw()
```

![](Class17_files/figure-commonmark/unnamed-chunk-9-1.png)

#### ggplot PC2 vs PC3

``` r
library(ggplot2)
ggplot(y)+
  aes(PC2,PC3,col=Condition)+
  labs(x="PC2",y="PC3",title="PC2 vs PC3 ")+
  geom_point()+
  geom_text_repel(label=rownames(y))+
  theme_bw()
```

![](Class17_files/figure-commonmark/unnamed-chunk-10-1.png)

## OPTIONAL: Differential-expression analysis

``` r
library(DESeq2)
sampleTable <- data.frame(condition = factor(rep(c("control", "treatment"), each = 2)))
rownames(sampleTable) <- colnames(txi.kallisto$counts)
dds <- DESeqDataSetFromTximport(txi.kallisto,
                                sampleTable, 
                                ~condition)
dds <- DESeq(dds)
res <- results(dds)
head(res)
```

    log2 fold change (MLE): condition treatment vs control 
    Wald test p-value: condition treatment vs control 
    DataFrame with 6 rows and 6 columns
                     baseMean log2FoldChange     lfcSE      stat    pvalue
                    <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENST00000539570  0.000000             NA        NA        NA        NA
    ENST00000576455  0.761453       3.155061   4.86052 0.6491203  0.516261
    ENST00000510508  0.000000             NA        NA        NA        NA
    ENST00000474471  0.484938       0.181923   4.24871 0.0428185  0.965846
    ENST00000381700  0.000000             NA        NA        NA        NA
    ENST00000445946  0.000000             NA        NA        NA        NA
                         padj
                    <numeric>
    ENST00000539570        NA
    ENST00000576455        NA
    ENST00000510508        NA
    ENST00000474471        NA
    ENST00000381700        NA
    ENST00000445946        NA
