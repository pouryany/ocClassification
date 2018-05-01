
library(Biobase)
library(GEOquery)
library(limma)

getwd()
gene.List <- read.csv(file = "bowenSignals26.csv")
gene.List <- colnames(gene.List)
gene.List <- gene.List[-c(1:8)]
library(stringr)

gene.List <- str_sub(gene.List,start = 2)


# get 24 gene panel across all data
{
# load series and platform data from GEO
# 30 ovarian LMP vs 60 ovarian cancer "GSE12172"
{
gset <- getGEO("GSE12172", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("11111010110111001100111111110101001010101111101111",
               "1101111010011110011110000011111010010111")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl

design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by= "none", number=Inf)

head(tT)

}



gene.exp  <- exprs(gset)
gene.exp <- gene.exp[row.names(gene.exp) %in% gene.List,]
gene.exp <- rbind(gene.exp,sml)
tail(gene.exp)
#
#
#
#
#
#
#
#
#
## 24 Samples from bowen "GSE14407"
# load series and platform data from GEO
{
gset <- getGEO("GSE14407", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "000000000000111111111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by= "none", number=Inf)

}

gene.exp2  <- exprs(gset)
gene.exp2 <- gene.exp2[row.names(gene.exp2) %in% gene.List,]
gene.exp2 <- rbind(gene.exp2,sml)
tail(gene.exp2)

#
#
#
#
#
#
#
#
#


### GSE 9899 data

{
gset <- getGEO("GSE9899", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("11111111110000000000000000001111111111111111111111",
               "111X11111111111111111111X111111111XX11111111111111",
               "11111111111111111111111111111111111111111111111111",
               "X11111111111X1X1111111X1111111X11111XX1111111X1111",
               "11X1111111111111X1111111111111111111111111X11111X1",
               "111111111X1111111X11X1111111X111111111X111111")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

}

gene.exp3  <- exprs(gset)
gene.exp3 <- gene.exp3[row.names(gene.exp3) %in% gene.List,]
gene.exp3 <- rbind(gene.exp3,sml)
tail(gene.exp3,1)
#
#
#
#
#
#
#
#
#


### 6 normal ovarian and 6 normal fallopian

# load series and platform data from GEO "GSE37648"
{
gset <- getGEO("GSE37648", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "222000222000"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G2-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
}
gene.exp4  <- exprs(gset)
gene.exp4 <- gene.exp4[row.names(gene.exp4) %in% gene.List,]
gene.exp4 <- rbind(gene.exp4,sml)
tail(gene.exp4,1)
#
#
#
#
#
#
#
#
# 185 malignancy and 10 normal from GSE26712
{gset <- getGEO("GSE26712", GSEMatrix =TRUE, AnnotGPL=TRUE)
    if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]

    # make proper column names to match toptable
    fvarLabels(gset) <- make.names(fvarLabels(gset))

    # group names for all samples
    gsms <- paste0("00000000001111111111111111111111111111111111111111",
                   "11111111111111111111111111111111111111111111111111",
                   "11111111111111111111111111111111111111111111111111",
                   "111111111111111111111111111111111111111111111")
    sml <- c()
    for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

    # log2 transform
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) }

    # set up the data and proceed with analysis
    sml <- paste("G", sml, sep="")    # set group names
    fl <- as.factor(sml)
    gset$description <- fl
    design <- model.matrix(~ description + 0, gset)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(G0-G1, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)}

gene.exp5  <- exprs(gset)
gene.exp5 <- gene.exp5[row.names(gene.exp5) %in% gene.List,]
gene.exp5 <- rbind(gene.exp5,sml)
tail(gene.exp5,1)


#
#
#
#
#
#
#


### GSE6008  contains a variety of samples

#
#
#
#
#
#
#
##
#
#
#
# GSE18521 contain 10 normals and a lot of cancer
{
    gset <- getGEO("GSE18521", GSEMatrix =TRUE, AnnotGPL=TRUE)
    if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]

    # make proper column names to match toptable
    fvarLabels(gset) <- make.names(fvarLabels(gset))

    # group names for all samples
    gsms <- paste0("XXXXXXXXXXXX11111111111111111111111111111111111111",
                   "1111111111111110000000000")
    sml <- c()
    for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

    # eliminate samples marked as "X"
    sel <- which(sml != "X")
    sml <- sml[sel]
    gset <- gset[ ,sel]

    # log2 transform
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) }

    # set up the data and proceed with analysis
    sml <- paste("G", sml, sep="")    # set group names
    fl <- as.factor(sml)
    gset$description <- fl
    design <- model.matrix(~ description + 0, gset)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(G1-G0, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
}

gene.exp6  <- exprs(gset)
gene.exp6 <- gene.exp6[row.names(gene.exp6) %in% gene.List,]
gene.exp6 <- rbind(gene.exp6,sml)
tail(gene.exp6)
#

#
#
#
#
#
#
#
#
#
#
# This is a fantastic sample related to our LCM Samples
{gset <- getGEO("GSE38666", GSEMatrix =TRUE, AnnotGPL=TRUE)
    if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]

    # make proper column names to match toptable
    fvarLabels(gset) <- make.names(fvarLabels(gset))

    # group names for all samples
    gsms <- "XXXXXXXX000000000000XXXXXXX111111111111111111"
    sml <- c()
    for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

    # eliminate samples marked as "X"
    sel <- which(sml != "X")
    sml <- sml[sel]
    gset <- gset[ ,sel]

    # log2 transform
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) }

    # set up the data and proceed with analysis
    sml <- paste("G", sml, sep="")    # set group names
    fl <- as.factor(sml)
    gset$description <- fl
    design <- model.matrix(~ description + 0, gset)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(G1-G0, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)}

gene.exp7  <- exprs(gset)
gene.exp7 <- gene.exp7[row.names(gene.exp7) %in% gene.List,]
gene.exp7 <- rbind(gene.exp7,sml)
tail(gene.exp7)
#
#
#
#
#
#
#
#
#
#
# Another fantastic data that contains LCM of OC and fallopian tubes
{gset <- getGEO("GSE10971", GSEMatrix =TRUE, AnnotGPL=TRUE)
    if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]

    # make proper column names to match toptable
    fvarLabels(gset) <- make.names(fvarLabels(gset))

    # group names for all samples
    gsms <- "0000000000000000000000001111111111111"
    sml <- c()
    for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

    # log2 transform
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) }

    # set up the data and proceed with analysis
    sml <- paste("G", sml, sep="")    # set group names
    fl <- as.factor(sml)
    gset$description <- fl
    design <- model.matrix(~ description + 0, gset)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(G1-G0, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)}
gene.exp8  <- exprs(gset)
gene.exp8 <- gene.exp8[row.names(gene.exp8) %in% gene.List,]
gene.exp8 <- rbind(gene.exp8,sml)
tail(gene.exp)

}


# # # # # # On with some analysis
#           #
#############
#           #
# # # # # # # # # # # # #Data Cleaning
library(forcats)
library(dplyr)
library(magrittr)
library(caret)

gene.exp.tot <- cbind(gene.exp2,gene.exp3,gene.exp4,gene.exp6,gene.exp7,gene.exp8)
gene.exp.tot <- t(gene.exp.tot)
gene.exp.tot <- gene.exp.tot[,-1]
gene.exp.tot <- rbind(gene.exp.tot,t(gene.exp))
gene.exp.tot <- as.data.frame(gene.exp.tot)
sapply(gene.exp.tot,class)


gene.exp.tot[,1:25] %<>% mutate_all(.,as.character)%>% mutate_all(.,as.numeric)

gene.exp.tot$sml  %<>% fct_collapse(G0 = c("G0","G2"))
#filter(., sml == "G2") #%>% mutate(., sml = "G0")
#gene.exp.tot[which(gene.exp.tot$sml== "G2"),]$sml <- "G0"

glimpse(gene.exp.tot)

validation_index <- createDataPartition(gene.exp.tot$sml, p=0.80, list=FALSE)

train.set   <- gene.exp.tot[validation_index,]
tests.set   <- gene.exp.tot[-validation_index,]

input.vals  <- train.set[,1:25]
output.vals <- train.set[,26]
#gene.exp.tot[,1:25] <- as.numeric((gene.exp.tot[,1:25]))

featurePlot(x= input.vals, y= output.vals, plot = "box")
plot(train.set)
library(ggbiplot)
library(ClusterR)


totSignals <- gene.exp.tot[,1:25]
totType    <- as.factor(gene.exp.tot[,26])
class(totSignals) <- "numeric"
tot.pca <- prcomp(totSignals, center = TRUE, scale. = TRUE)
print(tot.pca)

g <- ggbiplot(tot.pca , obs.scale = 1, var.scale =1 , var.axes = F, groups = totType, ellipse = T)
#g <- g+ ggbiplot(tot.pca , obs.scale = 1, var.scale =1 , var.axes = F,ellipse.prob = 0.95, groups = totType, ellipse = T)
#g <-  g+geom_point(aes(shape=as.factor(modTot$cluster) , colour = totType ), size = 3)
g <- g + scale_color_discrete(name = "Sample Types")
g <- g+ theme(legend.direction = "vertical", legend.position = "right")
print(g)





### From GSE 4122, it's in an old machinery
{# load series and platform data from GEO

    gset <- getGEO("GSE4122", GSEMatrix =TRUE, AnnotGPL=TRUE)
    if (length(gset) > 1) idx <- grep("GPL201", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]

    # make proper column names to match toptable
    fvarLabels(gset) <- make.names(fvarLabels(gset))

    # group names for all samples
    gsms <- "1111110121X10012102121020001202020120120200111X211211111012X1101011"
    sml <- c()
    for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

    # eliminate samples marked as "X"
    sel <- which(sml != "X")
    sml <- sml[sel]
    gset <- gset[ ,sel]

    # log2 transform
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) }

    # set up the data and proceed with analysis
    sml <- paste("G", sml, sep="")    # set group names
    fl <- as.factor(sml)
    gset$description <- fl
    design <- model.matrix(~ description + 0, gset)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(G2-G0, G1-G0, G2-G1, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)}
#
#
#
#
#
#
#
#
##
#
#
##
#
#
#
#


