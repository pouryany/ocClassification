library(Biobase)
library(GEOquery)
library(stringr)
library(forcats)
library(dplyr)
library(magrittr)
library(caret)
library(pROC)
library(ggbiplot)
library(ClusterR)
library(sva)
library(a4Base)
library(limma)
library(plyr)


rand.label.gen <- function(n=25){
    a <- do.call(paste0,replicate(5,sample(LETTERS,n,TRUE),FALSE))
    return(a)
}
#gene.List <-  sample(tT$ID,25)

# Some list preparation
{
    gene.List <- read.csv(file = "bowenSignals26.csv")
    gene.List <- colnames(gene.List)
    gene.List <- gene.List[-c(1:8)]
    gene.List <- str_sub(gene.List,start = 2)
}

# get 26 gene panel across all data

if(FALSE){
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
        fit2 <- eBayes(fit2, 0.01)
        tT   <- topTable(fit2, number = Inf)
    }
    gene.exp8  <- exprs(gset)
    gene.exp8 <- gene.exp8[row.names(gene.exp8) %in% gene.List,]
    gene.exp8 <- rbind(gene.exp8,sml)
    a <- tail(gene.exp8,1)
    length(a[which(a == "G1")])
    
}
####
#
#
#
#
#
#
#
#
#Normalizing  26 gene panel across all data
# Unifying all genes first.




{
    temp <- read.csv("GeneList.csv")
    all.genes.list <- list(unlist(as.character(temp$x)))
    
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
        gset$batch <- 0
        gset0 <- gset
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- limma::topTable(fit2, number = Inf, adjust="fdr", sort.by= "none")
        
        sum(is.na(tT$adj.P.Val))
        
        all.genes <- rownames(tT)
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
        gset$batch <- 2
        gset2 <- gset[unlist(all.genes.list),]
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT2 <- limma::topTable(fit2, adjust="fdr", sort.by= "none", number=Inf)
        
        nrow(tT2)
        all.genes2 <- rownames(tT2)
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
        gset$batch <- 3
        gset3 <- gset[unlist(all.genes.list),]
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT3 <- limma::topTable(fit2, adjust="fdr", sort.by= "none", number=Inf)
        
        
        all.genes3 <- rownames(tT3)
        
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
        gset$batch <- 4
        gset4 <- gset[unlist(all.genes.list),]
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G2-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT4 <- limma::topTable(fit2, adjust="fdr", sort.by= "none", number=Inf)
        
        
        all.genes4 <- rownames(tT4)
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
    if(FALSE){gset <- getGEO("GSE26712", GSEMatrix =TRUE, AnnotGPL=TRUE)
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
    fit2 <- eBayes(fit2, 0.01)
    tT5 <- topTable(fit2, adjust="fdr", sort.by= "none", number=Inf)
    
    
    all.genes5 <- rownames(tT5)}
    
    #gene.exp5  <- exprs(gset)
    #gene.exp5 <- gene.exp5[row.names(gene.exp5) %in% gene.List,]
    #gene.exp5 <- rbind(gene.exp5,sml)
    #tail(gene.exp5,1)
    
    
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
        gset$batch <- 6
        gset6 <- gset[unlist(all.genes.list),]
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT6 <- limma::topTable(fit2, adjust="fdr", sort.by= "none", number=Inf)
        
        
        all.genes6 <- rownames(tT6)
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
        gset$batch <- 7
        gset7 <- gset[unlist(all.genes.list),]
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT7 <- limma::topTable(fit2, adjust="fdr", sort.by= "none", number=Inf)
        
        
        all.genes7 <- rownames(tT7)}
    
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
        gset$batch <- 8
        gset8 <- gset[unlist(all.genes.list),]
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT8 <- limma::topTable(fit2, adjust="fdr", sort.by= "none", number=Inf)
        
        all.genes8 <- rownames(tT8)
    }
    gene.exp8  <- exprs(gset)
    gene.exp8 <- gene.exp8[row.names(gene.exp8) %in% gene.List,]
    gene.exp8 <- rbind(gene.exp8,sml)
    a <- tail(gene.exp8,1)
    
    ###Removing batch effects
    

    
    #pData(gset0)
    #varMetadata(gset0)
    
    #Getting phenodata right
    {
        p0 <- pData(gset0)
        p0 <- p0[,c("description","batch")]
        
        p2 <- pData(gset2)
        p2 <- p2[,c("description","batch")]
        
        p3 <- pData(gset3)
        p3 <- p3[,c("description","batch")]
        
        p4 <- pData(gset4)
        p4 <- p4[,c("description","batch")]
        
        p6 <- pData(gset6)
        p6 <- p6[,c("description","batch")]
        
        p7 <- pData(gset7)
        p7 <- p7[,c("description","batch")]
        
        p8 <- pData(gset8)
        p8 <- p8[,c("description","batch")]
        
        pdata.all <- list(p0,p2,p3,p4,p6,p7,p8)
        
        pdata.all <- Reduce(rbind,pdata.all)
        
        pdata.all$description  %<>% fct_collapse(G0 = c("G0","G2"))
    }
    # Getting edataright
    {
        edata0 <- exprs(gset0)
        edata2 <- exprs(gset2)
        edata3 <- exprs(gset3)
        edata4 <- exprs(gset4)
        edata6 <- exprs(gset6)
        edata7 <- exprs(gset7)
        edata8 <- exprs(gset8)
        
        edata.all <- list(edata0,edata2,edata3,edata4,edata6,edata7,edata8)
        
        edata.all <- Reduce(cbind,edata.all)
    }
    
    
    
    
    mod0 = model.matrix(~1,data = pdata.all)
    mod  = model.matrix(~description,data = pdata.all)
    
    edata.all <-  edata.all[rowSums(is.na(edata.all)) ==0 ,]
    
    #n.sv  = num.sv(edata.all,mod,method = "leek")
    #svojb = sva(edata.all,mod,mod0, n.sv= n.sv)
    
    
    modcombat = model.matrix(~1, data=pdata.all)
    combat_edata = ComBat(dat=edata.all, batch=pdata.all$batch, mod=modcombat,
                          par.prior=TRUE, prior.plots=FALSE)
    
    
    
    # all.genes.list <- list(all.genes,all.genes2,all.genes3,all.genes4,
    #                       all.genes6,all.genes7,all.genes8)
    #
    #
    # all.genes.list <-  Reduce(intersect,all.genes.list)
    # write.csv(all.genes.list,"GeneList.csv")
    #
    # gset[all.genes.list,]
    
}

#
#
#
#
#
# get a random 26 gene panel across all data
if(FALSE){
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
        fit2 <- eBayes(fit2, 0.01)
        tT   <- topTable(fit2, number = Inf)
    }
    gene.exp8  <- exprs(gset)
    gene.exp8 <- gene.exp8[row.names(gene.exp8) %in% gene.List,]
    gene.exp8 <- rbind(gene.exp8,sml)
    a <- tail(gene.exp8,1)
    length(a[which(a == "G1")])
    
}
#Using just the 26 genes without batch effect removal
if(FALSE){
    gene.exp.tot <- cbind(gene.exp2,gene.exp3,gene.exp4,gene.exp6,gene.exp7,gene.exp8)
    gene.exp.tot <- t(gene.exp.tot)
    gene.exp.tot <- gene.exp.tot[,-1]
    gene.exp.tot <- rbind(gene.exp.tot,t(gene.exp))
    ##gene.exp.tot <- gene.exp.tot[,-1]
    gene.exp.tot <- as.data.frame(gene.exp.tot)
    #gene.exp.tot[,1:25] %<>% mutate_all(.,as.character)%>% mutate_all(.,as.numeric)
    #gene.exp.tot$sml  %<>% fct_collapse(G0 = c("G0","G2"))
    #filter(., sml == "G2") #%>% mutate(., sml = "G0")
    #gene.exp.tot[which(gene.exp.tot$sml== "G2"),]$sml <- "G0"
    #glimpse(gene.exp.tot)
}

#Using SVA for batch effect Removal and subsequent learning
if(FALSE){
    set.seed(12354)
    validation_index <- createDataPartition(pdata.all$description,
                                            p=0.70, list=FALSE)
    trainData  <- edata.all[,validation_index]
    testData   <- edata.all[,-validation_index]
    trainPheno <- pdata.all[validation_index,]
    testPheno  <- pdata.all[-validation_index,]
    
    trainMod   = model.matrix(~description,data=trainPheno)
    trainMod0  = model.matrix(~1,data=trainPheno)
    trainSv    = sva(trainData,trainMod,trainMod0,vfilter=5000)
    fsvaobj    = fsva(trainData,trainMod,trainSv,testData)
    
    train.set  = fsvaobj$db
    #train.set  = as.data.frame(train.set)
    #train.set  = cbind(train.set,as.factor(trainPheno$description))
    tests.set  = fsvaobj$new
    
    train.set <- train.set[row.names(train.set) %in% gene.List,]
    train.set <- as.data.frame(t(train.set))
    train.set <- cbind(train.set,as.factor(trainPheno$description))
    colnames(train.set)[26] <- "sml"
    
    tests.set <- tests.set[row.names(tests.set) %in% gene.List,]
    tests.set <- as.data.frame(t(tests.set))
    tests.set <- cbind(tests.set,as.factor(testPheno$description))
    colnames(tests.set)[26] <- "sml"
    
    gene.exp.tot <- rbind(train.set,tests.set)
    
}



# # # # # # On with some analysis
#           #
#############
#           #
# # # # # # # # # # # # #Data Cleaning

#Using Combat for batch effect Removal
{
    gene.exp.tot <- combat_edata[row.names(combat_edata) %in% gene.List,]
    gene.exp.tot <- as.data.frame(t(gene.exp.tot))
    gene.exp.tot <- cbind(gene.exp.tot,as.factor(pdata.all$description))
    colnames(gene.exp.tot)[26] <- "sml"
}

{
    

    #Using SVA for batch effect Removal and subsequent learning
    if(FALSE){
        set.seed(12354)
        validation_index <- createDataPartition(pdata.all$description,
                                                p=0.70, list=FALSE)
        trainData  <- edata.all[,validation_index]
        testData   <- edata.all[,-validation_index]
        trainPheno <- pdata.all[validation_index,]
        testPheno  <- pdata.all[-validation_index,]
        
        trainMod   = model.matrix(~description,data=trainPheno)
        trainMod0  = model.matrix(~1,data=trainPheno)
        trainSv    = sva(trainData,trainMod,trainMod0,vfilter=5000)
        fsvaobj    = fsva(trainData,trainMod,trainSv,testData)
       
        train.set  = fsvaobj$db
        #train.set  = as.data.frame(train.set)
        #train.set  = cbind(train.set,as.factor(trainPheno$description))
        tests.set  = fsvaobj$new
        
        train.set <- train.set[row.names(train.set) %in% gene.List,]
        train.set <- as.data.frame(t(train.set))
        train.set <- cbind(train.set,as.factor(trainPheno$description))
        colnames(train.set)[26] <- "sml"
        
        tests.set <- tests.set[row.names(tests.set) %in% gene.List,]
        tests.set <- as.data.frame(t(tests.set))
        tests.set <- cbind(tests.set,as.factor(testPheno$description))
        colnames(tests.set)[26] <- "sml"
        
        gene.exp.tot <- rbind(train.set,tests.set)
        
    }
    
    
    
    #ncol(gene.exp.tot)
    #tail(gene.exp.tot,1)
    #sapply(gene.exp.tot,class)
    #sum(gene.exp.tot$sml =="G1")
    
    #Fixing Gene names
    {gene.names <- tT[tT$ID %in% gene.List,]
        gene.names[,c(1,3)]
        
        gene.names <- gene.names[gene.names$ID %in% names(gene.exp.tot),]
        gene.names <- gene.names[,c(1,3)]
        names(gene.exp.tot)
        gene.names <-gene.names[names(gene.exp.tot),]
        gene.names <- gene.names[!is.na(gene.names$ID),]
        gene.names <- gene.names$Gene.symbol
        gene.names <- c(gene.names, "sml")
    }
    
    
    #gene.names[c(8,10,15,22,24)] <- c("aa","bb","cc","dd","ee")
    colnames(gene.exp.tot) <- gene.names
    
    
    
    ########
    #
    #
    #  Preparing data for fitting models
    #
    #
    #
    ########
    
    set.seed(2)
    validation_index <- createDataPartition(gene.exp.tot$sml, p=0.70, list=FALSE)
    
    train.set   <- gene.exp.tot[validation_index,]
    tests.set   <- gene.exp.tot[-validation_index,]
    #which(train.set$sml == "G0")
    #names(train.set)
    input.vals  <- train.set[,1:25]
    output.vals <- train.set[,26]
    #control <- trainControl(method="cv", number=5)
    nrow(tests.set)
    control <- trainControl(method="repeatedcv",number = 5,repeats = 10, summaryFunction = twoClassSummary
                            , classProbs = TRUE)
    #?trainControl
    #metric <- "Accuracy"
    
    
    
    #######
    #
    #
    # Fitting models
    #
    #
    ######
    
  
    
    set.seed(7)
    fit.lda <- train(sml ~ .,train.set, method="lda", metric="ROC", trControl=control)
    probs.lda <- predict(fit.lda,tests.set[,1:25],type = "prob")
    ROC.lda <- roc(tests.set$sml , probs.lda[, "G0"], levels = levels(tests.set$sml))
    #plot(ROC.lda, type = "S", print.thres = .5)
    #plot(fit.lda)
    
    # CART
    set.seed(7)
    fit.cart <- train(sml ~ .,train.set, method="rpart", metric="ROC"
                      ,tuneLength = 15, trControl=control)
    probs.cart <- predict(fit.cart,tests.set[,1:25],type = "prob")
    ROC.cart <- roc(tests.set$sml , probs.cart[, "G0"], levels = levels(tests.set$sml))
    #plot(ROC.cart, type = "S", print.thres = .5)
   
    
    # kNN
    set.seed(7)
    fit.knn <- train(sml ~ .,train.set, method="knn", metric="ROC", trControl=control)
    probs.knn <- predict(fit.knn,tests.set[,1:25],type = "prob")
    ROC.knn <- roc(tests.set$sml , probs.knn[, "G0"], levels = levels(tests.set$sml))
    #plot(ROC.knn, type = "S", print.thres = .2)
    

    # Generalized Linear Model with Stepwise Feature Selection
    set.seed(7)
    fit.glm <- train(sml ~ .,train.set, method = "glmStepAIC", metric="ROC", trControl=control)
    probs.glm <- predict(fit.glm,tests.set[,1:25],type = "prob")
    ROC.glm <- roc(tests.set$sml , probs.glm[, "G0"], levels = levels(tests.set$sml))
    
    
    # Random Forest
    set.seed(7)
    #tunegrid <- expand.grid(.mtry=3)   ###, tuneGrid=tunegrid,importance=TRUE
    fit.rf <- train(sml ~ .,train.set, method = "rf", metric="ROC",
                    trControl=control)
    probs.rf <- predict(fit.rf,tests.set[,1:25],type = "prob")
    ROC.rf <- roc(tests.set$sml , probs.rf[, "G0"], levels = levels(tests.set$sml))
    
    
    # Gradient boosting
    set.seed(7)
    #tunegrid <- expand.grid(.mtry=3)   ###, tuneGrid=tunegrid,importance=TRUE
    fit.gbm <- train(sml ~ .,train.set, method = "gbm", metric="ROC",
                    trControl=control)
    probs.gbm <- predict(fit.gbm,tests.set[,1:25],type = "prob")
    ROC.gbm <- roc(tests.set$sml , probs.gbm[, "G0"], levels = levels(tests.set$sml))
    
    
  
  
    # Neural Networks 
    set.seed(7)
    fit.nnet <- train(sml ~ .,train.set, method = "nnet", metric="ROC",
                     trControl=control)
    probs.nnet <- predict(fit.nnet,tests.set[,1:25],type = "prob")
    ROC.nnet <- roc(tests.set$sml , probs.nnet[, "G0"], levels = levels(tests.set$sml))
    
    
    
      
    # SVM Radial basis 
    set.seed(7)
    SVMgrid <- expand.grid(sigma = c(0.05), C = c(2))
    fit.svm <- train(sml ~ .,train.set, method = "svmRadial", metric="ROC",
                     tuneGrid = SVMgrid,preProc = c("scale","YeoJohnson"),
                     trControl=control)
    probs.svm <- predict(fit.svm,tests.set[,1:25],type = "prob")
    ROC.svm <- roc(tests.set$sml , probs.svm[, "G0"], levels = levels(tests.set$sml))
    
    
    
    plot(ROC.rf, type = "S",col = "#80b1d3",lwd=3,
         cex.lab=2,print.auc = T,print.auc.x = 0.5, print.auc.y= 0.27,print.auc.cex = 1.3)
    #plot(ROC.gbm, add = TRUE, col = "red")
    plot(ROC.lda, add = TRUE, col = "#fb8072",lwd=3
         ,print.auc = T, print.auc.y= 0.21,print.auc.cex = 1.3)
    plot(ROC.knn, add = TRUE, col = "#bebada", lwd=3
         ,print.auc = T,print.auc.x = 0.5, print.auc.y= 0.15,print.auc.cex = 1.3)
    #plot(ROC.glm, add = TRUE, col = "#8dd3c7")
    plot(ROC.svm, add = TRUE, col = "#a65628", lwd=3
         ,print.auc = T,print.auc.x = 0.5, print.auc.y= 0.09,print.auc.cex = 1.3)
    plot(ROC.cart,add = TRUE, col = "#8dd3c7", lwd=3
         ,print.auc = T, print.auc.y= 0.03,print.auc.cex = 1.3)
    legend("bottomright", legend=c("Random Forest", "LDA", "KNN","SVM","CART" ),
           col=c("#80b1d3", "#fb8072","#bebada","#a65628","#8dd3c7"), lwd=4)
    
    
    
    require(randomForest)
    plot(fit.cart)

    
    
   
    rff <- fit.rf$finalModel
    importance(rff,plot = T)
    results <- resamples(list(lda=fit.lda, knn=fit.knn, rf=fit.rf, cart=fit.cart,
                              glm=fit.glm, svm=fit.svm))
    
    

    summary(results)
    dotplot(results)
    
    
    test.res <- predict.train(fit.gbm,tests.set[,1:25])
    a <- caret::confusionMatrix(test.res,tests.set$sml, positive = "G1", mode = "everything")
    a <- a$byClass
    cbind(a,a)
  
    

    
    
    summary(fit.knn)
    
    rf.fit <- train(sml ~ .,train.set, method = "rf")
    results <- resamples(list(lda=fit.lda,  knn=fit.knn, rf=fit.rf))
    
    
    
    
    
    #gene.exp.tot[,1:25] <- as.numeric((gene.exp.tot[,1:25]))
    

    
    
    res.rf  <- predict(fit.rf,tests.set[,1:25])
    res.lda <- predict(fit.lda,tests.set[,1:25])
    res.knn <- predict(fit.knn,tests.set[,1:25])
    res.svm <- predict(fit.svm,tests.set[,1:25])
    res.cart<- predict(fit.cart,tests.set[,1:25])
    
    
    
    
    res.rf  <-  caret::confusionMatrix(res.rf,tests.set$sml, positive = "G1", mode = "everything")
    res.lda <-  caret::confusionMatrix(res.lda,tests.set$sml, positive = "G1", mode = "everything")
    res.knn <-  caret::confusionMatrix(res.knn,tests.set$sml, positive = "G1", mode = "everything")
    res.svm <-  caret::confusionMatrix(res.svm,tests.set$sml, positive = "G1", mode = "everything")
    res.cart<-  caret::confusionMatrix(res.cart,tests.set$sml, positive = "G1", mode = "everything")
    
  
    a <- list( "Random Forest" = res.rf$byClass, "LDA" = res.lda$byClass,
               "KNN" = res.knn$byClass,"SVM" = res.svm$byClass ,"CART" = res.cart$byClass)
    a <- Reduce(rbind, a)
    rownames(a) <-c("Random Forest", "LDA", "KNN","SVM","CART" )
    a <- t(a)
    a <- a[-c(5,6,8,9,10),]
    
    library(xtable)
    options(xtable.floating = FALSE)
    options(xtable.timestamp = "")
    
    print(xtable(a), include.rownames = T)
    
    
    plot(res.rf)
    
    disp.features <- input.vals[,dispnames]
    disp.vals     <- revalue(output.vals, c("G0"="Non-Cancer", "G1"="Cancer"))
    featurePlot(x= disp.features, y= disp.vals, plot = "box",cex.lab = 30)
    plot(train.set)
    
}



#####
#####
#####
#####
#####
# Some random testing 
#####
b <- data_frame("rands")
for (i in 1:50){   
    #Preparing a random set of data
    
    set.seed(i)
    rand.index   <- sample(nrow(combat_edata),25,replace = F)
    gene.exp.tot <- combat_edata[rand.index,]
    gene.exp.tot <- as.data.frame(t(gene.exp.tot))
    gene.exp.tot <- cbind(gene.exp.tot,as.factor(pdata.all$description))
    #Fixing Gene names
    gene.names   <- rand.label.gen(n = 25)
    gene.names   <- c(gene.names, "sml")
    colnames(gene.exp.tot) <- gene.names
    

    set.seed(2)
    
    validation_index <- createDataPartition(gene.exp.tot$sml, p=0.70, list=FALSE)
    train.set        <- gene.exp.tot[validation_index,]
    tests.set        <- gene.exp.tot[-validation_index,]
    input.vals       <- train.set[,1:25]
    output.vals      <- train.set[,26]
    control          <- trainControl(method="repeatedcv",number= 5,repeats = 5,
                                     summaryFunction = twoClassSummary,
                                     classProbs = TRUE)

    # Random Forest
    set.seed(7)
    fit.rf   <- train(sml ~ .,train.set, method = "rf", metric="ROC",
                    trControl=control)
    probs.rf <- predict(fit.rf,tests.set[,1:25],type = "prob")
    ROC.rf   <- roc(tests.set$sml, probs.rf[, "G0"], levels=levels(tests.set$sml))
    #results  <- resamples(list(rf = fit.rf))
    plot(ROC.rf, type = "S", print.thres = .5)
    
    
    test.res <- predict.train(fit.rf,tests.set[,1:25])
    a <- caret::confusionMatrix(test.res,tests.set$sml, positive = "G1",
                                mode = "everything")
    a$overall
    a <- a$byClass
    b <- cbind(b,a)
    print(paste0("\n done",i))
}


tests.set$sml == ""
b <- read.csv("randomMeasure.csv")
b <- b[,-c(1)]
b <- as.data.frame(b)
b <- rowMeans(b)

names(res.rf$byClass)


final.results <- as.data.frame(cbind(res.rf$byClass,b))
final.results <- final.results[-c(5,6,8,9,10),]
colnames(final.results) <- c("RF Select 25 genes", "RF Random 25 genes Average")
print(xtable(final.results), include.rownames = T)





names(final.results) <-c('Proposed.method','Random.Features')
final.results <- round(100*final.results , 1)

fwrite(final.results , 'final.results.csv')




#####
#
#
#  Some plotting with PCA
#
#####
{
    gene.exp.tot <- edata.all
    totSignals <- t(as.matrix(gene.exp.tot))
    #totSignals <- as.matrix(gene.exp.tot[,-26])
    #totType    <- (gene.exp.tot[,26])
    totType    <- (pdata.all$description)
    
    totType    <- revalue(totType, c("G0"="Non-Cancer", "G1"="Cancer"))
    pcapdata   <- pdata.all
    pcapdata[pcapdata$batch ==0,]$batch <- 1
    batchType  <- as.factor(pcapdata$batch)
    batchType  <-  revalue(batchType, c("1"="GSE12172", "2"="GSE14407","3"="GSE9899",
                                        "4"="GSE37648","6"="GSE18521","7"="GSE38666",
                                        "8"="GSE10971"))
    
    #totSignals <- totSignals[nrow(totSignals):1,]
    #totType    <-  totType[length(totType):1]
    #batchType  <- batchType[length(batchType):1]
    
    Tissue.Type <- totType
    class(totSignals) <- "numeric"
    tot.pca <- prcomp(totSignals, center = TRUE, scale. = TRUE)
    #print(tot.pca)
    #dim(tot.pca)
    #dim(tot.pca[[5]])
    #tot.pca$x <- tot.pca$x[,2:20]
    #tot.pca$rotation <- tot.pca$rotation[,2:20]
    #tot.pca$sdev <- tot.pca$sdev[2:20]
    g <- ggbiplot(tot.pca , obs.scale = 1, var.scale =1 , var.axes = F, groups = batchType, ellipse = T)
    #g <- g+ ggbiplot(tot.pca , obs.scale = 1, var.scale =1 , var.axes = F,ellipse.prob = 0.95, groups = totType, ellipse = T)
    #g <-  g+geom_point(aes(shape=as.factor(modTot$cluster) , colour = totType ), size = 3)
    g <- g + scale_color_discrete(name = "Batch") + geom_point(size=3,aes(color =  batchType,shape = Tissue.Type))
    g <- g+ theme(legend.direction = "vertical", 
            legend.position = "right",
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = 16),
            legend.title=element_text(face = "bold", size = 14),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    print(g)
}

design <- model.matrix(~ sml+ 0, as.data.frame(gene.exp.tot))
colnames(design) <- levels(pdata.all$description)
fit <- lmFit( as.data.frame(t(gene.exp.tot[,1:25])), design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- limma::topTable(fit2, number = Inf, adjust="fdr", sort.by= "P")
tTdisp <- tT[,c(2,4)]
print(xtable(tTdisp), include.rownames = T)

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

b
