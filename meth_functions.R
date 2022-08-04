suppressPackageStartupMessages({
  library("plyr")
  library("R.utils")
  library("missMethyl")
  library("limma")
  library("DMRcate")
  library("DMRcatedata")
  library("topconfects")
  library("minfi")
  library("IlluminaHumanMethylation450kmanifest")
  library("RColorBrewer")
  library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  library("GEOquery")
  library("eulerr")
  library("plyr")
  library("gplots")
  library("reshape2")
  library("forestplot")
  library("beeswarm")
  library("RCircos")
  library("qqman")
  library("ENmix")
})


# scree plot shows the amount of variation in a dataset that is accounted
# for by the first N principal components
myscree <- function(mx,n=10,main="") {
  pc <- princomp(mx)$sdev
  pcp <- pc/sum(pc)*100
  pcp <- pcp[1:10]
  barplot(pcp,cex.names = 1,las=2,ylim=c(0,60),
      ylab="percent (%) variance explained", main=main)
  text((0.5:length(pcp)*1.2),pcp,label=signif(pcp,3),pos=3,cex=0.8)
}


# Here is a function to make a volcano plot
make_volcano <- function(dm,name,mx) {
    sig <- subset(dm,adj.P.Val<0.05)
    N_SIG=nrow(sig)
    N_UP=nrow(subset(sig,logFC>0))
    N_DN=nrow(subset(sig,logFC<0))
    HEADER=paste(N_SIG,"@5%FDR,", N_UP, "up", N_DN, "dn")
    plot(dm$logFC,-log10(dm$P.Val),cex=0.5,pch=19,col="darkgray",
        main=name, xlab="log FC", ylab="-log10 pval")
    mtext(HEADER)
    grid()
    points(sig$logFC,-log10(sig$P.Val),cex=0.5,pch=19,col="red")
}

# Metal volcano
volcano_metal <- function(metal,name) {
  sig <- subset(metal,FDR<0.05)
  N_SIG=nrow(sig)
  N_UP=nrow(subset(sig,Zscore>0))
  N_DN=nrow(subset(sig,Zscore<0))
  HEADER=paste(N_SIG,"@5%FDR,", N_UP, "up", N_DN, "dn")
  plot(metal$Zscore,-log10(metal$P.value),cex=0.5,pch=19,col="darkgray",
       main=name, xlab="log FC", ylab="-log10 pval")
  mtext(HEADER)
  grid()
  points(sig$Zscore,-log10(sig$P.value),cex=0.5,pch=19,col="red")
}

# Here is a function to make heatmaps based on smallest p-values
make_heatmap <- function(dm,name,mx,n, groups) {
  topgenes <-  rownames(head(dm[order(dm$P.Value),],n))
  ss <- mx[which(rownames(mx) %in% topgenes),]
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
  colCols <- as.numeric(as.factor(groups))
  colCols <- gsub("1","orange",colCols)
  colCols <- gsub("0","yellow",colCols)
  heatmap.2(ss,scale="row",margin=c(10, 10),cexRow=0.4,trace="none",cexCol=0.4,
    ColSideColors=colCols ,  col=my_palette, main=name)
}

# Here is a function to make heatmaps based on topconfects
make_heatmap2 <- function(confects,name,mx,n, groups) {
  topgenes <-  head(confects$table$name,n)
  ss <- mx[which(rownames(mx) %in% topgenes),]
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
  colCols <- as.numeric(as.factor(groups))
  colCols <- gsub("1","orange",colCols)
  colCols <- gsub("0","yellow",colCols)
  heatmap.2(ss,scale="row",margin=c(10, 10),cexRow=0.4,trace="none",cexCol=0.4,
    ColSideColors=colCols ,  col=my_palette, main=name)
}


# make beeswarm charts
# dm = a limma differential meth object
# name = character name of the limma dm object
# mx = matrix of normalised data
# groups = a vector of factors corresponding to the cols in mx
# n = the number of top significant genes to plot (default = 15) 
make_beeswarms <- function(dm,name,mx,groups,n=15) {
    par(mar=c(3,3,1,1))
    NCOLS=5
    NROWS=floor(n/NCOLS)
    if (n %% NCOLS > 0) { NROWS <- NROWS + 1 }
    par(mfrow=c(NROWS, NCOLS))
    topgenes <-  rownames(head(dm[order(dm$P.Value),],n))
    ss <- mx[which(rownames(mx) %in% topgenes),]
    n <- 1:n
    g1name=levels(groups)[1]
    g2name=levels(groups)[2]
    g1dat <- ss[n,which(groups == g1name)]
    g2dat <- ss[n,which(groups == g2name)]
    g1l <-lapply(split(g1dat, row.names(g1dat)), unlist)
    g2l <-lapply(split(g2dat, row.names(g2dat)), unlist)

    for (i in n) {
      mydat <- list(g1l[[i]],g2l[[i]])
        beeswarm(mydat,ylim=c(0,1),cex=0.2, pch=19,
        las=2, cex.lab=0.6, main=names( g1l )[i] , 
        ylab="",labels = c(g1name,g2name))
      grid()
    }
}

# make beeswarm charts for best confects
# dm = a limma differential meth object
# name = character name of the limma dm object
# mx = matrix of normalised data
# groups = a vector of factors corresponding to the cols in mx
# n = the number of top significant genes to plot (default = 15) 
make_beeswarms_confects <- function(confects,name,mx,groups,n=15) {
    par(mar=c(3,3,1,1))
    NCOLS=5
    NROWS=floor(n/NCOLS)
    if (n %% NCOLS > 0) { NROWS <- NROWS + 1 }
    par(mfrow=c(NROWS, NCOLS))
    topgenes <-  head(confects$table,n)$name
    ss <- mx[which(rownames(mx) %in% topgenes),]
    n <- 1:n
    g1name=levels(groups)[1]
    g2name=levels(groups)[2]
    g1dat <- ss[n,which(groups == g1name)]
    g2dat <- ss[n,which(groups == g2name)]
    g1l <-lapply(split(g1dat, row.names(g1dat)), unlist)
    g2l <-lapply(split(g2dat, row.names(g2dat)), unlist)

    for (i in n) {
      mydat <- list(g1l[[i]],g2l[[i]])
        beeswarm(mydat,ylim=c(0,1),cex=0.2, pch=19,
        las=2, cex.lab=0.6, main=names( g1l )[i] , 
        ylab="",labels = c(g1name,g2name))
      grid()
    }
}

## Rcircos
make_circos <- function(dmr) {
  par(mfrow=c(1,1))
  dmrbed <- data.frame(seqnames=seqnames(dmr),
                       starts=as.integer(start(dmr)-1),
                       ends=as.integer(end(dmr)),
                       names=as.character(elementMetadata(dmr)[,8]),
                       scores=elementMetadata(dmr)[,7],
                       strands=strand(dmr))
  
  dmrbed_up <- subset(dmrbed,scores>0)
  dmrbed_dn <- subset(dmrbed,scores<0)  
  
  # lets try a circos plot
  data(UCSC.HG19.Human.CytoBandIdeogram)
  #chr.exclude <- c("chrX","chrY")
  chr.exclude <- NULL
  cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
  tracks.inside <- 1
  tracks.outside <- 0
  RCircos.Set.Core.Components(cyto.info, chr.exclude,
                              tracks.inside, tracks.outside)
  par(mai=c(0.5, 0.5, 0.5, 0.5))
  RCircos.Set.Plot.Area(margins = 0)

  track.num <- 1;
  side <- "in";
  data.col <- 5;
  RCircos.Tile.Plot(dmrbed_up, track.num, side)
  # heatmap dont work
  #RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col=5, track.num, side)
  #RCircos.Heatmap.Plot(dmrbed, data.col=5, track.num, side)

  track.num <- 2;
  side <- "in";
  data.col <- 5;
  RCircos.Tile.Plot(dmrbed_dn, track.num, side)

  RCircos.Chromosome.Ideogram.Plot(tick.interval = 10000)

}


# this is a wrapper which creates three charts
# We will be adding more
make_dm_plots <- function(dm,name,mx,mxs,groups=groups,confects=confects,dmr,comp=comp,cgi=cgi) {
    make_volcano(dm,name,mx)
    make_beeswarms(dm ,name , mx , groups , n= 15)
    make_heatmap(dm , name , mxs ,n = 50, groups)
    make_beeswarms_confects(confects, name, mx, groups, n=15)
    make_heatmap2(confects, name, mxs, n = 50, groups)
    
    dm_up <- rownames(subset(dm,adj.P.Val<0.05 & logFC>0))
    dm_dn <- rownames(subset(dm,adj.P.Val<0.05 & logFC<0))
    sig <- min(length(dm_up),length(dm_dn))
    if (sig>10) {
      make_forest_plots(comp)
      make_forest_plots(cgi)
      make_circos(dmr)
    }  
}  

make_forest_plots <- function(comp) {

  comp_data <- 
    structure(list(
      "mean"  = comp$up_comp$OR , 
      "lower" = comp$up_comp$lowerCI ,
      "upper" = comp$up_comp$upperCI ,
      .Names = c("mean", "lower", "upper"), 
      row.names = c(NA, -11L), 
      class = "data.frame"))

  comp_data <- as.data.frame(comp_data[1:3],row.names = rownames(comp$up_comp) )

  forestplot(comp_data,title = "hypermethylated",
    labeltext = as.list(rownames(comp_data)),
    mean=mean,lower=lower,upper=upper)

comp_data <- 
  structure(list(
    "mean"  = comp$dn_comp$OR , 
    "lower" = comp$dn_comp$lowerCI ,
    "upper" = comp$dn_comp$upperCI ,
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame"))

comp_data <- as.data.frame(comp_data[1:3],row.names = rownames(comp$dn_comp) )

  forestplot(comp_data,title = "hypomethylated",
    labeltext = as.list(rownames(comp_data)),
    mean=mean,lower=lower,upper=upper)

}

# this is a function which will perform differential methylation analysis
# if you provide it with the right inputs
dm_analysis <- function(samplesheet,sex,groups,mx,name,myann,beta) {
    design <- model.matrix(~ sex + groups)
    mxs <- mx[,which( colnames(mx) %in% samplesheet$Basename )]
    fit.reduced <- lmFit(mxs,design)
    fit.reduced <- eBayes(fit.reduced)
    summary(decideTests(fit.reduced))
    dm <- topTable(fit.reduced,coef=3, number = Inf)
    dma <- merge(myann,dm,by=0)
    dma <- dma[order(dma$P.Value),]
    dm_up <- rownames(subset(dm,adj.P.Val<0.05 & logFC>0))
    dm_dn <- rownames(subset(dm,adj.P.Val<0.05 & logFC<0))
    sig <- min(length(dm_up),length(dm_dn))
    confects <- limma_confects(fit.reduced, coef=3, fdr=0.05)
    if (sig > 10) {
      dmr <- run_dmrcate(mx=mxs,design=design) 
      head(dmr)
      comp <- compartment_enrichment(dma)
      comp
      cgi <- cgi_enrichment(dma)
      cgi
    } else {
      dmr <- NULL
      comp <- NULL
      cgi <- NULL
    }
    make_dm_plots(dm = dm ,name=name , mx=beta,mxs=mxs, groups = groups, 
      confects=confects,dmr = dmr, comp=comp, cgi=cgi)
    dat <- list("dma"=dma, "dm_up"=dm_up, "dm_dn"=dm_dn, "confects"=confects,
      "dmr"= dmr, "comp"=comp, "cgi"=cgi, "fit"=fit.reduced)
    return(dat)
}

# this function performs DMRcate for peak calling
run_dmrcate <- function(mx,design) {
  fit.reduced <- lmFit(mx,design)
  fit.reduced <- eBayes(fit.reduced)
  
  dmg <- makeGenomicRatioSetFromMatrix(mx, rownames = NULL, pData = NULL ,
                                       array = "IlluminaHumanMethylation450k" ,
                                       mergeManifest = FALSE, what = "M")
  
  myannotation <- cpg.annotate("array", dmg, arraytype = "450K",
                               analysis.type="differential", design=design, coef=3)
  
  dmrcoutput <- dmrcate(myannotation, lambda=1000, C=3)
  
  dmr <- extractRanges(dmrcoutput, genome = "hg19")
  
  return(dmr)
}

compartment_enrichment <- function(dma) {
  up <- subset(dma,logFC>0 & adj.P.Val<0.05)
  dn <- subset(dma,logFC<0 & adj.P.Val<0.05)
  all <- table(unique(dma)$Regulatory_Feature_Group)
  up <- table(unique(up)$Regulatory_Feature_Group)
  dn <- table(unique(dn)$Regulatory_Feature_Group)
  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(up,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  rownames(xx)[1] <- "Intergenic"
  xx[,1] = NULL
  colnames(xx) <- c("all","up")
  xx[is.na(xx)] <- 0
  head(xx)
  x=xx$up
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$up)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$up)-x[2], sum(xx$all) - sum(xx$up) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  up_comp <- xx
  
  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(dn,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  rownames(xx)[1] <- "Intergenic"
  xx[,1] = NULL
  colnames(xx) <- c("all","dn")
  xx[is.na(xx)] <- 0
  x=xx$dn
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$dn)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$dn)-x[2], sum(xx$all) - sum(xx$dn) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  dn_comp <- xx
  list("up_comp"=up_comp,"dn_comp"=dn_comp)
}


cgi_enrichment <- function(dma) {
  up <- subset(dma,logFC>0 & adj.P.Val<0.05)
  dn <- subset(dma,logFC<0 & adj.P.Val<0.05)
  all <- table(unique(dma)$Relation_to_Island)
  up <- table(unique(up)$Relation_to_Island)
  dn <- table(unique(dn)$Relation_to_Island)
  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(up,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  xx[,1] = NULL
  colnames(xx) <- c("all","up")
  xx[is.na(xx)] <- 0
  head(xx)
  x=xx$up
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$up)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$up)-x[2], sum(xx$all) - sum(xx$up) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  up_comp <- xx
  
  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(dn,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  xx[,1] = NULL
  colnames(xx) <- c("all","dn")
  xx[is.na(xx)] <- 0
  x=xx$dn
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$dn)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$dn)-x[2], sum(xx$all) - sum(xx$dn) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  dn_comp <- xx
  list("up_comp"=up_comp,"dn_comp"=dn_comp)
}


# Spearman ranks
myranks <- function(x) {
  xx<-x[[1]]
  xx$score <- sign(xx$logFC)/log10(xx$adj.P.Val)
  y <- xx[,"score",drop=FALSE]
  y$rn <- xx$Row.names
  return(y)
}

# enrichment analysis functions
run_mitch_rank <-function(dma){
  dmap <- dma[grep("Promoter_Associated",dma$Regulatory_Feature_Group),]
  dmap[which(dmap$UCSC_RefGene_Name==""),2] <- "NA"
  dmap$genename <- sapply(strsplit(dmap$UCSC_RefGene_Name,";"),"[[",1)
  dmap2 <- dmap[,c("genename","t")]
  #rank <- aggregate(. ~ genename,dmap2,function(x) { 
  #  x[which.max(abs(x))]
  #  } )
  rank <- aggregate(. ~ genename,dmap2,mean)
  rownames(rank) <- rank$genename
  rank$genename=NULL
  return(rank)
}

run_mitch_1d <- function(dma,name) {
  library("mitch")
  rank <- run_mitch_rank(dma)
  capture.output(
    res <- mitch_calc(x = rank,genesets = genesets, priority = "significance",resrows=20)
    , file = "/dev/null", append = FALSE,
    type = c("output", "message"), split = FALSE)
  head(res$enrichment_result,20)
  capture.output(
    mitch_plots(res,outfile=paste(name,".pdf",sep=""))
    , file = "/dev/null", append = FALSE,
    type = c("output", "message"), split = FALSE)
  return(res$enrichment_result)
}


calc_sc <- function(dm) {
  gn <- unique(unlist(strsplit( dm$UCSC_RefGene_Name ,";")))
  gnl <- strsplit( dm$UCSC_RefGene_Name ,";")
  gnl <- mclapply(gnl,unique,mc.cores=CORES)
  dm$UCSC_RefGene_Name <- gnl
  l <- mclapply(1:nrow(dm), function(i) {
    a <- dm[i,]
    len <- length(a[[1]][[1]])
    tvals <- as.numeric(rep(a[2],len))
    genes <- a[[1]][[1]]
    data.frame(genes,tvals)
  },mc.cores=CORES)
  df <- do.call(rbind,l)
  gme_res <- mclapply( 1:length(gn), function(i) {
    g <- gn[i]
    tstats <- df[which(df$genes==g),"tvals"]
    myn <- length(tstats)
    mymean <- mean(tstats)
    mymedian <- median(tstats)
    wtselfcont <- wilcox.test(tstats)
    res <- c("gene"=g,"nprobes"=myn,"mean"=mymean,"median"=mymedian,
      "p-value(sc)"=wtselfcont$p.value)
  } , mc.cores=CORES )
  gme_res_df <- do.call(rbind, gme_res)
  rownames(gme_res_df) <- gme_res_df[,1]
  gme_res_df <- gme_res_df[,-1]
  tmp <- apply(gme_res_df,2,as.numeric)
  rownames(tmp) <- rownames(gme_res_df)
  gme_res_df <- as.data.frame(tmp)
  gme_res_df$sig <- -log10(gme_res_df[,4])
  gme_res_df <- gme_res_df[order(-gme_res_df$sig),]
  gme_res_df$`fdr(sc)` <- p.adjust(gme_res_df$`p-value(sc)`)
  out <- list("df"=df,"gme_res_df"=gme_res_df)
  return(out)
}

gmea_volc <- function(res) {
  sig <- subset(res,`fdr(sc)` < 0.05)
  plot(res$median , -log10(res$`p-value(sc)`) ,
    xlab="effect size (mean t-stat)", ylab="-log10(p-value)",
    pch=19, cex=0.5, col="gray",main="self contained test")
  grid()
  points(sig$median , -log10(sig$`p-value(sc)`) ,
    pch=19, cex=0.5, col="red")
}

gmea_boxplot <- function(res) {
  df <- res[[1]]
  res <- res[[2]]
  par(mfrow=c(1,2))
  n=50
  gs <- head(rownames(res),50)
  tstats <- lapply(gs, function(g) {
    df[which(df$genes==g),"tvals"]
  })
  names(tstats) <- gs
  tstats <- tstats[order(unlist(lapply(tstats,median)))]
  boxplot(tstats,horizontal=TRUE,las=1,
    main="smallest p-val(selfcont)",cex.axis=0.6,
    xlab="t-statistic")
  grid()
  sig <- subset(res,`fdr(sc)` < 0.05)
  gs <- head(rownames(sig[order(-abs(sig$median)),]),n)
  if ( length(gs) >2 ) {
    tstats <- lapply(gs, function(g) {
      df[which(df$genes==g),"tvals"]
    })
    names(tstats) <- gs
    tstats <- tstats[order(unlist(lapply(tstats,median)))]
    boxplot(tstats,horizontal=TRUE,las=1,
      main="biggest effect size(median)",cex.axis=0.6,
      xlab="t-statistic")
    grid()
  } else {
    plot(1)
    mtext("too few significant genes found")
  }
  par(mfrow=c(1,1))
}
