
randomGeneSets <- function(gene_catalog, lengths, seed){
  set.seed(seed) ; gene_catalog_half <- sample(gene_catalog,length(gene_catalog)/2)
  num_gsets <- length(lengths)
  set.seed(seed) ; seeds <- sample(1:1e6, num_gsets)
  gsets <- lapply(1:num_gsets,function(i) {
    set.seed(seeds[i]) ; gs <- sample(gene_catalog_half,lengths[i])
    return(gs)
  } )
  names(gsets)<-stringi::stri_rand_strings(length(gsets), 15, pattern = "[A-Za-z]")
  return(gsets)
}

incorp_dm <- function(genesets,myann,mval,seed, frac_genes,frac_probes,
  groupsize,delta=1,gene_catalog) {

  # divide gene sets between hyper and hypomethylated
  nset <- floor(length(genesets)/2)
  set.seed(seed) ; gtup <-sample(genesets,nset)
  set.seed(seed) ; gtdn <- sample(setdiff(genesets,gtup),nset)
  gup <- unname(unlist(gtup))
  gdn <- unname(unlist(gtdn))
  # add extra 10% DM genes
  gnon <- setdiff(gene_catalog,c(gup,gdn))
  gextra <- round(length(gnon)*0.1)
  set.seed(seed) ; gup <- c(gup,sample(gnon,gextra))
  gnon <- setdiff(gene_catalog,c(gup,gdn))
  set.seed(seed) ; gdn <- c(gdn,sample(gnon,gextra))
  # make probe-gene vector
  probe2gene <- strsplit(myann$UCSC_RefGene_Name,";")
  names(probe2gene) <- rownames(myann)
  probe2gene <- unlist(probe2gene)
  # select probes hypermethylated
  set.seed(seed) ; gup2 <- sample(gup,floor(length(gup)*frac_genes))
  pup <- names(probe2gene[which(probe2gene %in% gup2)])
  set.seed(seed) ; pup2 <- sample(pup,floor(length(pup)*frac_probes))
  # select probes hypomethylated
  set.seed(seed) ; gdn2 <- sample(gdn,floor(length(gdn)*frac_genes))
  pdn <- names(probe2gene[which(probe2gene %in% gdn2)])
  set.seed(seed) ; pdn2 <- sample(pdn,floor(length(pdn)*frac_probes))
  # add 10% DM probes as well
  probes <- rownames(myann)
  pnon <- setdiff(probes,c(pup,pdn))
  pextra <- round(length(pnon)*0.1)
  set.seed(seed) ; pup <- c(pup,sample(pnon,pextra))
  pnon <- setdiff(probes,c(pup,pdn))
  set.seed(seed) ; pdn <- c(pdn,sample(pnon,pextra))
  # divide samples between ctrl and case
  ncols <- ncol(mval)
  maxgroupsize=floor(ncols/2)
  if ( groupsize > maxgroupsize ) { stop("groupsize cannot be larger than half the ncols of mval") }
  set.seed(seed) ; ctrl <- sample(1:ncols,groupsize)
  set.seed(seed) ; case <- sample(setdiff(1:ncols,ctrl),groupsize)
  mval_ctrl <- mval[,ctrl]
  mval_case <- mval[,case]
  # incorporate altered signals - change by +1 or -1
  mval_case[rownames(mval_case) %in% pup2,] <-  mval_case[rownames(mval_case) %in% pup2,] + delta
  mval_case[rownames(mval_case) %in% pdn2,] <-  mval_case[rownames(mval_case) %in% pdn2,] - delta
  mval2 <- cbind(mval_ctrl,mval_case)
  result <- list("mval"=mval2,"probes up"=pup2,"probes down"=pdn2,
    "genes up"=gup2,"genes down"=gdn2,
    "genesets up"=gtup,"genesets down"=gtdn)
  return(result)
}

# limma
runlimma <- function(mval,design,myann) {
  fit.reduced <- lmFit(mval,design)
  fit.reduced <- eBayes(fit.reduced)
  dm <- topTable(fit.reduced,coef=ncol(design), number = Inf)
  dm <- merge(myann,dm,by=0)
  dm <- dm[order(dm$P.Value),]
  rownames(dm) <- dm$Row.names
  dm$Row.names=NULL
  return(dm)
}


simgsa <- function(genesetdatabase, myann, mval, seed, frac_genes, frac_probes, groupsize, delta=1, num_dm_sets=50) {
  # generate gene sets
  gene_catalog <- unique(unlist(strsplit(myann$UCSC_RefGene_Name,";")))
  lengths <- unname(unlist(lapply(genesetdatabase,length)))
  gsets <- randomGeneSets(gene_catalog,lengths,seed=seed)
  # select gene sets to alter
  set.seed(seed) ; gset_mod <- sample(gsets,num_dm_sets)
  # incorporate select changes
  sim <- incorp_dm(genesets=gset_mod, myann=myann, mval=mval, seed=seed,
    frac_genes=0.5,frac_probes=0.5,groupsize=groupsize,delta=delta,
    gene_catalog=gene_catalog)
  # set up limma
  mval2 <- sim$mval
  ncols <- ncol(mval2)
  groupsize <- ncols/2
  ss <- data.frame(colnames(mval2))
  colnames(ss) <- "sample"
  ss$case <- c(rep(0,groupsize),rep(1,groupsize))
  d <- model.matrix(~ ss$case )
  dm3 <- runlimma(mval=mval2,design=d,myann=myann)
  pup3 <- rownames(subset(dm3,adj.P.Val<0.05 & logFC>0))
  pdn3 <- rownames(subset(dm3,adj.P.Val<0.05 & logFC<0))
  if ( length(pup3) < 250 ) { pup3 <- head(rownames(subset(dm3, logFC > 0)), 250) }
  if ( length(pdn3) < 250 ) { pdn3 <- head(rownames(subset(dm3, logFC < 0)), 250) }
  # convert gene sets to entrez
  suppressWarnings(suppressMessages({ gene2entrez <- mapIds(org.Hs.eg.db, gene_catalog, 'ENTREZID', 'SYMBOL') }))
  gsets_entrez <- lapply(gsets,function(gs) {
    gs2 <- unique(gene2entrez[names(gene2entrez) %in% gs])
    gs2 <- gs2[!is.na(gs2)]
    return(gs2)
  })
  suppressWarnings(suppressMessages({
    gsaup3 <- gsameth(sig.cpg=pup3, all.cpg=rownames(dm3), collection=gsets_entrez, array.type="EPIC")
    gsadn3 <- gsameth(sig.cpg=pdn3, all.cpg=rownames(dm3), collection=gsets_entrez, array.type="EPIC")
  }))
  gsig_up3 <- rownames(subset(gsaup3,FDR<0.05))
  gsig_dn3 <- rownames(subset(gsadn3,FDR<0.05))
  gtup <- names(sim[[6]])
  gtdn <- names(sim[[7]])
  UPTP=length(intersect(gsig_up3 ,gtup))
  UPFP=length(setdiff(gsig_up3 ,gtup))
  UPFN=length(setdiff(gtup,gsig_up3))
  DNTP=length(intersect(gsig_dn3 ,gtdn))
  DNFP=length(setdiff(gsig_dn3 ,gtdn))
  DNFN=length(setdiff(gtdn,gsig_dn3))
  TP=UPTP+DNTP
  FP=UPFP+DNFP
  FN=UPFN+DNFN
  TN=nrow(gsadn3)-DNTP-DNFP-DNFN-UPTP-UPFP-UPFN
  PREC=TP/(TP+FP)
  REC=TP/(TP+FN)
  F1=TP/(TP+(0.5*(FP+FN)))
  result <- c("TP"=TP,"FP"=FP,"FN"=FN,"TN"=TN,"PREC"=PREC,"REC"=REC)
  return(result)
}

# enrich parametric
ttenrich <- function(m,genesets,cores=1,testtype="selfcontained") {
  res <- mclapply( 1:length(genesets), function(i) {
    scores <- m[,1]
    gs <- genesets[i]
    name <- names(gs)
    n_members <- length(which(rownames(m) %in% gs[[1]]))
    if ( n_members > 4 ) {
      tstats <- m[which(rownames(m) %in% gs[[1]]),]
      myn <- length(tstats)
      mymean <- mean(tstats)
      mymedian <- median(tstats)
      if ( testtype == "selfcontained" ) { wt <- t.test(tstats) }
      if ( testtype == "competitive" ) { wt <- t.test(tstats,scores) }
      res <- c(name,myn,mymean,mymedian,wt$p.value)
    }
  } , mc.cores = cores)
  res_df <- do.call(rbind, res)
  rownames(res_df) <- res_df[,1]
  res_df <- res_df[,-1]
  colnames(res_df) <- c("n_genes","t_mean","t_median","pval")
  tmp <- apply(res_df,2,as.numeric)
  rownames(tmp) <- rownames(res_df)
  res_df <- tmp
  res_df <- as.data.frame(res_df)
  res_df <- res_df[order(res_df$pval),]
  res_df$logp <- -log10(res_df$pval )
  res_df$fdr <- p.adjust(res_df$pval,method="fdr")
  res_df[order(abs(res_df$pval)),]
  return(res_df)
}

# enrich non-parametric
wtenrich <- function(m,genesets,cores=1,testtype="selfcontained") {
  res <- mclapply( 1:length(genesets), function(i) {
    scores <- m[,1]
    gs <- genesets[i]
    name <- names(gs)
    n_members <- length(which(rownames(m) %in% gs[[1]]))
    if ( n_members > 4 ) {
      tstats <- m[which(rownames(m) %in% gs[[1]]),]
      myn <- length(tstats)
      mymean <- mean(tstats)
      mymedian <- median(tstats)
      if ( testtype == "selfcontained" ) { wt <- wilcox.test(tstats) }
      if ( testtype == "competitive" ) { wt <- wilcox.test(tstats,scores) }
      res <- c(name,myn,mymean,mymedian,wt$p.value)
    }
  } , mc.cores = cores)
  res_df <- do.call(rbind, res)
  rownames(res_df) <- res_df[,1]
  res_df <- res_df[,-1]
  colnames(res_df) <- c("n_genes","t_mean","t_median","pval")
  tmp <- apply(res_df,2,as.numeric)
  rownames(tmp) <- rownames(res_df)
  res_df <- tmp
  res_df <- as.data.frame(res_df)
  res_df <- res_df[order(res_df$pval),]
  res_df$logp <- -log10(res_df$pval )
  res_df$fdr <- p.adjust(res_df$pval,method="fdr")
  res_df[order(abs(res_df$pval)),]
  return(res_df)
}

# LA parametric competitive
simlac <- function(genesetdatabase, myann, mval, seed, frac_genes, frac_probes, groupsize, delta=1, num_dm_sets=50) {
  # generate gene sets
  gene_catalog <- unique(gt$gene)
  lengths <- unname(unlist(lapply(genesetdatabase,length)))
  gsets <- randomGeneSets(gene_catalog,lengths,seed=seed)
  # select gene sets to alter
  set.seed(seed) ; gset_mod <- sample(gsets,num_dm_sets)
  # incorporate select changes
  sim <- incorp_dm(genesets=gset_mod, myann=myann, mval=mval, seed=seed,
    frac_genes=0.5,frac_probes=0.5,groupsize=groupsize,delta=delta,
    gene_catalog=gene_catalog)
  # set up limma
  mval2 <- sim$mval
  ncols <- ncol(mval2)
  groupsize <- ncols/2
  ss <- data.frame(colnames(mval2))
  colnames(ss) <- "sample"
  ss$case <- c(rep(0,groupsize),rep(1,groupsize))
  d <- model.matrix(~ ss$case )
  dm3 <- runlimma(mval=mval2,design=d,myann=myann)
  dd <- merge(dm3,gt,by.x=0,by.y="probe")
  m1 <- aggregate(t ~ gene,dd,mean)
  rownames(m1) <- m1$gene
  m1$gene=NULL
  lares1 <- ttenrich(m=m1,genesets=gsets,cores=2,testtype="competitive")
  gsig_up3 <- rownames(subset(lares1, fdr < 0.05 & t_mean > 0))
  gsig_dn3 <- rownames(subset(lares1, fdr < 0.05 & t_mean < 0))
  gtup <- names(sim[[6]])
  gtdn <- names(sim[[7]])
  UPTP=length(intersect(gsig_up3 ,gtup))
  UPFP=length(setdiff(gsig_up3 ,gtup))
  UPFN=length(setdiff(gtup,gsig_up3))
  DNTP=length(intersect(gsig_dn3 ,gtdn))
  DNFP=length(setdiff(gsig_dn3 ,gtdn))
  DNFN=length(setdiff(gtdn,gsig_dn3))
  TP=UPTP+DNTP
  FP=UPFP+DNFP
  FN=UPFN+DNFN
  TN=nrow(lares1)-DNTP-DNFP-DNFN-UPTP-UPFP-UPFN
  PREC=TP/(TP+FP)
  REC=TP/(TP+FN)
  F1=TP/(TP+(0.5*(FP+FN)))
  result <- c("TP"=TP,"FP"=FP,"FN"=FN,"TN"=TN,"PREC"=PREC,"REC"=REC)
  return(result)
}

# LA parametric competitive top
simlactop <- function(genesetdatabase, myann, mval, seed, frac_genes, frac_probes, groupsize, delta=1, num_dm_sets=50) {
  # generate gene sets
  gene_catalog <- unique(gt$gene)
  lengths <- unname(unlist(lapply(genesetdatabase,length)))
  gsets <- randomGeneSets(gene_catalog,lengths,seed=seed)
  # select gene sets to alter
  set.seed(seed) ; gset_mod <- sample(gsets,num_dm_sets)
  # incorporate select changes
  sim <- incorp_dm(genesets=gset_mod, myann=myann, mval=mval, seed=seed,
    frac_genes=0.5,frac_probes=0.5,groupsize=groupsize,delta=delta,
    gene_catalog=gene_catalog)
  # set up limma
  mval2 <- sim$mval
  ncols <- ncol(mval2)
  groupsize <- ncols/2
  ss <- data.frame(colnames(mval2))
  colnames(ss) <- "sample"
  ss$case <- c(rep(0,groupsize),rep(1,groupsize))
  d <- model.matrix(~ ss$case )
  dm3 <- runlimma(mval=mval2,design=d,myann=myann)
  dd <- merge(dm3,gt,by.x=0,by.y="probe")
  m1 <- aggregate(t ~ gene,dd, function(x) {
    if (abs(max(x)) > abs(min(x))) { max(x) } else { min(x) }
  })
  rownames(m1) <- m1$gene
  m1$gene=NULL
  lares1 <- ttenrich(m=m1,genesets=gsets,cores=2,testtype="competitive")
  gsig_up3 <- rownames(subset(lares1, fdr < 0.05 & t_mean > 0))
  gsig_dn3 <- rownames(subset(lares1, fdr < 0.05 & t_mean < 0))
  gtup <- names(sim[[6]])
  gtdn <- names(sim[[7]])
  UPTP=length(intersect(gsig_up3 ,gtup))
  UPFP=length(setdiff(gsig_up3 ,gtup))
  UPFN=length(setdiff(gtup,gsig_up3))
  DNTP=length(intersect(gsig_dn3 ,gtdn))
  DNFP=length(setdiff(gsig_dn3 ,gtdn))
  DNFN=length(setdiff(gtdn,gsig_dn3))
  TP=UPTP+DNTP
  FP=UPFP+DNFP
  FN=UPFN+DNFN
  TN=nrow(lares1)-DNTP-DNFP-DNFN-UPTP-UPFP-UPFN
  PREC=TP/(TP+FP)
  REC=TP/(TP+FN)
  F1=TP/(TP+(0.5*(FP+FN)))
  result <- c("TP"=TP,"FP"=FP,"FN"=FN,"TN"=TN,"PREC"=PREC,"REC"=REC)
  return(result)
}

# LA nonparametric competitive
simnlac <- function(genesetdatabase, myann, mval, seed, frac_genes, frac_probes, groupsize, delta=1, num_dm_sets=50) {
  # generate gene sets
  gene_catalog <- unique(gt$gene)
  lengths <- unname(unlist(lapply(genesetdatabase,length)))
  gsets <- randomGeneSets(gene_catalog,lengths,seed=seed)
  # select gene sets to alter
  set.seed(seed) ; gset_mod <- sample(gsets,num_dm_sets)
  # incorporate select changes
  sim <- incorp_dm(genesets=gset_mod, myann=myann, mval=mval, seed=seed,
    frac_genes=0.5,frac_probes=0.5,groupsize=groupsize,delta=delta,
    gene_catalog=gene_catalog)
  # set up limma
  mval2 <- sim$mval
  ncols <- ncol(mval2)
  groupsize <- ncols/2
  ss <- data.frame(colnames(mval2))
  colnames(ss) <- "sample"
  ss$case <- c(rep(0,groupsize),rep(1,groupsize))
  d <- model.matrix(~ ss$case )
  dm3 <- runlimma(mval=mval2,design=d,myann=myann)
  dd <- merge(dm3,gt,by.x=0,by.y="probe")
  m1 <- aggregate(t ~ gene,dd,mean)
  rownames(m1) <- m1$gene
  m1$gene=NULL
  lares1 <- wtenrich(m=m1,genesets=gsets,cores=2,testtype="competitive")
  gsig_up3 <- rownames(subset(lares1, fdr < 0.05 & t_mean > 0))
  gsig_dn3 <- rownames(subset(lares1, fdr < 0.05 & t_mean < 0))
  gtup <- names(sim[[6]])
  gtdn <- names(sim[[7]])
  UPTP=length(intersect(gsig_up3 ,gtup))
  UPFP=length(setdiff(gsig_up3 ,gtup))
  UPFN=length(setdiff(gtup,gsig_up3))
  DNTP=length(intersect(gsig_dn3 ,gtdn))
  DNFP=length(setdiff(gsig_dn3 ,gtdn))
  DNFN=length(setdiff(gtdn,gsig_dn3))
  TP=UPTP+DNTP
  FP=UPFP+DNFP
  FN=UPFN+DNFN
  TN=nrow(lares1)-DNTP-DNFP-DNFN-UPTP-UPFP-UPFN
  PREC=TP/(TP+FP)
  REC=TP/(TP+FN)
  F1=TP/(TP+(0.5*(FP+FN)))
  result <- c("TP"=TP,"FP"=FP,"FN"=FN,"TN"=TN,"PREC"=PREC,"REC"=REC)
  return(result)
}

# LA competitive mitch
simlacm <- function(genesetdatabase, myann, mval, seed, frac_genes, frac_probes, groupsize, delta=1, num_dm_sets=50) {
  # generate gene sets
  gene_catalog <- unique(gt$gene)
  lengths <- unname(unlist(lapply(genesetdatabase,length)))
  gsets <- randomGeneSets(gene_catalog,lengths,seed=seed)
  # select gene sets to alter
  set.seed(seed) ; gset_mod <- sample(gsets,num_dm_sets)
  # incorporate select changes
  sim <- incorp_dm(genesets=gset_mod, myann=myann, mval=mval, seed=seed,
    frac_genes=0.5,frac_probes=0.5,groupsize=groupsize,delta=delta,
    gene_catalog=gene_catalog)
  # set up limma
  mval2 <- sim$mval
  ncols <- ncol(mval2)
  groupsize <- ncols/2
  ss <- data.frame(colnames(mval2))
  colnames(ss) <- "sample"
  ss$case <- c(rep(0,groupsize),rep(1,groupsize))
  d <- model.matrix(~ ss$case )
  dm3 <- runlimma(mval=mval2,design=d,myann=myann)
  dd <- merge(dm3,gt,by.x=0,by.y="probe")
  m1 <- aggregate(t ~ gene,dd,mean)
  rownames(m1) <- m1$gene
  m1$gene=NULL
  lamres1 <- runmitch(m=m1,genesets=gsets,cores=2)
  gsig_up3 <- rownames( subset( lamres1, p.adjustANOVA < 0.05 & s.dist > 0 ) )
  gsig_dn3 <- rownames( subset( lamres1, p.adjustANOVA < 0.05 & s.dist < 0 ) )
  gtup <- names(sim[[6]])
  gtdn <- names(sim[[7]])
  UPTP=length(intersect(gsig_up3 ,gtup))
  UPFP=length(setdiff(gsig_up3 ,gtup))
  UPFN=length(setdiff(gtup,gsig_up3))
  DNTP=length(intersect(gsig_dn3 ,gtdn))
  DNFP=length(setdiff(gsig_dn3 ,gtdn))
  DNFN=length(setdiff(gtdn,gsig_dn3))
  TP=UPTP+DNTP
  FP=UPFP+DNFP
  FN=UPFN+DNFN
  TN=nrow(lamres1)-DNTP-DNFP-DNFN-UPTP-UPFP-UPFN
  PREC=TP/(TP+FP)
  REC=TP/(TP+FN)
  F1=TP/(TP+(0.5*(FP+FN)))
  result <- c("TP"=TP,"FP"=FP,"FN"=FN,"TN"=TN,"PREC"=PREC,"REC"=REC)
  return(result)
}

# LA rank competitive mitch
simlrm <- function(genesetdatabase, myann, mval, seed, frac_genes, frac_probes, groupsize, delta=1, num_dm_sets=50) {
  # generate gene sets
  gene_catalog <- unique(gt$gene)
  lengths <- unname(unlist(lapply(genesetdatabase,length)))
  gsets <- randomGeneSets(gene_catalog,lengths,seed=seed)
  # select gene sets to alter
  set.seed(seed) ; gset_mod <- sample(gsets,num_dm_sets)
  # incorporate select changes
  sim <- incorp_dm(genesets=gset_mod, myann=myann, mval=mval, seed=seed,
    frac_genes=0.5,frac_probes=0.5,groupsize=groupsize,delta=delta,
    gene_catalog=gene_catalog)
  # set up limma
  mval2 <- sim$mval
  ncols <- ncol(mval2)
  groupsize <- ncols/2
  ss <- data.frame(colnames(mval2))
  colnames(ss) <- "sample"
  ss$case <- c(rep(0,groupsize),rep(1,groupsize))
  d <- model.matrix(~ ss$case )
  dm3 <- runlimma(mval=mval2,design=d,myann=myann)
  # rank probes first
  dm3$rank <-  rank(dm3$t) - nrow(subset(dm3,t<0))
  mm <- merge(dm3,gt,by.x=0,by.y="probe")
  head(mm)
  mma <-aggregate(mm$rank ~ gene,mm,mean)
  rownames(mma) <- mma$gene
  mma$gene = NULL
  colnames(mma) <- "meanrank"
  lrmres <- runmitch(m=mma,genesets=gsets,cores=2)
  gsig_up3 <- rownames( subset( lrmres, p.adjustANOVA < 0.05 & s.dist > 0 ) )
  gsig_dn3 <- rownames( subset( lrmres, p.adjustANOVA < 0.05 & s.dist < 0 ) )
  gtup <- names(sim[[6]])
  gtdn <- names(sim[[7]])
  UPTP=length(intersect(gsig_up3 ,gtup))
  UPFP=length(setdiff(gsig_up3 ,gtup))
  UPFN=length(setdiff(gtup,gsig_up3))
  DNTP=length(intersect(gsig_dn3 ,gtdn))
  DNFP=length(setdiff(gsig_dn3 ,gtdn))
  DNFN=length(setdiff(gtdn,gsig_dn3))
  TP=UPTP+DNTP
  FP=UPFP+DNFP
  FN=UPFN+DNFN
  TN=nrow(lrmres)-DNTP-DNFP-DNFN-UPTP-UPFP-UPFN
  PREC=TP/(TP+FP)
  REC=TP/(TP+FN)
  F1=TP/(TP+(0.5*(FP+FN)))
  result <- c("TP"=TP,"FP"=FP,"FN"=FN,"TN"=TN,"PREC"=PREC,"REC"=REC)
  return(result)
}

# AL approach parametric competitive
simalc <- function(genesetdatabase, myann, mval, seed, frac_genes, frac_probes, groupsize, delta=1, num_dm_sets=50) {
  # generate gene sets
  gene_catalog <- unique(gt$gene)
  lengths <- unname(unlist(lapply(genesetdatabase,length)))
  gsets <- randomGeneSets(gene_catalog,lengths,seed=seed)
  # select gene sets to alter
  set.seed(seed) ; gset_mod <- sample(gsets,num_dm_sets)
  # incorporate select changes
  sim <- incorp_dm(genesets=gset_mod, myann=myann, mval=mval, seed=seed,
    frac_genes=0.5,frac_probes=0.5,groupsize=groupsize,delta=delta,
    gene_catalog=gene_catalog)
  # set up limma
  mval2 <- sim$mval
  ncols <- ncol(mval2)
  groupsize <- ncols/2
  ss <- data.frame(colnames(mval2))
  colnames(ss) <- "sample"
  ss$case <- c(rep(0,groupsize),rep(1,groupsize))
  d <- model.matrix(~ ss$case )
  # al pipeline
  mm <- merge(mval2,gt,by.x=0,by.y="probe")
  mm$Row.names = NULL
  a <- aggregate(. ~ gene,mm,mean)
  rownames(a) <- a$gene
  a$gene=NULL
  fit.reduced <- lmFit(a,d)
  fit.reduced <- eBayes(fit.reduced)
  al <- topTable(fit.reduced,coef=ncol(d), number = Inf)
  m1 <- as.data.frame(al$t)
  rownames(m1) <- rownames(al)
  colnames(m1) <- "t"
  alres1 <- ttenrich(m=m1,genesets=gsets,cores=2,testtype="competitive")
  # summarise results
  gsig_up3 <- rownames(subset(alres1, fdr < 0.05 & t_mean > 0))
  gsig_dn3 <- rownames(subset(alres1, fdr < 0.05 & t_mean < 0))
  gtup <- names(sim[[6]])
  gtdn <- names(sim[[7]])
  UPTP=length(intersect(gsig_up3 ,gtup))
  UPFP=length(setdiff(gsig_up3 ,gtup))
  UPFN=length(setdiff(gtup,gsig_up3))
  DNTP=length(intersect(gsig_dn3 ,gtdn))
  DNFP=length(setdiff(gsig_dn3 ,gtdn))
  DNFN=length(setdiff(gtdn,gsig_dn3))
  TP=UPTP+DNTP
  FP=UPFP+DNFP
  FN=UPFN+DNFN
  TN=nrow(alres1)-DNTP-DNFP-DNFN-UPTP-UPFP-UPFN
  PREC=TP/(TP+FP)
  REC=TP/(TP+FN)
  F1=TP/(TP+(0.5*(FP+FN)))
  result <- c("TP"=TP,"FP"=FP,"FN"=FN,"TN"=TN,"PREC"=PREC,"REC"=REC)
  return(result)
}

# AL approach nonparametric competitive
simnalc <- function(genesetdatabase, myann, mval, seed, frac_genes, frac_probes, groupsize, delta=1, num_dm_sets=50) {
  # generate gene sets
  gene_catalog <- unique(gt$gene)
  lengths <- unname(unlist(lapply(genesetdatabase,length)))
  gsets <- randomGeneSets(gene_catalog,lengths,seed=seed)
  # select gene sets to alter
  set.seed(seed) ; gset_mod <- sample(gsets,num_dm_sets)
  # incorporate select changes
  sim <- incorp_dm(genesets=gset_mod, myann=myann, mval=mval, seed=seed,
    frac_genes=0.5,frac_probes=0.5,groupsize=groupsize,delta=delta,
    gene_catalog=gene_catalog)
  # set up limma
  mval2 <- sim$mval
  ncols <- ncol(mval2)
  groupsize <- ncols/2
  ss <- data.frame(colnames(mval2))
  colnames(ss) <- "sample"
  ss$case <- c(rep(0,groupsize),rep(1,groupsize))
  d <- model.matrix(~ ss$case )
  # al pipeline
  mm <- merge(mval2,gt,by.x=0,by.y="probe")
  mm$Row.names = NULL
  a <- aggregate(. ~ gene,mm,mean)
  rownames(a) <- a$gene
  a$gene=NULL
  fit.reduced <- lmFit(a,d)
  fit.reduced <- eBayes(fit.reduced)
  al <- topTable(fit.reduced,coef=ncol(d), number = Inf)
  m1 <- as.data.frame(al$t)
  rownames(m1) <- rownames(al)
  colnames(m1) <- "t"
  alres1 <- wtenrich(m=m1,genesets=gsets,cores=2,testtype="competitive")
  # summarise results
  gsig_up3 <- rownames(subset(alres1, fdr < 0.05 & t_mean > 0))
  gsig_dn3 <- rownames(subset(alres1, fdr < 0.05 & t_mean < 0))
  gtup <- names(sim[[6]])
  gtdn <- names(sim[[7]])
  UPTP=length(intersect(gsig_up3 ,gtup))
  UPFP=length(setdiff(gsig_up3 ,gtup))
  UPFN=length(setdiff(gtup,gsig_up3))
  DNTP=length(intersect(gsig_dn3 ,gtdn))
  DNFP=length(setdiff(gsig_dn3 ,gtdn))
  DNFN=length(setdiff(gtdn,gsig_dn3))
  TP=UPTP+DNTP
  FP=UPFP+DNFP
  FN=UPFN+DNFN
  TN=nrow(alres1)-DNTP-DNFP-DNFN-UPTP-UPFP-UPFN
  PREC=TP/(TP+FP)
  REC=TP/(TP+FN)
  F1=TP/(TP+(0.5*(FP+FN)))
  result <- c("TP"=TP,"FP"=FP,"FN"=FN,"TN"=TN,"PREC"=PREC,"REC"=REC)
  return(result)
}

runmitch <- function(m,genesets,cores=1) {
  suppressMessages({ mres <- mitch_calc(m,genesets,minsetsize=5,cores=cores) })
  mres <- mres$enrichment_result
  rownames(mres) <- mres$set
  mres$set=NULL
  return(mres)
}

simalm <- function(genesetdatabase, myann, mval, seed, frac_genes, frac_probes, groupsize, delta=1, num_dm_sets=50) {
  # generate gene sets
  gene_catalog <- unique(gt$gene)
  lengths <- unname(unlist(lapply(genesetdatabase,length)))
  gsets <- randomGeneSets(gene_catalog,lengths,seed=seed)
  # select gene sets to alter
  set.seed(seed) ; gset_mod <- sample(gsets,num_dm_sets)
  # incorporate select changes
  sim <- incorp_dm(genesets=gset_mod, myann=myann, mval=mval, seed=seed,
    frac_genes=0.5,frac_probes=0.5,groupsize=groupsize,delta=delta,
    gene_catalog=gene_catalog)
  # set up limma
  mval2 <- sim$mval
  ncols <- ncol(mval2)
  groupsize <- ncols/2
  ss <- data.frame(colnames(mval2))
  colnames(ss) <- "sample"
  ss$case <- c(rep(0,groupsize),rep(1,groupsize))
  d <- model.matrix(~ ss$case )
  # alm pipeline
  mm <- merge(mval2,gt,by.x=0,by.y="probe")
  mm$Row.names = NULL
  a <- aggregate(. ~ gene,mm,mean)
  rownames(a) <- a$gene
  a$gene=NULL
  fit.reduced <- lmFit(a,d)
  fit.reduced <- eBayes(fit.reduced)
  al <- topTable(fit.reduced,coef=ncol(d), number = Inf)
  m1 <- as.data.frame(al$t)
  rownames(m1) <- rownames(al)
  colnames(m1) <- "t"
  almres1 <- runmitch(m=m1,genesets=gsets,cores=2)
  # summarise results
  gsig_up3 <- rownames( subset( almres1, p.adjustANOVA < 0.05 & s.dist > 0 ) )
  gsig_dn3 <- rownames( subset( almres1, p.adjustANOVA < 0.05 & s.dist < 0 ) )
  gtup <- names(sim[[6]])
  gtdn <- names(sim[[7]])
  UPTP=length(intersect(gsig_up3 ,gtup))
  UPFP=length(setdiff(gsig_up3 ,gtup))
  UPFN=length(setdiff(gtup,gsig_up3))
  DNTP=length(intersect(gsig_dn3 ,gtdn))
  DNFP=length(setdiff(gsig_dn3 ,gtdn))
  DNFN=length(setdiff(gtdn,gsig_dn3))
  TP=UPTP+DNTP
  FP=UPFP+DNFP
  FN=UPFN+DNFN
  TN=nrow(almres1)-DNTP-DNFP-DNFN-UPTP-UPFP-UPFN
  PREC=TP/(TP+FP)
  REC=TP/(TP+FN)
  F1=TP/(TP+(0.5*(FP+FN)))
  result <- c("TP"=TP,"FP"=FP,"FN"=FN,"TN"=TN,"PREC"=PREC,"REC"=REC)
  return(result)
}

gsagg <- function(x,genesets,cores=1) {
  meds <- mclapply(1:length(genesets), function(i) {
    gs = genesets[[i]]
    xx <- x[rownames(x) %in% gs,]
    med <- apply(xx,2,mean)
  },mc.cores=cores)
  mymed <- do.call(rbind,meds)
  rownames(mymed) <- names(genesets)
  as.data.frame(mymed)
}

aalimma <- function(agag,design) {
  fit.reduced <- lmFit(agag,design)
  fit.reduced <- eBayes(fit.reduced)
  dmagg <- topTable(fit.reduced,coef=ncol(design), number = Inf)
  return(dmagg)
}

aal <- function(mval,myann,genesets,design,cores=1) {
  medf <- chragg(mval,myann,cores=cores)
  agag <- gsagg(x=medf,genesets=genesets,cores=cores)
  aalres <- aalimma(agag=agag,design=design)
  return(aalres)
}

simaa <- function(genesetdatabase, myann, mval, seed, frac_genes, frac_probes, groupsize, delta=1, num_dm_sets=50) {
  # generate gene sets
  gene_catalog <- unique(gt$gene)
  lengths <- unname(unlist(lapply(genesetdatabase,length)))
  gsets <- randomGeneSets(gene_catalog,lengths,seed=seed)
  # select gene sets to alter
  set.seed(seed) ; gset_mod <- sample(gsets,num_dm_sets)
  # incorporate select changes
  sim <- incorp_dm(genesets=gset_mod, myann=myann, mval=mval, seed=seed,
    frac_genes=0.5,frac_probes=0.5,groupsize=groupsize,delta=delta,
    gene_catalog=gene_catalog)
  # set up limma
  mval2 <- sim$mval
  ncols <- ncol(mval2)
  groupsize <- ncols/2
  ss <- data.frame(colnames(mval2))
  colnames(ss) <- "sample"
  ss$case <- c(rep(0,groupsize),rep(1,groupsize))
  d <- model.matrix(~ ss$case )
  mm <- merge(mval2,gt,by.x=0,by.y="probe")
  mm$Row.names = NULL
  a <- aggregate(. ~ gene,mm,mean)
  rownames(a) <- a$gene
  a$gene=NULL
  mystack <- stack(gsets)
  mmm <- merge(a,mystack,by.x=0,by.y="values")
  mmm$Row.names=NULL
  aa <- aggregate(. ~ ind,mmm,mean)
  rownames(aa) <- aa$ind
  aa$ind=NULL
  fit.reduced <- lmFit(aa,d)
  fit.reduced <- eBayes(fit.reduced)
  aares1 <- topTable(fit.reduced,coef=ncol(d), number = Inf)
  # summarise results
  gsig_up3 <- rownames(subset(aares1, adj.P.Val < 0.05 & logFC > 0))
  gsig_dn3 <- rownames(subset(aares1, adj.P.Val < 0.05 & logFC < 0))
  gtup <- names(sim[[6]])
  gtdn <- names(sim[[7]])
  UPTP=length(intersect(gsig_up3 ,gtup))
  UPFP=length(setdiff(gsig_up3 ,gtup))
  UPFN=length(setdiff(gtup,gsig_up3))
  DNTP=length(intersect(gsig_dn3 ,gtdn))
  DNFP=length(setdiff(gsig_dn3 ,gtdn))
  DNFN=length(setdiff(gtdn,gsig_dn3))
  TP=UPTP+DNTP
  FP=UPFP+DNFP
  FN=UPFN+DNFN
  TN=nrow(aares1)-DNTP-DNFP-DNFN-UPTP-UPFP-UPFN
  PREC=TP/(TP+FP)
  REC=TP/(TP+FN)
  F1=TP/(TP+(0.5*(FP+FN)))
  result <- c("TP"=TP,"FP"=FP,"FN"=FN,"TN"=TN,"PREC"=PREC,"REC"=REC)
  return(result)
}

F1 <- function(x,y) {
  ( 2 * x * y ) / ( x + y )
}
