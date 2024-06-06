library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
myann <- data.frame(anno[,c("UCSC_RefGene_Name","UCSC_RefGene_Group","Islands_Name","Relation_to_Island")])
gp <- myann[,"UCSC_RefGene_Name",drop=FALSE]

x <- table(unlist(lapply(strsplit(gp$UCSC_RefGene_Name,";"), function(x) { length(unique(x)) })))
y <- x[1:11]
y <- c(y,sum(x[12:length(x)]))
names(y)[length(y)] = ">10"

pdf("r2bars.pdf")
barplot(y,ylim=c(0,6e5), xlab="no. of annotated genes", ylab="no. probes" )
 text((1:12*1.2)-0.54,y+20000,labels=y)
dev.off()

ul <- unlist(lapply(strsplit(gp$UCSC_RefGene_Name,";"), function(x) { length(unique(x)) }))
length(which(ul>1))
#1gene=559007;>1gene=57591

ut <- table(table(unlist(strsplit( gp[which(ul>1),1] , ";"))))
ut <- c(ut[1:20],sum(ut[21:length(ut)]))
names(ut)[length(ut)] <- ">20"

gp$shared <- lapply(strsplit(gp$UCSC_RefGene_Name,";"), function(x) { length(unique(x)) }) > 1

genes <- unique(unlist(strsplit(gp$UCSC_RefGene_Name,";")))
genelist <- strsplit(gp$UCSC_RefGene_Name,";")
genelist <- lapply(genelist,unique)
names(genelist) <- rownames(gp)

gs <-stack(genelist)

m <- merge(gs,gp,by.x="ind",by.y=0,all.x)

tab <- mclapply(genes,function(g) {
  m[m$values==g,"shared"]
} ,mc.cores=32)

names(tab) <- genes

uni <- lapply(tab, function(x) { SHARED=length(which(x==TRUE)) ; UNIQUE=length(which(x==FALSE)) ; c(UNIQUE,SHARED) } )

uni <- as.data.frame(do.call(rbind,uni))
colnames(uni) <- c("unique","shared")
uni$total <- rowSums(uni)

uni$unique_proportion <- uni$unique / uni$total
hist(uni$unique_proportion,xlab="proportion of unique probes",ylab="no. genes",main="probe uniqueness")

library(mitch)
sets <- mitch::gmt_import("c2.cp.reactome.v2023.1.Hs.symbols.gmt")
setgenes <- unique(unlist(sets))
length(setgenes)

length(which(uni$unique_proportion>0.8))

length(intersect(rownames(uni[which(uni$unique_proportion>0.8),]), setgenes))

8262/18300
#45%

length(which(uni$unique_proportion<0.2))
length(intersect(rownames(uni[which(uni$unique_proportion<0.2),]), setgenes))
647/4811
#13%
