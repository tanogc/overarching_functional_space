#################################################################################################################
#                                                                                                               #
# Supplementary Data 1                                                                                          #
#                                                                                                               #
# R script to create regional and continental functional spaces (FS),                                           # 
# and assess latitudinal patterns of FD metrics                                                                 # 
#                                                                                                               #
#                                                                                                               #
# A trait space at an overarching-scale yields more conclusive macroecological patterns of functional diversity #
#                                                                                                               #
#                                                                                                               #
# Cesc Múrria, Gabone Iturrarte and Cayetano Gutiérrez-Cánovas                                                  #
#                                                                                                               #
# Code written by Cayetano Gutiérrez-Cánovas and Gabone Iturralde                                               #
# email for queries: cayeguti@um.es                                                                             #
#                                                                                                               #
#################################################################################################################

# working directory
setwd("your_folder")
res_plot<-paste(getwd(),"/res/",sep="")
plot_assum<-paste(getwd(),"/res/assum/",sep="")

# Loading packages
library(ade4)
library(vegan)
library(plyr)
library(FD)
library(MuMIn)
library(adegraphics)

# Additional functions
source("0_FD_functions.R")
source("0_quality_funct_space_fromdist.R")


### Loading matrices

# Taxonomic data
site.mat <- read.table("site_taxa.txt", header = TRUE, sep = "\t") # Loading site x taxa matrix
site.mat1 <- site.mat[-9,-1]
rownames(site.mat1) <- site.mat$Genus[-9]
reg.mat<-read.table("regional_taxa.txt", header = TRUE, sep = "\t",row.names=1)

# Trait data
traits_raw <- read.table("traits.txt", header = TRUE, sep = "\t")
traits <- subset(traits_raw, select = SIZE1:PERATT)
rownames(traits) <- traits_raw$Genus

# Trait category blocks
traits <- prep.fuzzy.df(traits, c(7,2,3,4,8,5,4,4,8))

# Checking the genus in traits matrix are the same of that in community matrix
traits.tax.site <- traits[(intersect(rownames(traits), colnames(site.mat1))),]
tax.mat.site2 <- subset(site.mat1, select = rownames(traits.tax.site))

##############################################################################
# Building the continental-scale FS
##############################################################################

# Calculate species x species (Pavoine's) Gower dissimilarity matrix based on effect traits
trait_df <- ktab.list.df(list(traits[which(rowSums(traits)==9),]))
tr.dist <- dist.ktab(trait_df, type= c("F"))

# This data set has missing values
is.euclid(tr.dist)

# Estimating the optimum number of dimensions
ov_qual_fs <- quality_funct_space_fromdist(tr.dist, nbdim=15) 
ov_qual_fs$meanSD # results suggest 10D

# Continental-scale FS (PCoA)
ov.pco <- dudi.pco(tr.dist, scan = F, nf = 10)

# Explained variance by each axis
round(ov.pco$eig[1:10]/sum(ov.pco$eig),2)

# Cumulative fariance explained by all FS axes (10 axes) #7D, 57.9%
cumsum(ov.pco$eig)[7]/sum(ov.pco$eig) 

# Spearman rank correlations between original trait categories 
# and functional space axes
cor.res<-round(cor(traits[which(rowSums(traits)==9),],ov.pco$li, method="spearman"),2)
cor.res
write.table(cor.res, "cor_res.txt", sep="\t")

# Selecting taxa ocurring in the empirical datasets (site.mat1) & (reg.mat)
funspace.site <- ov.pco$li[(intersect(rownames(ov.pco$li), colnames(site.mat1))),]
site.mat1<-site.mat1[,(intersect(rownames(ov.pco$li), colnames(site.mat1)))]
reg.mat<-reg.mat[,(intersect(rownames(ov.pco$li), colnames(reg.mat)))]

###########################################################################
# Coninental scale (study-specific) FS (PCoA)

## Contiental-based FS

# Matrix to store correlations between continental and regional FS axes
cor.b<-data.frame(matrix(NA, 1, 10))
rownames(cor.b)<-c("Continental-biased")
colnames(cor.b)<-paste("Axis", 1:10)

# Selecting taxa from the regional pool
traits.con<-traits[(intersect(rownames(traits), colnames(reg.mat))),]

# Gower dissimilarity
trait_con_df <- ktab.list.df(list(traits.con[which(rowSums(traits.con)==9),]))
tr.dist_con <- dist.ktab(trait_con_df, type= c("F"))

# Estimating the optimum number of dimensions
con_qual_fs <- quality_funct_space_fromdist(tr.dist_con, nbdim=15) 
con_qual_fs$meanSD # 0.006 at 7D

# FS (PCoA)
con.pco <- dudi.pco(tr.dist_con, scan = F, nf = 10)

cumsum(con.pco$eig)[7]/sum(con.pco$eig) # Variance explained by FS (PCoA); 67.39

# Correlation between regional and continental-scale FS axes
con.cfs<-ov.pco$li[(intersect(rownames(ov.pco$li), rownames(con.pco$li))),]

cor.b[1,]<-round(diag(cor(con.pco$li , con.cfs, method="pearson")),3)

###########################################################################
# Building regional FSs (PCoA)

# Regions
rownames(reg.mat)

# Matrix to store correlations between continental and regional FS axes
cor.reg<-data.frame(matrix(NA, 6, 10))
rownames(cor.reg)<-c("Rif", "Beticas", "Picos", "Jura", "Carpathians", "Sweeden")
colnames(cor.reg)<-paste("Axis", 1:10)

rownames(reg.mat)

###########################################################################
### Rif ###

# Selecting taxa from the regional pool
rif.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[1,]>0)]))),]
rif.gr<-data.frame(group2=traits_raw$Group2,traits)[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[1,]>0)]))),]

# Gower dissimilarity
traits_rif <- ktab.list.df(list(rif.tr[which(rowSums(rif.tr)==9),]))
tr.dist_rif <- dist.ktab(traits_rif, type= c("F"))

# FS (PCoA)
rif.pco <- dudi.pco(tr.dist_rif, scan = F, nf = 10)
rif.tax <- site.mat1[1:8,(intersect(rownames(rif.pco$li), colnames(site.mat1)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_rif, nbdim=15) 
qual_fs$meanSD # 0.006 at 7D

cumsum(rif.pco$eig)[7]/sum(rif.pco$eig) # Variance explained by FS (PCoA); 68.9

# Correlation betwewen regional and continental-scale FS axes
rif.cfs<-ov.pco$li[(intersect(rownames(ov.pco$li), rownames(rif.pco$li))),]

cor.reg[1,]<-round(diag(cor(rif.pco$li , rif.cfs, method="pearson")),3)

###########################################################################
### Betic ###

# Selecting taxa from the regional pool
bet.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[2,]>0)]))),]

# Gower dissimilarity
traits_bet <- ktab.list.df(list(bet.tr[which(rowSums(bet.tr)==9),]))
tr.dist_bet <- dist.ktab(traits_bet, type= c("F"))

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_bet, nbdim=15) 
qual_fs$meanSD # 0.008 at 7D

# FS (PCoA)
bet.pco <- dudi.pco(tr.dist_bet, scan = F, nf = 10)
bet.tax <- site.mat1[9:18,(intersect(rownames(bet.pco$li), colnames(site.mat1)))]

cumsum(bet.pco$eig)[7]/sum(bet.pco$eig) # Variance explained by FS (PCoA); 69.1

# Correlation betwewen regional and continental-scale FS axes
bet.cfs<-ov.pco$li[(intersect(rownames(ov.pco$li), rownames(bet.pco$li))),]

cor.reg[2,]<-round(diag(cor(bet.pco$li , bet.cfs, method="pearson")),3)

###########################################################################
### Picos ###

# Selecting taxa from the regional pool
pic.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[3,]>0)]))),]

# Gower dissimilarity
traits_pic <- ktab.list.df(list(pic.tr[which(rowSums(pic.tr)==9),]))
tr.dist_pic <- dist.ktab(traits_pic, type= c("F"))

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_pic, nbdim=15) 
qual_fs$meanSD # 0.006 at 7D

# FS (PCoA)
pic.pco <- dudi.pco(tr.dist_pic, scan = F, nf = 10)
pic.tax <- site.mat1[19:28,(intersect(rownames(pic.pco$li), colnames(site.mat1)))]

cumsum(pic.pco$eig)[7]/sum(pic.pco$eig) # Variance explained by FS (PCoA); 69.5

# Correlation betwewen regional and continental-scale FS axes
pic.cfs<-ov.pco$li[(intersect(rownames(ov.pco$li), rownames(pic.pco$li))),]

cor.reg[3,]<-round(diag(cor(pic.pco$li , pic.cfs, method="pearson")),3)

###########################################################################
### Jura ###

# Selecting taxa from the regional pool
jura.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[4,]>0)]))),]

# Gower dissimilarity
traits_jura <- ktab.list.df(list(jura.tr[which(rowSums(jura.tr)==9),]))
tr.dist_jura <- dist.ktab(traits_jura, type= c("F"))

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_jura, nbdim=15) 
qual_fs$meanSD # 0.009 at 6D

# FS (PCoA)
jura.pco <- dudi.pco(tr.dist_jura, scan = F, nf = 10)
jura.tax <- site.mat1[29:39,(intersect(rownames(jura.pco$li), colnames(site.mat1)))]

cumsum(jura.pco$eig)[7]/sum(jura.pco$eig) # Variance explained by FS (PCoA); 69.2

# Correlation betwewen regional and continental-scale FS axes
jura.cfs<-ov.pco$li[(intersect(rownames(ov.pco$li), rownames(jura.pco$li))),]

cor.reg[4,]<-round(diag(cor(jura.pco$li , jura.cfs, method="pearson")),3)

###########################################################################
### Carpathian ###

# Selecting taxa from the regional pool
car.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[5,]>0)]))),]

# Gower dissimilarity
traits_car <- ktab.list.df(list(car.tr[which(rowSums(car.tr)==9),]))
tr.dist_car <- dist.ktab(traits_car, type= c("F"))

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_car, nbdim=15) 
qual_fs$meanSD # 0.005 at 6D

# FS (PCoA)
car.pco <- dudi.pco(tr.dist_car, scan = F, nf = 10)
car.tax <- site.mat1[40:50,(intersect(rownames(car.pco$li), colnames(site.mat1)))]

cumsum(car.pco$eig)[7]/sum(car.pco$eig) # Variance explained by FS (PCoA); 70.5

# Correlation betwewen regional and continental-scale FS axes
car.cfs<-ov.pco$li[(intersect(rownames(ov.pco$li), rownames(car.pco$li))),]

cor.reg[5,]<-round(diag(cor(car.pco$li , car.cfs, method="pearson")),3)

###########################################################################
### Sweden ###

# Selecting taxa from the regional pool
swe.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[6,]>0)]))),]

# Gower dissimilarity
traits_swe <- ktab.list.df(list(swe.tr[which(rowSums(swe.tr)==9),]))
tr.dist_swe <- dist.ktab(traits_swe, type= c("F"))

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_swe, nbdim=15) 
qual_fs$meanSD # 0.008 at 6D

# FS (PCoA)
swe.pco <- dudi.pco(tr.dist_swe, scan = F, nf = 10)
swe.tax <- site.mat1[51:61,(intersect(rownames(swe.pco$li), colnames(site.mat1)))]

cumsum(swe.pco$eig)[7]/sum(swe.pco$eig) # Variance explained by FS (PCoA); 70.9

# Correlation betwewen regional and continental-scale FS axes
swe.cfs<-ov.pco$li[(intersect(rownames(ov.pco$li), rownames(swe.pco$li))),]

cor.reg[6,]<-round(diag(cor(swe.pco$li , swe.cfs, method="pearson")),3)

###########################################################################
# Plot showing regional taxon pool representation within the continental-scale FSs
###########################################################################

pdf(file="func_space_no_points.pdf",onefile=T,width=12,height=12)  
par(mfrow=c(3,3),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))

c("gold", "darkorange","red", "deeppink", "steelblue1", "royalblue4")->col.reg

# Loop to plot the 6 FSs

for (i in 1:6)
  
{
  plot(range(ov.pco$li[1]), range(ov.pco$li[2]),type = "n", main=row.names(reg.mat)[i],xlab="Functional axis 1", ylab="Functional axis 2")
  points(ov.pco$li[,c(1,2)],col="#4D4D4D",pch=".",cex=7)
  plot_chull2D(ov.pco$li[,c(1,2)],col_ch="#1E90FF50",border_ch="#4D4D4D")
  plot_chull2D(ov.pco$li[(intersect(rownames(ov.pco$li), colnames(reg.mat[which(reg.mat[i,]>0)]))),c(1,2)],col_ch=adjustcolor(col.reg[i], alpha.f = 0.4),border_ch=adjustcolor(col.reg[i], alpha.f = 0.4))
  colMeans(ov.pco$li[,c(1,2)])->cent_ov
  points(cent_ov[1],cent_ov[2],col="black",pch="+",cex=4)
  colMeans(ov.pco$li[(intersect(rownames(ov.pco$li), colnames(reg.mat[which(reg.mat[i,]>0)]))),c(1,2)])->cent_r
  points(cent_r[1],cent_r[2],col="red",pch="+",cex=2.5)
  
 }

plot(range(ov.pco$li[1]), range(ov.pco$li[2]),type = "n", main="Continental",xlab="Functional axis 1", ylab="Functional axis 2")
points(ov.pco$li[,c(1,2)],col="#4D4D4D",pch=".",cex=7)
plot_chull2D(ov.pco$li[,c(1,2)],col_ch="#1E90FF50",border_ch="#4D4D4D")
plot_chull2D(ov.pco$li[(intersect(rownames(ov.pco$li), colnames(reg.mat[which(colSums(reg.mat)>0)]))),c(1,2)],col_ch=adjustcolor("blue", alpha.f = 0.4),border_ch=adjustcolor("blue", alpha.f = 0.4))
colMeans(ov.pco$li[,c(1,2)])->cent_ov
points(cent_ov[1],cent_ov[2],col="black",pch="+",cex=4)
colMeans(ov.pco$li[(intersect(rownames(ov.pco$li), colnames(reg.mat))),c(1,2)])->cent_r
points(cent_r[1],cent_r[2],col="red",pch="+",cex=2.5)

dev.off()

### Barplot showing the correlations between continental and regional FS axes

# Mean absolute correlation between regional and continental-scale FS axes
round(colMeans(abs(rbind(cor.b, cor.reg))),2)

range(round(colMeans(abs(rbind(cor.b, cor.reg))),2)[3:7])

pdf(file="FS_cor.pdf",onefile=T,width=6,height=6)

par(cex.lab=1.25, cex.axis=1.25, mfrow=c(1,1))
colnames(cor.reg)<-1:10
barplot(colMeans(abs(cor.reg))[1:7], col="black", ylab = "Mean absolute correlation", xlab="Axis", ylim=c(0,1))

dev.off()

### Plotting regional-scale FS

con.gr<-data.frame(group2=traits_raw$Group2,traits)[(intersect(rownames(traits), colnames(reg.mat))),]
rif.gr<-data.frame(group2=traits_raw$Group2,traits)[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[1,]>0)]))),]
bet.gr<-data.frame(group2=traits_raw$Group2,traits)[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[2,]>0)]))),]
pic.gr<-data.frame(group2=traits_raw$Group2,traits)[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[3,]>0)]))),]
jura.gr<-data.frame(group2=traits_raw$Group2,traits)[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[4,]>0)]))),]
car.gr<-data.frame(group2=traits_raw$Group2,traits)[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[5,]>0)]))),]
swe.gr<-data.frame(group2=traits_raw$Group2,traits)[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[6,]>0)]))),]

pdf(file="continental_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(con.pco$li, fac = con.gr$group2[which(rowSums(traits.con)==9)], plines.col = 1:18, col = TRUE)
dev.off()

pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(rif.pco$li, fac = rif.gr$group2[which(rowSums(rif.tr)==9)], plines.col = 1:18, col = TRUE)
dev.off()

pdf(file="bet_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(bet.pco$li, fac = bet.gr$group2[which(rowSums(bet.tr)==9)], plines.col = 1:18, col = TRUE)
dev.off()

pdf(file="pic_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(pic.pco$li, fac = pic.gr$group2[which(rowSums(pic.tr)==9)], plines.col = 1:18, col = TRUE)
dev.off()

pdf(file="jura_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(jura.pco$li, fac = jura.gr$group2[which(rowSums(jura.tr)==9)], plines.col = 1:18, col = TRUE)
dev.off()

pdf(file="car_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(car.pco$li, fac = car.gr$group2[which(rowSums(car.tr)==9)], plines.col = 1:18, col = TRUE)
dev.off()

pdf(file="swe_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(swe.pco$li, fac = swe.gr$group2[which(rowSums(swe.tr)==9)], plines.col = 1:18, col = TRUE)
dev.off()

#######################################################################################

# Computing FD metrics (FRic, FDis, FEve)

# Using Overarching-scale FS
chull_3d(ov.pco$li,m=6,prec="Qt") -> ch_ov
fric_3d(site.mat1,ov.pco$li,m=6,prec="QJ", fric.3d.max = ch_ov) -> FRic_ov 
fdisp_k_sub(tr.dist, site.mat1,tax_sub=colnames(site.mat1), m=10)$FDis -> FDis_ov
feve_k(ov.pco$li,site.mat1,m=10) -> FEve_ov 

# Using Contiental-scale FS
chull_3d(con.pco$li,m=6,prec="Qt") -> ch_st
fric_3d(site.mat1,con.pco$li,m=6,prec="Qt", fric.3d.max = ch_st) -> FRic_co
fdisp_k_sub(tr.dist_con, site.mat1,tax_sub=colnames(site.mat1), m=7)$FDis -> FDis_co
feve_k(con.pco$li,site.mat1,m=7) -> FEve_co

# Substituting NAs by 0 (only if necessary)
FRic_co[which(is.na(FRic_co)==T)]<-0
FEve_co[which(is.na(FEve_co)==T)]<-0

# Comparing FD metrics for the study and regional pools (gamma diversity)
reg_pool<-rbind(reg.mat, rep(1,ncol(reg.mat)))
rownames(reg_pool)[7]<-"study pool"

# FD metrics of study and regional pools  based on overarching FS
#2D
fric_3d(reg_pool,ov.pco$li,m=2,prec="Qt") -> FRic_co_2d
write.table(FRic_co_2d, "FRic_co_2d.txt", sep="\t")
#6D
fric_3d(reg_pool,ov.pco$li,m=6,prec="Qt", fric.3d.max =ch_ov) -> FRic_ov_hull


# FD metrics of study and regional pools based on study-specific FS 
fric_3d(reg.mat,con.pco$li,m=6,prec="Qt", fric.3d.max = ch_st) -> FRic_st
fdisp_k_sub(tr.dist_con, reg.mat,tax_sub=colnames(reg.mat), m=7)$FDis -> FDis_st
feve_k(con.pco$li,reg.mat,m=7) -> FEve_st

# Using Regional-scale FSs

m_fric=6

# Rif
fric_3d(rif.tax,rif.pco$li,m=m_fric,prec="Qt") -> FRic_rif
fdisp_k(tr.dist_rif, rif.tax, m=7) -> FDis_rif
feve_k(rif.pco$li,rif.tax,m=7) -> FEve_rif

# Betic
fric_3d(bet.tax,bet.pco$li,m=m_fric,prec="Qt") -> FRic_bet
fdisp_k(tr.dist_bet, bet.tax, m=7) -> FDis_bet
feve_k(bet.pco$li,bet.tax,m=7) -> FEve_bet

# Picos
fric_3d(pic.tax,pic.pco$li,m=m_fric,prec="QJ") -> FRic_pic
fdisp_k(tr.dist_pic, pic.tax, m=7) -> FDis_pic
feve_k(pic.pco$li,pic.tax,m=7) -> FEve_pic

# Jura
fric_3d(jura.tax,jura.pco$li,m=m_fric,prec="Qt") -> FRic_jura
fdisp_k(tr.dist_jura, jura.tax, m=7) -> FDis_jura
feve_k(jura.pco$li,jura.tax,m=7) -> FEve_jura

# Carpathian
fric_3d(car.tax,car.pco$li,m=m_fric,prec="Qt") -> FRic_car
fdisp_k(tr.dist_car, car.tax, m=7) -> FDis_car
feve_k(car.pco$li,car.tax,m=7) -> FEve_car

# Sweden
fric_3d(swe.tax,swe.pco$li,m=m_fric,prec="Qt") -> FRic_swe
fdisp_k(tr.dist_swe, swe.tax, m=7) -> FDis_swe
feve_k(swe.pco$li,swe.tax,m=7) -> FEve_swe

# Merging regionally-estimated FD data
FRic_re<-c(FRic_rif,FRic_bet,FRic_pic,FRic_jura,FRic_car,FRic_swe)
FDis_re<-c(FDis_rif$FDis,FDis_bet$FDis,FDis_pic$FDis,FDis_jura$FDis,FDis_car$FDis,FDis_swe$FDis)
FEve_re<-c(FEve_rif,FEve_bet,FEve_pic,FEve_jura,FEve_car,FEve_swe)

# Region
reg<-factor(c(rep("R",8),rep("B",10),rep("P",10),rep("J",11),rep("C",11),rep("S",11)))
reg<-factor(reg, levels=c("R","B","P","J","C","S"))

##########################################
#### Correlating FD metrics with latitude
##########################################

dat <- read.table("latitude.txt", header = TRUE, sep = "\t")

dat$lat_wgs84[-9]->lat

# Checking variable distributions
par(mfrow=c(1,1))
hist(lat)

hist(FRic_ov)
hist(FRic_co)
hist(FRic_re)
hist(FDis_ov)
hist(FDis_re)
hist(FEve_ov)
hist(FEve_re)

hist(sqrt(FRic_ov))
hist(sqrt(FRic_re))
hist(sqrt(FRic_co))

# Saving un-transformed data
fric_df<-data.frame(FRic_ov, FRic_co, FRic_re)

# Transformations to improve linearity and reduce skewness
sqrt(FRic_ov)->FRic_ov
sqrt(FRic_co)->FRic_co

hist(FRic_ov)
hist(FRic_co)
hist(FRic_re)

# Ranges of values
range(fric_df$FRic_ov)
range(fric_df$FRic_co)
range(fric_df$FRic_re)

################################################################
# Models relating genus richnes and FD components with latitude
################################################################

################
# Genus richness
################

# Total richness
specnumber(site.mat1)->ric
mod_ric0<-lm(ric~lat)
mod_ric1<-lm(ric~lat+I(lat^2))

model.sel(mod_ric0, mod_ric1)

mod_ric<-lm(ric~lat)

################
# FRic
################

# overarching

mod_fric_ov0<-lm(FRic_ov~lat)
mod_fric_ov1<-lm(FRic_ov~lat+I(lat^2))

model.sel(mod_fric_ov0, mod_fric_ov1)

# continental

mod_fric_co0<-lm(FRic_co~lat)
mod_fric_co1<-lm(FRic_co~lat+I(lat^2))

model.sel(mod_fric_co0, mod_fric_co1)

# regional

mod_FRic_re0<-lm(FRic_re~lat)
mod_FRic_re1<-lm(FRic_re~lat+I(lat^2))

model.sel(mod_FRic_re0, mod_FRic_re1)

# final

mod_fric_ov<-lm(FRic_ov~lat)
mod_fric_co<-lm(FRic_co~lat+I(lat^2))
mod_FRic_re<-lm(FRic_re~lat+I(lat^2))

################
# FDis
################

# overarching

mod_fdis_ov0<-lm(FDis_ov~lat)
mod_fdis_ov1<-lm(FDis_ov~lat+I(lat^2))

model.sel(mod_fdis_ov0, mod_fdis_ov1)

# continental

mod_fdis_co0<-lm(FDis_co~lat)
mod_fdis_co<-lm(FDis_co~lat+I(lat^2))

model.sel(mod_fdis_co0, mod_fdis_co)

# regional

mod_FDis_re0<-lm(FDis_re~lat)
mod_FDis_re1<-lm(FDis_re~lat+I(lat^2))

model.sel(mod_FDis_re0, mod_FDis_re1)

# final

mod_fdis_ov<-lm(FDis_ov~lat+I(lat^2))
mod_fdis_co<-lm(FDis_co~lat+I(lat^2))
mod_FDis_re<-lm(FDis_re~lat+I(lat^2))

################
# FEve
################

# overarching

mod_feve_ov0<-lm(FEve_ov~lat)
mod_feve_ov1<-lm(FEve_ov~lat+I(lat^2))

model.sel(mod_feve_ov0, mod_feve_ov1)

# continental

mod_feve_co0<-lm(FEve_co~lat)
mod_feve_co<-lm(FEve_co~lat+I(lat^2))

model.sel(mod_feve_co0, mod_feve_co)

# regional

mod_FEve_re0<-lm(FEve_re~lat)
mod_FEve_re1<-lm(FEve_re~lat+I(lat^2))

model.sel(mod_FEve_re0, mod_FEve_re1)

# final

mod_feve_ov<-lm(FEve_ov~lat+I(lat^2))
mod_feve_co<-lm(FEve_co~lat+I(lat^2))
mod_FEve_re<-lm(FEve_re~lat+I(lat^2))

###############
# Model results
##############

##############################
### overarching
##############################

res_ric<-summary(mod_ric)

# Checking residuals on continental model
par(mfrow=c(1,2))
plot(fitted(mod_ric),resid(mod_ric),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod_ric),main="Residual normality")
par(mfrow=c(1,1))

res_fric_ov<-summary(mod_fric_ov)

# Checking residuals on continental model
par(mfrow=c(1,2))
plot(fitted(mod_fric_ov),resid(mod_fric_ov),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod_fric_ov),main="Residual normality")
par(mfrow=c(1,1))

res_fdis_ov<-summary(mod_fdis_ov)

# Checking residuals on continental model
par(mfrow=c(1,2))
plot(fitted(mod_fdis_ov),resid(mod_fdis_ov),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod_fdis_ov),main="Residual normality")
par(mfrow=c(1,1))

res_feve_ov<-summary(mod_feve_ov)

# Checking residuals on continental model
par(mfrow=c(1,2))
plot(fitted(mod_feve_ov),resid(mod_feve_ov),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod_feve_ov),main="Residual normality")
par(mfrow=c(1,1))

##############################
### continental
##############################

res_fric_co<-summary(mod_fric_co)

# Checking residuals on continental model
par(mfrow=c(1,2))
plot(fitted(mod_fric_co),resid(mod_fric_co),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod_fric_co),main="Residual normality")
par(mfrow=c(1,1))

res_fdis_co<-summary(mod_fdis_co)

# Checking residuals on continental model
par(mfrow=c(1,2))
plot(fitted(mod_fdis_ov),resid(mod_fdis_ov),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod_fdis_ov),main="Residual normality")
par(mfrow=c(1,1))

res_feve_co<-summary(mod_feve_co)

# Checking residuals on continental model
par(mfrow=c(1,2))
plot(fitted(mod_feve_co),resid(mod_feve_co),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod_feve_co),main="Residual normality")
par(mfrow=c(1,1))

##############################
### regional
##############################

res_FRic_re<-summary(mod_FRic_re)

# Checking residuals on continental model
par(mfrow=c(1,2))
plot(fitted(mod_FRic_re),resid(mod_FRic_re),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod_FRic_re),main="Residual normality")
par(mfrow=c(1,1))

res_FDis_re<-summary(mod_FDis_re)

# Checking residuals on continental model
par(mfrow=c(1,2))
plot(fitted(mod_fdis_ov),resid(mod_fdis_ov),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod_fdis_ov),main="Residual normality")
par(mfrow=c(1,1))

res_FEve_re<-summary(mod_FEve_re)

# Checking residuals on continental model
par(mfrow=c(1,2))
plot(fitted(mod_FEve_re),resid(mod_FEve_re),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod_FEve_re),main="Residual normality")
par(mfrow=c(1,1))

############################
# storing model results
############################

metric_name<-c(rep("ric", 2),
  rep("FRic_ov", 2), rep("FRic_co", 3),rep("FRic_re", 3),
  rep("FDis_ov", 3), rep("FDis_co", 3),rep("FDis_re", 3),
  rep("FEve_ov", 3), rep("FEve_co", 3),rep("FEve_re", 3))

mod_terms<-c("Intercept", "Latitude","Intercept", "Latitude", rep(c("Intercept", "Latitude", "Latitude^2"), 8))

r2<-c(rep(res_ric$r.squared, 2),
      rep(res_fric_ov$r.squared, 2), rep(res_fric_co$r.squared, 3), rep(res_FRic_re$r.squared, 3),
      rep(res_fdis_ov$r.squared, 3), rep(res_fdis_co$r.squared, 3), rep(res_FDis_re$r.squared, 3),
      rep(res_feve_ov$r.squared, 3), rep(res_feve_co$r.squared, 3), rep(res_FEve_re$r.squared, 3))

res<-data.frame(metric=metric_name,mod_terms,
      rbind(res_ric$coef,
      res_fric_ov$coef, res_fric_co$coef, res_FRic_re$coef,
      res_fdis_ov$coef, res_fdis_co$coef, res_FDis_re$coef,
      res_feve_ov$coef, res_feve_co$coef, res_FEve_re$coef), r2)

write.table(res, "model_res_all.txt", sep="\t", row.names = F)

############################
# Null models
############################

# Null model for overarching-FS data
null_res_ov<-null.func.test(taxa=site.mat1,
               traits=traits,
               runs=999,
               m=10,
               m_fric=6,
               fric_type = "lin",
               prec="Qt")

# Null model for continental-FS data
null_res_co<-null.func.test(taxa=site.mat1,
                         traits=traits.con,
                         runs=999,
                         m=7,
                         m_fric=6,
                         prec="Qt")

# Null model for regional-FS data
null_res_re<-null.func.test.reg(taxa=site.mat1,
                           traits=traits,
                           runs=999,
                           m=7,
                           m_fric=6,
                           prec="Qt",
                           reg.mat)
# Aggregating data
null.res<-data.frame(FS=rep(c("Overarching-FS", "Continental-FS", "Regional-FS"), each=nrow(null_res_ov$pval.res)),
                     rbind(null_res_ov$pval.res, null_res_co$pval.res, null_res_re$pval.res))

# saving results
write.table(null.res,"null_mod_res.txt",sep="\t")

################################
# Correlation and model plots
################################

# Latitude sequence
lat.seq<-seq(min(lat),max(lat),length.out = 1000)

pdf(file="cor_plots.pdf",onefile=T,width=10,height=6.5)

par(mfrow=c(2,3),cex=1.35,cex.axis=0.8,cex.lab=0.9,mar=c(3, 5, 3, 1))

c("gold", "darkorange","red", "deeppink", "steelblue1", "royalblue4")->col.reg
line.let= 0.95
adj.let=-0.4
cex.let=1.4
line.cor=-6

#### Continental vs. biased

m<-lm((FRic_ov)^2~(FRic_co)^2)

m.coef<-as.numeric(round(coef(m),3))
f<-paste("y= ", m.coef[1],if(m.coef[2]>0) " + " else " - ", abs(m.coef[2]),"x", sep="")
cor1 <- bquote(italic(r^2) == .(round(summary(m)$r.squared, 2)))

plot((FRic_ov)^2~(FRic_co)^2, ylab="FRic - overarching", xlab="FRic - continental",
     pch=16,col=col.reg[reg], ylim=c(0,0.4))
abline(m,lwd=3,col="blue")
abline(coef=c(0,1),lwd=2,lty=2)
mtext("a)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(cor1,  line = -1.1, adj = 0.1, cex = 1.2, font = 2)
mtext(f)

# FDis cor

m<-lm(FDis_ov~FDis_co)

m.coef<-as.numeric(round(coef(m),3))
f<-paste("y= ", m.coef[1],if(m.coef[2]>0) " + " else " - ", abs(m.coef[2]),"x", sep="")
cor1 <- bquote(italic(r^2) == .(round(summary(m)$r.squared, 2)))

plot(FDis_ov~FDis_co, ylab="FDis - overarching", xlab="FDis - continental",
     pch=16,col=col.reg[reg])
abline(m,lwd=3,col="blue")
abline(coef=c(0,1),lwd=2,lty=2)
mtext("b)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(cor1,  line=line.cor, adj = 0.9, cex = 1.2, font = 2)
mtext(f)

# FEve cor

m<-lm(FEve_ov~FEve_co)

m.coef<-as.numeric(round(coef(m),3))
f<-paste("y= ", m.coef[1],if(m.coef[2]>0) " + " else " - ", abs(m.coef[2]),"x", sep="")
cor1 <- bquote(italic(r^2) == .(round(summary(m)$r.squared, 2)))

plot(FEve_ov~FEve_co, ylab="FEve - overarching", xlab="FEve - continental",
     pch=16,col=col.reg[reg])
abline(m,lwd=3,col="blue")
abline(coef=c(0,1),lwd=2,lty=2)
mtext("c)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(cor1,  line = line.cor, adj = 0.9, cex = 1.2, font = 2)
mtext(f)

#### Continental vs. regional

m<-lm((FRic_ov)^2~(FRic_re))

m.coef<-as.numeric(round(coef(m),3))
f<-paste("y= ", m.coef[1],if(m.coef[2]>0) " + " else " - ", abs(m.coef[2]),"x", sep="")
cor1 <- bquote(italic(r^2) == .(round(summary(m)$r.squared, 2)))

plot((FRic_re),(FRic_ov)^2, ylab="FRic - overarching", xlab="FRic - reg.",
     pch=16,col=col.reg[reg], ylim=c(0,0.4))
abline(m,lwd=3,col="blue")
abline(coef=c(0,1),lwd=2,lty=2)
mtext("d)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(cor1,  line = -1.1, adj = 0.1, cex = 1.2, font = 2)
mtext(f)

# FDis cor

m<-lm(FDis_ov~FDis_re)

m.coef<-as.numeric(round(coef(m),3))
f<-paste("y= ", m.coef[1],if(m.coef[2]>0) " + " else " - ", abs(m.coef[2]),"x", sep="")
cor1 <- bquote(italic(r^2) == .(round(summary(m)$r.squared, 2)))

plot(FDis_ov~FDis_re, ylab="FDis - overarching", xlab="FDis - reg.",
     pch=16,col=col.reg[reg])
abline(m,lwd=3,col="blue")
abline(coef=c(0,1),lwd=2,lty=2)
mtext("e)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(cor1,  line = line.cor, adj = 0.9, cex = 1.2, font = 2)
mtext(f)

# FEve cor

m<-lm(FEve_ov~FEve_re)

m.coef<-as.numeric(round(coef(m),3))
f<-paste("y= ", m.coef[1],if(m.coef[2]>0) " + " else " - ", abs(m.coef[2]),"x", sep="")
cor1 <- bquote(italic(r^2) == .(round(summary(m)$r.squared, 2)))

plot(FEve_ov~FEve_re, ylab="FEve - overarching", xlab="FEve - reg.",
     pch=16,col=col.reg[reg])
abline(m,lwd=3,col="blue")
abline(coef=c(0,1),lwd=2,lty=2)
mtext("f)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(cor1,  line = line.cor, adj = 0.9, cex = 1.2, font = 2)
mtext(f)

dev.off()

################### Models

pdf(file="mod_plots_all.pdf",onefile=T,width=9,height=11)

par(mfrow=c(3,3),cex=1.35,cex.axis=0.8,cex.lab=0.9,mar=c(4, 4, 4, 1))

c("gold", "darkorange","red", "deeppink", "steelblue1", "royalblue4")->col.reg
line.let= 0.95
adj.let=-0.4
cex.let=1.4
line.cor=-6

# FRic - overarching
plot(lat,(FRic_ov)^2, ylab="FRic", xlab = "", pch=16, col=col.reg[reg],
     xaxt = "n")
axis(1,c(40,50,60))
# Fitted values
lines(lat.seq,(predict(mod_fric_ov,data.frame(lat=lat.seq)))^2, lwd=3, col="blue")
#Null model
lines(lat.seq,(pred.lat(null_res_ov$null.fric, lat.seq, type="linear"))^2, lwd=3, col="grey", lty=2)
r2 <- bquote(italic(r)^2 == .(round(summary(mod_fric_ov)$r.squared,2)))
mtext("a)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)

# FDis - overarching
plot(lat,FDis_ov,ylab="FDis", xlab = "", pch=16, col=col.reg[reg], xaxt = "n")
axis(1,c(40,50,60))
lines(lat.seq,predict(mod_fdis_ov,data.frame(lat=lat.seq)), lwd=3, col="blue")
#Null model
lines(lat.seq,pred.lat(null_res_ov$null.fdis, lat.seq), lwd=3, col="grey", lty=2)
r2 <- bquote(italic(r)^2 == .(round(summary(mod_fdis_ov)$r.squared,2)))
mtext("b)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("Overarching", line = 2, adj = 0.5, cex = 1.75, font = 2)

# FEve - overarching
plot(lat,FEve_ov, ylab="FEve", xlab = "", pch=16, col=col.reg[reg], xaxt = "n")
axis(1,c(40,50,60))
lines(lat.seq,predict(mod_feve_ov,data.frame(lat=lat.seq)), lwd=3, col="blue")
#Null model
lines(lat.seq,pred.lat(null_res_ov$null.feve, lat.seq), lwd=3, col="grey", lty=2)

r2 <- bquote(italic(r)^2 == .(round(summary(mod_feve_ov)$r.squared,2)))
mtext("c)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(r2, line = -5.5, adj = 0.9, cex = 1.2, font = 2)

# FRic - continental
plot(lat,(FRic_co)^2, ylab="FRic", xlab = "", pch=16, col=col.reg[reg],
     xaxt = "n")
axis(1,c(40,50,60))
# Fitted values
lines(lat.seq,(predict(mod_fric_co,data.frame(lat=lat.seq)))^2, lwd=3, col="blue")
#Null model
lines(lat.seq,(pred.lat(null_res_co$null.fric, lat.seq))^2, lwd=3, col="grey", lty=2)

r2 <- bquote(italic(r)^2 == .(round(summary(mod_fric_co)$r.squared,2)))
mtext("d)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)

# FDis - continental
plot(lat,FDis_co, ylab="FDis", xlab = "", pch=16, col=col.reg[reg], xaxt = "n")
axis(1,c(40,50,60))
lines(lat.seq,predict(mod_fdis_co,data.frame(lat=lat.seq)), lwd=3, col="blue")
#Null model
lines(lat.seq,pred.lat(null_res_co$null.fdis, lat.seq), lwd=3, col="grey", lty=2)
r2 <- bquote(italic(r)^2 == .(round(summary(mod_fdis_co)$r.squared,2)))
mtext("e)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(r2, line = -1.1, adj = 0.9, cex = 1.2, font = 2)
mtext("Continental", line = 2, adj = 0.5, cex = 1.75, font = 2)

# FEve - continental
plot(lat,FEve_co, ylab="FEve", xlab = "", pch=16, col=col.reg[reg], xaxt = "n")
axis(1,c(40,50,60))
lines(lat.seq,predict(mod_feve_co,data.frame(lat=lat.seq)), lwd=3, col="blue")
#Null model
lines(lat.seq,pred.lat(null_res_co$null.feve, lat.seq), lwd=3, col="grey", lty=2)
r2 <- bquote(italic(r)^2 == .(round(summary(mod_feve_co)$r.squared,2)))
mtext("f)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(r2, line = -5.5, adj = 0.9, cex = 1.2, font = 2)

# FRic - regional
plot(lat,(FRic_re), ylab="FRic", xlab = "", pch=16, col=col.reg[reg],
     xaxt = "n", ylim=c(0,0.7))
axis(1,c(40,50,60))
# Fitted values
lines(lat.seq,(predict(mod_FRic_re,data.frame(lat=lat.seq))), lwd=3, col="blue")
#Null model
lines(lat.seq,(pred.lat(null_res_re$null.fric, lat.seq)), lwd=3, col="grey", lty=2)

r2 <- bquote(italic(r)^2 == .(round(summary(mod_FRic_re)$r.squared,2)))
mtext("g)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(r2, line = -5.5, adj = 0.9, cex = 1.2, font = 2)

# FDis - regional
plot(lat,FDis_re, ylab="FDis", xlab = "Latitude (degrees)", pch=16, col=col.reg[reg], xaxt = "n")
axis(1,c(40,50,60))
lines(lat.seq,predict(mod_FDis_re,data.frame(lat=lat.seq)), lwd=3, col="blue")
#Null model
lines(lat.seq,pred.lat(null_res_re$null.fdis, lat.seq), lwd=3, col="grey", lty=2)
r2 <- bquote(italic(r)^2 == .(round(summary(mod_FDis_re)$r.squared,2)))
mtext("h)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(r2, line = -5.5, adj = 0.9, cex = 1.2, font = 2)
mtext("Regional", line = 2, adj = 0.5, cex = 1.75, font = 2)

# FEve - regional
plot(lat,FEve_re, ylab="FEve", xlab = "", pch=16, col=col.reg[reg], xaxt = "n")
axis(1,c(40,50,60))
lines(lat.seq,predict(mod_FEve_re,data.frame(lat=lat.seq)), lwd=3, col="blue")
#Null model
lines(lat.seq,pred.lat(null_res_re$null.feve, lat.seq), lwd=3, col="grey", lty=2)
r2 <- bquote(italic(r)^2 == .(round(summary(mod_FEve_re)$r.squared,2)))
mtext("i)", line = line.let, adj = adj.let, cex = cex.let, font = 2)
mtext(r2, line = -5.5, adj = 0.9, cex = 1.2, font = 2)

dev.off()

