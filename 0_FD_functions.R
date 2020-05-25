#################################################################################################################
#                                                                                                               #
# Supplementary Data 1                                                                                          #
#                                                                                                               #
# R script to estimate Functional Diversity (FD) metrics                                                        # 
#                                                                                                               #
#                                                                                                               #
# A trait space at an overarching-scale yields more conclusive macroecological patterns of functional diversity  #
#                                                                                                               #
#                                                                                                               #
# Cesc Múrria, Gabone Iturrarte and Cayetano Gutiérrez-Cánovas                                                  #
#                                                                                                               #
# Code written by Cayetano Gutiérrez-Cánovas                                                                    #
# email for queries: cayeguti@um.es                                                                             #
#                                                                                                               #
#################################################################################################################

# prep.fuzzy.df() Functions to prepare matrices, construct the functional space and to estimate functional diversity components
# 
# Inputs:
# traits: fuzzy coding traits
# col.blocks: blocks indicating the number of trait categories per trait
#
# Outputs
# traits: fuzzy coding traits as percentages

prep.fuzzy.df<-function (traits, col.blocks) 
{
  if (!is.data.frame(traits)) 
    stop("Data.frame expected")
  if (sum(col.blocks) != ncol(traits)) {
    stop("Non convenient data in col.blocks")
  }
  if (is.null(names(col.blocks))) {
    names(col.blocks) <- paste("FV", as.character(1:length(col.blocks)), sep = "")
  }
  f1 <- function(x) {
    a <- sum(x)
    if (is.na(a)) 
      return(rep(0, length(x)))
    if (a == 0) 
      return(rep(0, length(x)))
    return(x/a)
  }
  k2 <- 0
  col.w <- rep(1, ncol(traits))
  
  for (k in 1:(length(col.blocks))) {
    k1 <- k2 + 1
    if (col.blocks[k]==1) k2<-k1 else k2 <- k2 + col.blocks[k]
    X <- as.matrix(traits[, k1:k2])
    if (col.blocks[k]==1) X[which(X[,1]>0),]<-1 else X <- t(apply(X, 1, f1))
    X.marge <- apply(X, 1, sum)
    X.marge <- X.marge
    X.marge <- X.marge/sum(X.marge)
    X.mean <- apply(X * X.marge, 2, sum)
    nr <- sum(X.marge == 0)
    cat(nr, "missing data found in block", k, "\n")
    traits[, k1:k2] <- X
    col.w[k1:k2] <- X.mean
  }
  attr(traits, "col.blocks") <- col.blocks
  attr(traits, "col.freq") <- col.w
  col.num <- factor(rep((1:length(col.blocks)), col.blocks))
  attr(traits, "col.num") <- col.num
  return(traits)
}

# chull_3d() estimates the convex hull of a Functional Space
#
# Inputs:
# fpc: functional space
# m: number of axes to select
# prec: convex hull precision ("Qt" or "QJ")
#
# Output:
# a vector with the Functional Richness of each community

chull_3d<-function(fpc,m,prec=c("Qt","QJ")){
  
  convhulln(fpc[,1:m], c("FA",prec))$vol->fric.3d.max
  
  return(fric.3d.max)
  }
  
# fric_3d() estimates the Functional Richness of a set of communties
# This function computes the hypervolume to estimate how each community fills
# the functional space
#
# Inputs:
# taxa: community data
# fpc: functional space
# m: number of axes to select
# prec: convex hull precision ("Qt" or "QJ")
# fric.3d.max: volume of the regional convex hull (set fric.3d.max=1 to have the absolute FRic value)
#
# Output:
# a vector with the Functional Richness of each community

fric_3d<-function(taxa,fpc,m,prec=c("Qt","QJ"), fric.3d.max=NULL){
  fric.3d<-rep(NA,nrow(taxa))
  names(fric.3d)<-rownames(taxa)
  
  # Convex hull of regional pool
  if(is.null(fric.3d.max)==T) convhulln(fpc[,1:m], c("FA",prec))$vol->fric.3d.max
  
  specnumber(taxa)->ric
  
  for (com in 1:nrow(taxa)){
    fpc[which(unlist(rep(taxa[com,]))>0),1:m]->tr.com
    if (ric[com]>=m+1) convhulln(tr.com, c("FA",prec))$vol/fric.3d.max->fric.3d[com] else NA->fric.3d[com]
  }
  return(fric.3d)
}

# fric_1d() estimates the Functional Richness of a set of communties
# This function computes the hypervolume to estimate how each community fills
# the functional space
#
# Inputs:
# taxa: community data
# fpc: functional space
# m: number of axes to select
#
# Output:
# a vector with the Functional Richness of each community

fric_1d<-function(taxa,fpc,m){
  
  fric.1d<-rep(NA,nrow(taxa))
  names(fric.1d)<-rownames(taxa)
  fric.1d.max<-rep(NA,m)
  
  for (i in 1:m) sum(abs(range(fpc[,i])))->fric.1d.max[i]
  
  specnumber(taxa)->ric
  for (com in 1:nrow(taxa)){
    fpc[which(unlist(rep(taxa[com,]))>0),1:m]->tr.com
    
    if (ric[com]>=1) mean(sum(abs(range(tr.com)))/fric.1d.max)->fric.1d[com] else NA->fric.1d[com]
  }
  return(fric.1d)
}

# fdisp_k() estimates the Functional Dispersion of a set of communties
# This function computes the weighted mean distance to community centroid in
# the functional space
# 
# Function modified from Laliberté & Legendre (2010) Ecology
#
# Inputs:
# d: trait dissimilarity matrix
# a: community data
# m: number of axes to select
# tol: tolerance threshold to test whether the distance matrix is Euclidean
#
# Output:
# FDis: a vector with the Functional Dispersion of each community
# eig: eigenvectors of each functional axis
# vectors: functional axes

fdisp_k<-function (d, a, m, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  if (ncol(a) != n) 
    stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  if (any(sn.d != sn.a)) 
    stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
         "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    warning("At least one community has zero-sum abundances (no species).", 
            "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    warning("At least one species does not occur in any community (zero total abundance across all communities).", 
            "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  dimnames(vectors) <- list(colnames(a), NULL)
  pos <- eig > 0
  if (m>0) pos<-c(pos[1:m],rep(F,length(pos)-m))
  
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}

# fdisp_k_sub() estimates the Functional Dispersion of a set of communties
#
# This function computes the weighted mean distance to community centroid in
# the functional space for a subset of taxa (taxonomic matric has less taxa than
# trait dissimilarity matrix)
# 
# Function modified from Laliberté & Legendre (2010) Ecology
#
# Inputs:
# d: trait dissimilarity matrix
# a: community data
# m: number of axes to select
# tol: tolerance threshold to test whether the distance matrix is Euclidean
# tax_sub: taxonomic subset for which FDis will be computed
#
# Output:
# FDis: a vector with the Functional Dispersion of each community
# eig: eigenvectors of each functional axis
# vectors: functional axes

fdisp_k_sub<-function (d, a, tax_sub=NULL, m, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  #if (ncol(a) != (n-length(tax_sub))) 
  #  stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  #if (any(sn.d != sn.a)) 
  #  stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
  #       "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    warning("At least one community has zero-sum abundances (no species).", 
            "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    warning("At least one species does not occur in any community (zero total abundance across all communities).", 
            "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  rownames(vectors) <- attr(d, "Labels")
  vectors <- vectors[(intersect(rownames(vectors), tax_sub)), 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  pos <- eig > 0
  if (m>0) pos<-c(pos[1:m],rep(F,length(pos)-m))
  
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}

# feve_k() estimates the Functional Evenness of a set of communties
# This function computes regularity in the distribution of taxa or abundances
# across the Functional Space, using the Spanning Tree Method
# 
# Function modified from Laliberté & Legendre (2010) Ecology
#
# Inputs:
# fpc: functional space axes
# taxa: community data
# m: number of axes to select
#
# Output:
# FEve: a vector with the Functional Evenness of each community

feve_k<-function(fpc,taxa,m){
  
  rdt=1
  
  # Creating variables
  nrow(taxa)->c
  FEve <- rep(NA, c) ; names(FEve) <- row.names(taxa)
  nb.sp<-specnumber(taxa)
  
  # generating the taxonomic matrix arranged according to the replicated trait values
  tax.pool<-ncol(taxa)
  taxa.rep<-data.frame(matrix(NA,nrow(taxa),tax.pool*rdt))
  spp.list<-c(1:(tax.pool*rdt))
  
  for (spp in 1:tax.pool) {paste(rep("spp",rdt),spp,sep="")->spp.list[((spp-1)*rdt+1):(spp*rdt)]}
  
  colnames(taxa.rep)<-spp.list
  
  for (spp in 1:tax.pool){taxa.rep[,((spp-1)*rdt+1):(spp*rdt)]<-taxa[,spp]/rdt}                             
  
  # Estimating Functional Evenness for each community
  
  for (i in 1:c) {
    sppres <- which(taxa.rep[i, ] > 0)
    # number of species in the community
    S <- length(sppres)
    ab <- as.matrix(taxa.rep[i, sppres])
    # scaling of abundances
    abundrel <- ab / sum(ab)
    
    # selecting the c
    tr <- data.frame(fpc[sppres,1:m ])
    
    if (nb.sp[i] > 2) {
      tr.dist <- dist(tr)
      linkmst <- mst(tr.dist)
      mstvect <- as.dist(linkmst)
      abund2 <- matrix(0, nrow = S, ncol = S)
      for (q in 1:S) for (r in 1:S) abund2[q, r] <- abundrel[q] +  abundrel[r]
      abund2vect <- as.dist(abund2)
      EW <- rep(0, S - 1)
      flag <- 1
      for (mv in 1:((S - 1) * S/2)) {
        if (mstvect[mv] != 0) {
          EW[flag] <- tr.dist[mv]/(abund2vect[mv])
          flag <- flag + 1
        }
      }
      minPEW <- rep(0, S - 1)
      OdSmO <- 1/(S - 1)
      for (l in 1:(S - 1)) minPEW[l] <- min((EW[l]/sum(EW)), 
                                            OdSmO)
      FEve[i] <- ((sum(minPEW)) - OdSmO)/(1 - OdSmO)
    } else FEve[i] <- NA
  }
  return(FEve)
}

# Function to plot 2D functional spaces
#
# coord: coordinates of two FS axes
# col_ch: colour of the FS filling
# col_points: colour of the FS vertices (points)

plot_chull2D<-function(coord,col_ch="#458B0050",border_ch="#458B00"){
  vert0<-convhulln( coord ,"Fx TO 'vert.txt'")
  vert1<-scan("vert.txt",quiet=T)
  vert_ij<-(vert1+1)[-1]
  points(coord[vert_ij,],pch=15, cex=1.5,col=border_ch)
  polygon(coord[vert_ij,], col=col_ch, border=border_ch)
}

# p.val() estimates p-value and the z-score for a given empirical value and null distribution
#
# Inputs:
# obs: empirical value
# null.dist: null distribution of values
# alternative: how empirical and null values were compared. If "two.sided" is selected a bilateral
# comparison is done. "less" provides the probability of obs to be greater than the null values,
# whilst when "greater" the probability of obs to be lower than the null values is estimated
#
# Outputs:
# $z.score: z-score
# $p.value: p-value
# $alternative: how empirical and null values were compared

p.val<-function (obs, null.dist, alternative = c("two.sided", "less", "greater"))
{ 
  alternative <- match.arg(alternative)
  nsimul=length(null.dist)
  z <- (obs - mean(null.dist))/sd(null.dist) # z-statistic
  pless <- sum(obs <= null.dist, na.rm = TRUE) # null values less or equal than real observation
  pmore <- sum(obs >= null.dist, na.rm = TRUE) # null values greater or equal than real observation
  
  p <- switch(alternative, two.sided = 2 * pmin(pless, pmore), less = pless, greater = pmore)
  p <- pmin(1, (p + 1)/(nsimul + 1)) 
  res<-list(z.score=z,p.value=p,alternative=alternative)
  
  return(res)
}

# null.func.test() is an ad hoc function that estimates the null distribution for 
# FRic, FDis and FEve
#
# Inputs:
# taxa: community data
# traits: a matrix with the functional trait types
# runs: number of randomizations
# m: functional axes to select (For FDis and FEve)
# m: functional axes to select For FRic
# prec: method to estimate the convex hull (Qt or QJ)

# Outputs:
#
# ses and pvalue for each model coefficient

null.func.test<-function(taxa,traits,runs=9,m, m_fric, prec="Qt",fric_type="quad"){
  
  res<-list()
  
  # Calculation of the GLOBAL functional distances between species
  
  ktab1_wf <- ktab.list.df(list(traits[which(rowSums(traits)==9),]))
  tr.dist <- dist.ktab(ktab1_wf, type= c("F"))
  
  dudi.pco(tr.dist,nf=m,scannf = F)->tr.pco
  
  funspace <- tr.pco$li
  
  # Using Global PCoA
  chull_3d(funspace,m=m_fric,prec=prec) -> ch_ov
  fric_3d(taxa,funspace,m_fric,prec=prec, fric.3d.max = ch_ov) -> FRic_g # global
  fdisp_k_sub(tr.dist, taxa,tax_sub=colnames(taxa), m)$FDis -> FDis_g # global
  feve_k(funspace,taxa,m) -> FEve_g # global
  
  #if(min(FRic_g, na.rm=T)==0) log(FRic_g+0.01)->FRic_g else log(FRic_g)->FRic_g
  
  sqrt(FRic_g)->FRic_g
  
  # Observed models
  if(fric_type=="quad") mod_fric_g<-lm(FRic_g~lat+I(lat^2)) else mod_fric_g<-lm(FRic_g~lat)
  mod_fdis_g<-lm(FDis_g~lat+I(lat^2))
  mod_feve_g<-lm(FEve_g~lat+I(lat^2))
  
  # Saving observed results
  res$obs.fric<-res$obs.fdis<-res$obs.feve<-rep(NA,3)
  
  names(res$obs.fric)<-names(res$obs.fdis)<-names(res$obs.feve)<-c("intercept","lat","lat^2")
  
  coef(mod_fric_g)->res$obs.fric
  coef(mod_fdis_g)->res$obs.fdis
  coef(mod_feve_g)->res$obs.feve
  
  # Creating an empty matrix for null model coefficients
  res$null.fric<-res$null.fdis<-res$null.feve<-data.frame(matrix(NA,runs,3))
  
  colnames(res$null.fric)<-colnames(res$null.fdis)<-colnames(res$null.feve)<-c("intercept","lat","lat^2")
  
  for (i in 1:runs){
    sample(1:nrow(traits))->ran
    traits[ran,]->traits.r # randomizing traits
    rownames(traits.r)<-rownames(traits)
    
    # Calculation of the GLOBAL functional distances between species
    ktab1_wf <- ktab.list.df(list(traits.r[which(rowSums(traits.r)==9),]))
    tr.dist.r <- dist.ktab(ktab1_wf, type= c("F"))
    
    dudi.pco(tr.dist.r,nf=m,scannf = F)->tr.pco.r
    
    funspace.r <- tr.pco.r$li
    
    fric_3d(taxa, funspace.r,m_fric,prec=prec, fric.3d.max = ch_ov) -> FRic_g_null # global null
    fdisp_k_sub(tr.dist.r, taxa,tax_sub=colnames(taxa), m) -> FDis_g_null # global null
    feve_k(funspace.r,taxa,m) -> FEve_g_null # global null
    
    FRic_g_null[which(is.na(FRic_g_null)==T)]<-0
    FEve_g_null[which(is.na(FEve_g_null)==T)]<-0
    
    sqrt(FRic_g_null)->FRic_g_null
    
    # Null models
    if(fric_type=="quad") mod_fric_g_null<-lm(FRic_g_null~lat+I(lat^2)) else mod_fric_g_null<-lm(FRic_g_null~lat)
    mod_fdis_g_null<-lm(FDis_g_null$FDis~lat+I(lat^2))
    mod_feve_g_null<-lm(FEve_g_null~lat+I(lat^2))
    
    coef(mod_fric_g_null)->res$null.fric[i,1:length(coef(mod_fric_g_null))]
    coef(mod_fdis_g_null)->res$null.fdis[i,]
    coef(mod_feve_g_null)->res$null.feve[i,]
    
    cat("Randomization number", i, "\n")
  }
  
  pval.res<-matrix(NA,9,ncol=2)
  colnames(pval.res)<-c("SES","pval")
  rownames(pval.res)<-c("FRic-intercept", "FRic-lat", "FRic-lat^2",
                        "FDis-intercept","FDis-lat","FDis-lat^2",
                        "FEve-intercept","FEve-lat","FEve-lat^2")
  
  # SES
  p.val(res$obs.fric[1],res$null.fric[,1])$z.score->pval.res[1,1]
  p.val(res$obs.fric[2],res$null.fric[,2])$z.score->pval.res[2,1]
  if(fric_type=="quad") p.val(res$obs.fric[3],res$null.fric[,3])$z.score->pval.res[3,1]
  
  p.val(res$obs.fdis[1],res$null.fdis[,1])$z.score->pval.res[4,1]
  p.val(res$obs.fdis[2],res$null.fdis[,2])$z.score->pval.res[5,1]
  p.val(res$obs.fdis[3],res$null.fdis[,3])$z.score->pval.res[6,1]
  
  p.val(res$obs.feve[1],res$null.feve[,1])$z.score->pval.res[7,1]
  p.val(res$obs.feve[2],res$null.feve[,2])$z.score->pval.res[8,1]
  p.val(res$obs.feve[3],res$null.feve[,3])$z.score->pval.res[9,1]
  
  # p-value
  p.val(res$obs.fric[1],res$null.fric[,1])$p.value->pval.res[1,2]
  p.val(res$obs.fric[2],res$null.fric[,2])$p.value->pval.res[2,2]
  if(fric_type=="quad") p.val(res$obs.fric[3],res$null.fric[,3])$p.value->pval.res[3,2]
  
  p.val(res$obs.fdis[1],res$null.fdis[,1])$p.value->pval.res[4,2]
  p.val(res$obs.fdis[2],res$null.fdis[,2])$p.value->pval.res[5,2]
  p.val(res$obs.fdis[3],res$null.fdis[,3])$p.value->pval.res[6,2]
  
  p.val(res$obs.feve[1],res$null.feve[,1])$p.value->pval.res[7,2]
  p.val(res$obs.feve[2],res$null.feve[,2])$p.value->pval.res[8,2]
  p.val(res$obs.feve[3],res$null.feve[,3])$p.value->pval.res[9,2]
  
  res$pval.res<-pval.res
  
  return(res)
}

# pred.lat() ad-hoc function to predict mean null model values
#
# Input
#
# x: model coefficients (either observed or null)
# seq.lat: sequence of latitude values (independent values)
# type: type of model (linear or quadratic)

pred.lat <- function(x, seq.lat, type="quadratic") {
  
  x <- colMeans(x)
  
  if (type=="linear") res <- (x[1] + x[2]*seq.lat) else res <- (x[1] + x[2]*seq.lat + x[3]*seq.lat^2)
  
  return(as.numeric(res))
}

# null.func.test.reg() is an ad hoc function that estimates the null distribution for 
# FRic, FDis and FEve for regional datasets
#
# Inputs:
# taxa: community data
# traits: a matrix with the functional trait types
# m: functional axes to select (For FDis and FEve)
# m: functional axes to select For FRic
# prec: method to estimate the convex hull (Qt or QJ)
# reg.mat: data.frame with taxon occurrence for each region
# runs: number of randomizations

# Outputs:
#
# ses and pvalue for each model coefficient

null.func.test.reg<-function(taxa,traits,m, m_fric, prec="Qt", reg.mat, runs=9){
  
  res<-list()
  
  # Calculating regional convex hulls
  chull_3d(rif.pco$li,m=m_fric,prec=prec) -> ch_rif
  chull_3d(bet.pco$li,m=m_fric,prec=prec) -> ch_bet
  chull_3d(pic.pco$li,m=m_fric,prec=prec) -> ch_pic
  chull_3d(jura.pco$li,m=m_fric,prec=prec) -> ch_jura
  chull_3d(car.pco$li,m=m_fric,prec=prec) -> ch_car
  chull_3d(swe.pco$li,m=m_fric,prec=prec) -> ch_swe
  
  ch<-c(ch_rif,ch_bet,ch_pic,ch_jura,ch_car,ch_swe)
  
  # Observed regional FD values
  obs.fd<-reg_FD(taxa,traits,m, m_fric, prec=prec, reg.mat, fric.3d.max = ch)
  
  # Assigning to variables
  obs.fd$FRic->FRic_r
  obs.fd$FDis->FDis_r
  obs.fd$FEve->FEve_r
  
  #if(min(FRic_r, na.rm=T)==0) log(FRic_r+0.01)->FRic_r else log(FRic_r)->FRic_r
  sqrt(FRic_r)->FRic_r
  
  # Observed models
  mod_fric_r<-lm(FRic_r~lat+I(lat^2))
  mod_fdis_r<-lm(FDis_r~lat+I(lat^2))
  mod_feve_r<-lm(FEve_r~lat+I(lat^2))
  
  # Saving observed results
  res$obs.fric<-res$obs.fdis<-res$obs.feve<-rep(NA,3)
  
  names(res$obs.fric)<-names(res$obs.fdis)<-names(res$obs.feve)<-c("intercept","lat","lat^2")
  
  coef(mod_fric_r)->res$obs.fric
  coef(mod_fdis_r)->res$obs.fdis
  coef(mod_feve_r)->res$obs.feve
  
  # Creating an empty matrix for null model coefficients
  res$null.fric<-res$null.fdis<-res$null.feve<-data.frame(matrix(NA,runs,3))
  
  colnames(res$null.fric)<-colnames(res$null.fdis)<-colnames(res$null.feve)<-c("intercept","lat","lat^2")
  
  for (i in 1:runs){
    sample(1:nrow(traits))->ran
    traits[ran,]->traits.r # randomizing traits
    rownames(traits.r)<-rownames(traits)
    
    # Observed regional FD values
    null.fd<-reg_FD(taxa,traits.r,m, m_fric, prec=prec, reg.mat, fric.3d.max = ch)
    
    # Assigning to variables
    null.fd$FRic->FRic_r_null
    null.fd$FDis->FDis_r_null
    null.fd$FEve->FEve_r_null
    
    #if(min(FRic_r_null, na.rm=T)==0) log(FRic_r_null+0.01)->FRic_r_null else log(FRic_r_null)->FRic_r_null
    
    sqrt(FRic_r_null)->FRic_r_null
    
    # Null models
    mod_fric_r_null<-lm(FRic_r_null~lat+I(lat^2))
    mod_fdis_r_null<-lm(FDis_r_null~lat+I(lat^2))
    mod_feve_r_null<-lm(FEve_r_null~lat+I(lat^2))
    
    coef(mod_fric_r_null)->res$null.fric[i,]
    coef(mod_fdis_r_null)->res$null.fdis[i,]
    coef(mod_feve_r_null)->res$null.feve[i,]
    
    cat("Randomization number", i, "\n")
  }
  
  pval.res<-matrix(NA,9,ncol=2)
  colnames(pval.res)<-c("SES","pval")
  rownames(pval.res)<-c("FRic-intercept", "FRic-lat", "FRic-lat^2",
                        "FDis-intercept","FDis-lat","FDis-lat^2",
                        "FEve-intercept","FEve-lat","FEve-lat^2")
  
  # SES
  p.val(res$obs.fric[1],res$null.fric[,1])$z.score->pval.res[1,1]
  p.val(res$obs.fric[2],res$null.fric[,2])$z.score->pval.res[2,1]
  p.val(res$obs.fric[3],res$null.fric[,3])$z.score->pval.res[3,1]
  
  p.val(res$obs.fdis[1],res$null.fdis[,1])$z.score->pval.res[4,1]
  p.val(res$obs.fdis[2],res$null.fdis[,2])$z.score->pval.res[5,1]
  p.val(res$obs.fdis[3],res$null.fdis[,3])$z.score->pval.res[6,1]
  
  p.val(res$obs.feve[1],res$null.feve[,1])$z.score->pval.res[7,1]
  p.val(res$obs.feve[2],res$null.feve[,2])$z.score->pval.res[8,1]
  p.val(res$obs.feve[3],res$null.feve[,3])$z.score->pval.res[9,1]
  
  # p-value
  p.val(res$obs.fric[1],res$null.fric[,1])$p.value->pval.res[1,2]
  p.val(res$obs.fric[2],res$null.fric[,2])$p.value->pval.res[2,2]
  p.val(res$obs.fric[3],res$null.fric[,3])$p.value->pval.res[3,2]
  
  p.val(res$obs.fdis[1],res$null.fdis[,1])$p.value->pval.res[4,2]
  p.val(res$obs.fdis[2],res$null.fdis[,2])$p.value->pval.res[5,2]
  p.val(res$obs.fdis[3],res$null.fdis[,3])$p.value->pval.res[6,2]
  
  p.val(res$obs.feve[1],res$null.feve[,1])$p.value->pval.res[7,2]
  p.val(res$obs.feve[2],res$null.feve[,2])$p.value->pval.res[8,2]
  p.val(res$obs.feve[3],res$null.feve[,3])$p.value->pval.res[9,2]
  
  res$pval.res<-pval.res
  
  return(res)
}

# reg_FD() is an ad hoc function that estimates FRic, FDis and FEve 
# for regional datasets
#
# Inputs:
# taxa: community data
# traits: a matrix with the functional trait types
# m: functional axes to select (For FDis and FEve)
# m: functional axes to select For FRic
# prec: method to estimate the convex hull (Qt or QJ)
# reg.mat: data.frame with taxon occurrence for each region

# Outputs:
#
# FRic: functional richness
# FDis: functional dispersion
# FEve: functional evenness

reg_FD<-function(taxa,traits,m, m_fric, prec="Qt", reg.mat, fric.3d.max){
  
  ###########################################################################
  ### Rif ###
  
  # Selecting taxa from the regional pool
  rif.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[1,]>0)]))),]
  
  # Gower dissimilarity
  ktab1_rif <- ktab.list.df(list(rif.tr[which(rowSums(rif.tr)==9),]))
  distrait_rif <- dist.ktab(ktab1_rif, type= c("F"))
  
  # FS (PCoA)
  rif.pco <- dudi.pco(distrait_rif, scan = F, nf = m)
  rif.tax <- taxa[1:8,(intersect(rownames(rif.pco$li), colnames(taxa)))]
  
  ###########################################################################
  ### Betic ###
  
  # Selecting reg.mat from the regional pool
  bet.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[2,]>0)]))),]
  
  # Gower dissimilarity
  ktab1_bet <- ktab.list.df(list(bet.tr[which(rowSums(bet.tr)==9),]))
  distrait_bet <- dist.ktab(ktab1_bet, type= c("F"))
  
  # FS (PCoA)
  bet.pco <- dudi.pco(distrait_bet, scan = F, nf = m)
  bet.tax <- taxa[9:18,(intersect(rownames(bet.pco$li), colnames(taxa)))]
  
  ###########################################################################
  ### Picos ###
  
  # Selecting taxa from the regional pool
  pic.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[3,]>0)]))),]
  
  # Gower dissimilarity
  ktab1_pic <- ktab.list.df(list(pic.tr[which(rowSums(pic.tr)==9),]))
  distrait_pic <- dist.ktab(ktab1_pic, type= c("F"))
  
  # FS (PCoA)
  pic.pco <- dudi.pco(distrait_pic, scan = F, nf = m)
  pic.tax <- taxa[19:28,(intersect(rownames(pic.pco$li), colnames(taxa)))]
  
  ###########################################################################
  ### Jura ###
  
  # Selecting taxa from the regional pool
  jura.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[4,]>0)]))),]
  
  # Gower dissimilarity
  ktab1_jura <- ktab.list.df(list(jura.tr[which(rowSums(jura.tr)==9),]))
  distrait_jura <- dist.ktab(ktab1_jura, type= c("F"))
  
  # FS (PCoA)
  jura.pco <- dudi.pco(distrait_jura, scan = F, nf = m)
  jura.tax <- taxa[29:39,(intersect(rownames(jura.pco$li), colnames(taxa)))]
  
  ###########################################################################
  ### Carpathian ###
  
  # Selecting taxa from the regional pool
  car.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[5,]>0)]))),]
  
  # Gower dissimilarity
  ktab1_car <- ktab.list.df(list(car.tr[which(rowSums(car.tr)==9),]))
  distrait_car <- dist.ktab(ktab1_car, type= c("F"))
  
  # FS (PCoA)
  car.pco <- dudi.pco(distrait_car, scan = F, nf = m)
  car.tax <- taxa[40:50,(intersect(rownames(car.pco$li), colnames(taxa)))]
  
  ###########################################################################
  ### Sweden ###
  
  # Selecting taxa from the regional pool
  swe.tr<-traits[(intersect(rownames(traits), colnames(reg.mat[which(reg.mat[6,]>0)]))),]
  
  # Gower dissimilarity
  ktab1_swe <- ktab.list.df(list(swe.tr[which(rowSums(swe.tr)==9),]))
  distrait_swe <- dist.ktab(ktab1_swe, type= c("F"))
  
  # FS (PCoA)
  swe.pco <- dudi.pco(distrait_swe, scan = F, nf = m)
  swe.tax <- taxa[51:61,(intersect(rownames(swe.pco$li), colnames(taxa)))]
  
  # Using Regional-scale FS
  
  # Regional convex hulls (gamma)
  ch_rif <- fric.3d.max[1]
  ch_bet <- fric.3d.max[2]
  ch_pic <- fric.3d.max[3]
  ch_jura <- fric.3d.max[4]
  ch_car <- fric.3d.max[5]
  ch_swe <- fric.3d.max[6]
  
  # Rif
  fric_3d(rif.tax,rif.pco$li,m=m_fric,prec=prec, fric.3d.max = ch_rif) -> FRic_rif
  fdisp_k(distrait_rif, rif.tax, m) -> FDis_rif
  feve_k(rif.pco$li,rif.tax,m) -> FEve_rif
  
  # Betic
  fric_3d(bet.tax,bet.pco$li,m=m_fric,prec=prec, fric.3d.max = ch_bet) -> FRic_bet
  fdisp_k(distrait_bet, bet.tax, m) -> FDis_bet
  feve_k(bet.pco$li,bet.tax,m) -> FEve_bet
  
  # Picos
  fric_3d(pic.tax,pic.pco$li,m,prec=prec, fric.3d.max = ch_pic) -> FRic_pic
  fdisp_k(distrait_pic, pic.tax, m) -> FDis_pic
  feve_k(pic.pco$li,pic.tax,m) -> FEve_pic
  
  # Jura
  fric_3d(jura.tax,jura.pco$li,m=m_fric,prec=prec, fric.3d.max = ch_jura) -> FRic_jura
  fdisp_k(distrait_jura, jura.tax, m) -> FDis_jura
  feve_k(jura.pco$li,jura.tax,m) -> FEve_jura
  
  # Carpathian
  fric_3d(car.tax,car.pco$li,m=m_fric,prec=prec, fric.3d.max = ch_car) -> FRic_car
  fdisp_k(distrait_car, car.tax, m) -> FDis_car
  feve_k(car.pco$li,car.tax,m) -> FEve_car
  
  # Sweden
  fric_3d(swe.tax,swe.pco$li,m=m_fric,prec=prec, fric.3d.max = ch_swe) -> FRic_swe
  fdisp_k(distrait_swe, swe.tax, m) -> FDis_swe
  feve_k(swe.pco$li,swe.tax,m) -> FEve_swe
  
  # Joining regionally-estimated data
  FRic_r<-c(FRic_rif,FRic_bet,FRic_pic,FRic_jura,FRic_car,FRic_swe)
  FDis_r<-c(FDis_rif$FDis,FDis_bet$FDis,FDis_pic$FDis,FDis_jura$FDis,FDis_car$FDis,FDis_swe$FDis)
  FEve_r<-c(FEve_rif,FEve_bet,FEve_pic,FEve_jura,FEve_car,FEve_swe)
  reg.tr<-list(rif.tr,bet.tr,pic.tr,jura.tr,car.tr,swe.tr)
  
  return(list(FRic=FRic_r, FDis=FDis_r, FEve=FEve_r, reg.tr=reg.tr))
  
}

