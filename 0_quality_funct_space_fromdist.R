#################################################################################################################################
## R function for computing the quality of functional dendrogramm and multidimensional functional spaces                       ##
##                                                                                                                             ##
##    This function is a simplified version of the Appendix S1 associated to Maire et al. 2015 (Global Ecol. and Biogeogr.)    ##
##        (i.e. functional distance matrix provided by user and only UPGMA tree is computed)                                   ##
##                  Given a functional distance matrix, the function computes the quality                                      ## 
##                  (i.e. mean squared-deviation between initial functional distance and standardized distance in the          ##
##                  functional space) for UPGMA dendrogram and all the multidimensional functional spaces                      ##
##                  from 2 to N dimensions (N selected by the user).                                                           ## 
##                  A graphical output illustrating the quality of each functional space is also provided.                     ##
##                                                                                                                             ##
##                                                                                                                             ##
##     INPUTS:                                                                                                                 ##
##                                                                                                                             ##
##      - "dist_funct" : a species x species functional distance matrix. Object should be of class 'dist',                     ##
##                       for instance output of function 'daisy' (package 'cluster') or dist.ktab' from ade4 package.          ##
##                                                                                                                             ##                 
##      - "nbdim" : maximum number of dimensions for multidimensional functional spaces. By default, nbdim=7                   ##
##      		Final number of dimensions depends on the number of positive eigenvalues (after correction) obtained with PCoA     ##
##                                                                                                                             ##
##                                                                                                                             ##
##      - "plot" : character string to set the name of the jpeg file for plots illustrating the quality of functional spaces   ##
##                  NA means no plot                                                                                           ##
##                                                                                                                             ##
##     NB: 1/  high value for 'nbdim' can increase computation time                                                            ##
##         2/ if at least one trait is not numeric, 'metric' should be set as 'Gower'                                          ##
##         3/ if metric=Euclidean, functional traits are scaled (mean=0, sd=1) before computing functional distances           ##
##         4/ R libraries 'ape', 'clue', 'cluster', 'geometry', 'gtools' are required                                          ##
##                                                                                                                             ##
##                                                                                                                             ##
##      OUTPUTS: a list with                                                                                                   ##
##                                                                                                                             ##
##       - $ meanSD : a vector with mean squared deviation values for all the functional spaces tested							           ##
##            names are "t_clustering algorithm" for the best tree and 'm_kD' (with k=2:nbdim) for multidimensional spaces     ##
##                                                                                                                             ##
##       - $ details_funct_space : a list with details about functional spaces                                                 ##
##                                                                                                                             ##
##          - $ mat_coord : coordinates of species in the nbdim multidimensional functional space (PCoA)                       ##
##                                                                                                                             ##
##          - $ upgma_tree : dendrogram built with the UPGMA clustering algorithm                                             ##
##                                                                                                                             ##
##          - $dist_raw and $dist_st : lists with raw and standardized distances between species in each functional space      ##
##                 distance matrices of functional spaces are as names in meanSD (e.g. $dist_raw$t_UPGMA or $dist_raw$m_3D)    ##
##                                                                                                                             ##
##       - a jpeg file in the current working directory with :                                                                 ##
##                      - a barplot showing the meanSD for all functional spaces                                               ##
##                      - 'nbdim' panels (or only 15 if nbdim>15) illustrating the quality of each functional space            ##
##                        Points represent species pairs. Mean squared deviation (mSD) is provided at the top of each panel.	 ##
##                                                                                                                             ##
#################################################################################################################################

quality_funct_space_fromdist <- function( dist_funct,  nbdim=7,   plot="quality_funct_space") 
{

#loading required libraries
require(ape)
require(clue)
require(cluster)
require(geometry)
require(gtools)

################################################################################################################################# 

  # checking data
  if ( ("dist" %in% class(dist_funct) ) ==FALSE )   {  stop(" 'dist_funct' must be of class 'dist'")     }
  if (length(dist_funct)<3)   {  stop(" there must be at least 3 species in 'dist_funct' ")     }
  if (sum(is.na(dist_funct))!=0)   {  stop(" NA are not allowed in 'dist_funct' ")     }
  if (nbdim<2)   {  stop(" 'nbdim' must be higher than 1")     }
  
  # functional distance 
  mat_dissim<-dist_funct
  
  # species names
  nm_sp<-row.names(as.matrix(dist_funct))
  
################################################################
# computing PCoA using Caillez correction
mat_pcoa<-pcoa(mat_dissim, correction="none")

# changing number of dimensions given number of positive eigenvalues
nbdim<-min(nbdim,ncol(mat_pcoa$vectors) )

# keeping species coordoinates on the 'nbdim' axes
mat_coord<-mat_pcoa$vectors[,1:nbdim]
row.names(mat_coord)<-nm_sp
colnames(mat_coord)<-paste("PC",1:nbdim,sep="")

# lists to store distance matrices
dist_raw<-list()
dist_st<-list()

# computing Euclidean distances between species in the (nbdim-1) multidimensionnal functional spaces 
for (k in 2:nbdim) 
    {
    eval(parse(text=paste("dist_",k,"D<-dist(mat_coord[,1:",k,"],method='euclidean')", sep="")))
    eval(parse(text=paste("dist_raw$m_",k,"D<-dist_",k,"D", sep="")))
    } # end of k
  
################################################################ 
# computing mean squared deviation between initial distance and standardized final distance in the functional space
meanSD<-rep(NA,nbdim) ; names(meanSD)<-paste("m_",2:nbdim,"D",sep="")

x<-mat_dissim # initial distance
S<-length(nm_sp) # species richness

 # for muldimensionnal spaces
for (k in 2:nbdim)  
  {
  eval(parse(text=paste("y<-dist_",k,"D",sep="")))
  yst<- y/max(y) * max(x)
  eval(parse(text=paste("dist_st$m_",k,"D<-dist_",k,"D", sep="")))
  meanSD[paste("m_",k,"D",sep="")]<-round( ( (sum((x-yst)^2)) / (S*(S-1)/2) ) ,6)
  }  # end of k

# list of outputs
res<-list(meanSD=meanSD, details_funct_space=list( mat_coord=mat_coord, dist_raw=dist_raw, dist_st=dist_st )  )
  
# plotting change in meanSD with increasing number of dimensions  
par(mfrow=c(1,1),mar=c(5,5,1,1))
barplot(height=meanSD,names.arg=names(meanSD), xlab="Functional space", ylab= "Quality (Mean SD)", 
        space=0, cex.names=0.7, col=c("red", rep("blue",nbdim-1) ) )

abline(h=0.01,lwd=3,lty=2)

# plotting quality of each functional space

# functional distances
x<-as.dist( mat_dissim) 

################################################################################################################################
################################################################################################################################
   
invisible(res)

} # end of function quality_funct_space_fromdist

################################################################################################################################
################################################################################################################################


