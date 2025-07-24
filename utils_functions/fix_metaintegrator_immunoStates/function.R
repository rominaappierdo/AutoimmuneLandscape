# DecLinearRegression <- function(xMat,yVector){  #basisMatrixSample, geneExprsample
#   print(head(yVector))
#   #run linear model without intercept
#   model    <- stats::lm(yVector ~ xMat -1)
#   print(summary(model))
#   #go from coefficients to proportions
#   coeff    <- model$coefficients
#   coeff[coeff < 0] <- 0
#   props <- coeff/sum(coeff)
# 
#   #get RMSE and Correlation for deconvolution
#   dec_r    <- stats::cor(model$fitted.values,yVector)
#   dec_rmse <- Metrics::rmse(model$fitted.values,yVector)
# 
#   #return final list in output
#   newList <- list("props" = props, "rmse" = dec_rmse, "r" = dec_r)
# 
# }
# 
# iSdeconvolution<- function(basisMatrix, exprMat){
#   #format basis matrix and expression data
#   basisMatrix <- data.matrix(basisMatrix)
#   exprMat <- data.matrix(exprMat)
#   print(head(exprMat))
#   #order by row-names
#   basisMatrix <- basisMatrix[order(rownames(basisMatrix)),] #genes as row names
#   #cell types as col names
#   exprMat <- exprMat[order(rownames(exprMat)),]
# 
#   #convert into real [non-log] space if in log space [be wary of NAs]
#   if(max(exprMat, na.rm=T) < 50){
#     print('is log scale')
#     exprMat <- 2^exprMat
#   }
#   #run quantile normalization
#   colN <- colnames(exprMat)
#   rowN <- rownames(exprMat)
# 
#   exprMat <- preprocessCore::normalize.quantiles(exprMat)
#   colnames(exprMat) <- colN
#   rownames(exprMat) <- rowN
#   #only pick genes in common between data and basis matrix
#   bMgenes  <- row.names(basisMatrix)
#   gMgenes  <- row.names(exprMat)
#   GMintBM  <- gMgenes %in% bMgenes
# 
#   if(sum(GMintBM)==0) {
#     stop("None of the gene names are present in the the basis matrix. Make sure your gene names are correct and standard.")
#   }
# 
#   exprMat <- exprMat[GMintBM,]
#   BMintGM <- bMgenes %in% row.names(exprMat)
#   basisMatrix <- basisMatrix[BMintGM,]
# 
#   #standardize basis matrix [rescale globally]
#   basisMatrix <- (basisMatrix-mean(basisMatrix,na.rm=T))/stats::sd(as.vector(basisMatrix),na.rm=T)
# 
#   #declare header for the output matrix
#   header <- c('Sample',colnames(basisMatrix),"P-value","Correlation","RMSE")
# 
#   #declare empty output matrix
#   output   <- matrix()
#   pvalue   <- 9999
# 
#   #run for every sample
#   sampleIndex <- 1
#   while(sampleIndex <= ncol(exprMat)){
#     #iterate one variable @ the time
#     geneExpressionSample <- exprMat[,sampleIndex]
# 
#     #remove NAs [this is necessary to avoid issues downstream]
#     basisMatrixSample    <- basisMatrix[which(!is.na(geneExpressionSample)),]   # sample-specific basis matrix
#     geneExpressionSample <- geneExpressionSample[which(!is.na(geneExpressionSample))]
#     #scale cell-mixture sample
#     geneExpressionSample <- scale(geneExpressionSample)
# 
#     #run linear regression on a single sample
#     decOut  <- DecLinearRegression(basisMatrixSample,geneExpressionSample)
# 
#     #create output vector
#     out <- c(colnames(exprMat)[sampleIndex],
#              decOut$props,
#              pvalue,
#              decOut$r,
#              decOut$rmse)
# 
# 
#     #update output matrix
#     if(sampleIndex == 1){
#       output <- out
#     }else{
#       output <- rbind(output, out)
#     }
# 
#     #update sample index
#     sampleIndex <- sampleIndex + 1
#   }
#   #format matrix object containing all results
#   outObj <- rbind(header,output)
#   outObj <- outObj[,-1]
#   outObj <- outObj[-1,]
#   outObj <- matrix(as.numeric(unlist(outObj)),nrow=nrow(outObj))
# 
#   #object is a matrix of samples X cells [+ some output]
#   rownames(outObj) <- colnames(exprMat)
#   colnames(outObj) <- c(colnames(basisMatrix),"P-value","Correlation","RMSE")
# 
#   return(outObj)
# 
# }




immunoStatesDecov<- function(metaObject){
  invisible("immunoStatesMatrix")
  
  #figure out if this is not a unix machine, in which case mc.cores must be 1
  numCores = 10
  get_os <- function() {
    if (.Platform$OS.type == "windows") { 
      "win"
    } else if (Sys.info()["sysname"] == "Darwin") {
      "mac" 
    } else if (.Platform$OS.type == "unix") { 
      "unix"
    } else {
      stop("Unknown OS")
    }
  }
  if(!get_os() %in% c("mac","unix")){numCores = 1}
  
  #compute expression matrices
  iSExpMat <- mclapply(metaObject$originalData,
                       function(i)
                         .extractDataForGenesDT(i,promiscProbes = F),
                       mc.cores = min(numCores,length(metaObject$originalData)))
  #saveRDS(iSExpMat, file = "/mnt/data2/disk_4T/appierdo/dottorato/data/goodness_of_fit/stuff.RDS")  
  #run immunoStates and save output into a data.table
  outPut        <- lapply(iSExpMat[which(!sapply(iSExpMat,
                                                 function(i)
                                                   all(is.na(i))))],
                          function(j){
                            outDT <- as.data.table(iSdeconvolution(immunoStatesMatrix,j),
                                                   keep.rownames = T)
                            
                            outDT[,natural_killer_cell:=CD56bright_natural_killer_cell+CD56dim_natural_killer_cell]
                            outDT[,           monocyte:=CD14_positive_monocyte+CD16_positive_monocyte]
                            outDT[,             B_cell:=naive_B_cell+memory_B_cell+plasma_cell]
                            outDT[,             T_cell:=CD8_positive_alpha_beta_T_cell+CD4_positive_alpha_beta_T_cell+gamma_delta_T_cell]
                            outDT[,        granulocyte:=neutrophil+eosinophil+basophil]
                            
                            return(outDT)
                          })
  names(outPut) <- names(iSExpMat[which(!sapply(iSExpMat,
                                                function(i)
                                                  all(is.na(i))))])
  
  #return output into meta-object
  metaObject$immunoStates <- outPut
  return(metaObject)
}


environment(immunoStatesDecov) <- asNamespace('MetaIntegrator')
assignInNamespace("immunoStatesDecov", immunoStatesDecov, ns = "MetaIntegrator")
