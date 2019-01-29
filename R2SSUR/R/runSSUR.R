#' R2SSUR -- Bayesian Sparse Seemingly Unrelated Regression
#' @title runSSUR
#' @description
#' Run a SUR Bayesian sampler
#' @name runSSUR
#' @param data path to data file
#' @param outFilePath path to where the output is to be written
#' @param nIter number of iterations
#' @param nChains number of parallel chains to run
#' @examples
#' \donttest{
#' data(example_data, package = "R2SSUR")
#' 
#' R2SSUR::runSSUR(example_data[["data"]],outFilePath = "results/",
#'                 blockList = example_data[["blockList"]], structureGraph = example_data[["structureGraph"]],
#'                 nIter = 20000,nChains = 2, method = "SSUR")
#' 
#' ## check output
#' greyscale = grey((1000:0)/1000)
#' data(example_ground_truth, package = "R2SSUR")
#' 
#' est_gamma = as.matrix( read.table("results/data_SSUR_gamma_out.txt") )
#' est_G = as.matrix( read.table("results/data_SSUR_G_out.txt") )
#' s = ncol(est_G)
#' 
#' par(mfrow=c(2,2))
#' image(est_gamma,col=greyscale); image(example_ground_truth[["gamma"]],col=greyscale)
#' image((est_G+diag(s))[s:1,],col=greyscale); image(example_ground_truth[["G"]][s:1,],col=greyscale)
#' par(mfrow=c(1,1))
#' }
#' 
#' @export
runSSUR = function(data, blockList, varType=NULL, structureGraph=NULL, outFilePath="", 
                nIter=10,  nChains=1, method="SSUR", gammaSampler="Bandit", gammaInit="MLE", usingGPrior=FALSE)
{
  
  dir.create("tmp/")
  
  # is the data given as matrix?
  if( ! (is.character(data) & length(data) == 1) )
  {
    if( is.null(dim(data)) )
      stop("Please provide either the path to a plain-text (.txt) file or a matrix for 'data'")
    
    ## If it's valid matrix, simply write it and re-assign the variable data to hold its path
    write.table(data,"tmp/data.txt", row.names = FALSE, col.names = FALSE)
    data = "tmp/data.txt"
  }
  
  # cleanup file PATHS
  dataLength = nchar(data)
  if( dataLength == 0 )
    stop("Please provide a correct path to a plain-text (.txt) file")
  
  if( substr(data,1,1) == "~" )
    data = path.expand(data)
  
  outFilePathLength = nchar(outFilePath)
  if( outFilePathLength > 0 )
  {
    if( substr(outFilePath,outFilePathLength,outFilePathLength) != "/" )
      paste( outFilePath , "/" , sep="" )
  }
  
  # blockList
  if( ! (is.character(blockList) & length(blockList) == 1) )
  {
    if( length(blockList) < 2 )
      stop("Need at least 2 blocks!")
    
    blockLabels = rep(NA, max(unlist(blockList)))
    for( i in 1:length(blockList))
      blockLabels[blockList[[i]]] = i-1
    
    blockLabels[is.na(blockLabels)] = -1
    
    # now try and read from given data
    if( !file.exists( data ) )
      stop("Input file doesn't exists!")      
    
    dataHeader = read.table(data,header = FALSE,nrows = 1)
    
    # check the indexes are in the correct range
    if( max(unlist(blockList) ) > length(dataHeader) ){
      stop("blockList indexes provided outside the range in the data matrix!")
    }else{
      if( max(unlist(blockList) ) < length(dataHeader) ){
        blockLabels = c(blockLabels,rep(-1, length(dataHeader) - length(blockLabels) ) )
      }
      # else is fine, it means they're equal
    }
    
    blockListLength = length(blockList)
    write.table(blockLabels,"tmp/blockLabels.txt", row.names = FALSE, col.names = FALSE)
    blockList = "tmp/blockLabels.txt"
    
  }else{ # else assume it's already a valid path, should we do checks on this as well?
  
    bl = read.table(blockList)
    blockListLength = length(bl)
    rm(bl)
    
  }
####### CARE ########## ZOMBIE CODE! WILL NEED AD SOME POINT WHEN WE INTRODUCE IMPUTATION
#   # default varType
#   if(is.null(varType)){
#     varType = rep( 0,length(dataHeader) )
#   }else{
    
#     if( length(varType) > length(dataHeader) ){
#       stop("more varible types provided than columns in data matrix")
#     }else{
#       if( length(varType) < length(dataHeader) ){
        
#         if( length(varType) != sum(blockLabels!=-1) ){
#           stop("less varible types provided than used columns in data matrix")
#         }
        
#       }# else is fine
#     }
#   }
#   write.table(varType,"tmp/varType.txt", row.names = FALSE, col.names = FALSE)
###### END ######### ZOMBIE CODE! WILL NEED AD SOME POINT WHEN WE INTRODUCE IMPUTATION
  
  ## structure graph
  if( is.null(structureGraph) )
  {
    cat("structureGraph is null, so I'll assume that the first block correspond to Ys, the second to Xs
        and that there's no other block. Everything else produces an error!\n")
    
    structureGraph = structureGraph = matrix(c(0,0,1,0),2,2,byrow=TRUE)
    
  }
  
  if( ! (is.character(structureGraph) & length(structureGraph) == 1) )
  {
    if( is.null(dim(structureGraph)) | !is.matrix(structureGraph) )
      stop("The graph structure must be given as an adjacency matrix")

    ## check no loops and at least one covariate-only group
    if( any( ((structureGraph!=0) + (t(structureGraph)!=0))>1 )  )
      stop("The graph should contain no loops nor diagonal elements")
    
    if( !any( rowSums(structureGraph) == 0 ) )
      stop("At least one set of varaible needs to be covariates (right-hand-side) only")
    
    if( ncol(structureGraph) != blockListLength )
      stop("The graph structure must have the same number of blocks as blockList")
    
    ## else is fine       
    write.table(structureGraph,"tmp/structureGraph.txt", row.names = FALSE, col.names = FALSE)
    structureGraph = "tmp/structureGraph.txt"
  
  }# else assume it's already a valid path, should we do checks on this as well?
  
  dir.create(outFilePath)
  
  status = R2SSUR_internal(data, blockList, structureGraph, outFilePath, nIter,  nChains, method, gammaSampler, gammaInit, usingGPrior)

  if(outFilePath != "tmp/")
    unlink("tmp",recursive = TRUE)
  
  return(status)
}