#' R2SSUR -- Bayesian Sparse Seemingly Unrelated Regression
#' @title runSSUR
#' @description
#' Run a SUR Bayesian sampler
#' @name runSSUR
#' @param dataFile path to data file
#' @param outFilePath path to where the output is to be written
#' @param nIter number of iterations
#' @param nChains number of parallel chains to run
#' @examples
#' require(utils)
#' dir.create("tmp")
#' # *** code to call the function on sample data here ***
#' unlink("tmp", recursive=TRUE)
#' 
#' @export
runSSUR = function(dataFile, blockList, varType=NULL, structureGraph, outFilePath="", 
                nIter=10,  nChains=1, method="SSUR", gammaSampler="Bandit", gammaInit="MLE", usingGPrior=FALSE)
{
  
  # cleanup file PATHS
  dataFileLength = nchar(dataFile)
  if( dataFileLength == 0 )
    stop("Please provide a correct path to a plain-text (.txt) file")
  
  if( substr(dataFile,1,1) == "~" )
    dataFile = path.expand(dataFile)
  
  outFilePathLength = nchar(outFilePath)
  if( outFilePathLength > 0 )
  {
    if( substr(outFilePath,outFilePathLength,outFilePathLength) != "/" )
      paste( outFilePath , "/" , sep="" )
  }
  
  # blockList
  if( length(blockList) < 2 ){
    
    stop("Need at least 2 blocks!")
  
  }else{
    
    blockLabels = rep(NA, max(unlist(blockList)))
    for( i in 1:length(blockList))
      blockLabels[blockList[[i]]] = i-1
    
    blockLabels[is.na(blockLabels)] = -1
    
    # now try and read from given dataFile
    if( !file.exists( dataFile ) ){
      stop("Input file doesn't exists!")      
    }else{
      
      dataHeader = read.table(dataFile,header = FALSE,nrows = 1)
      
      # check the indexes are in the correct range
      if( max(unlist(blockList) ) > length(dataHeader) ){
        stop("blockList indexes provided outside the range in the data matrix!")
      }else{
        if( max(unlist(blockList) ) < length(dataHeader) ){
          blockLabels = c(blockLabels,rep(-1, length(dataHeader) - length(blockLabels) ) )
        }
        # else is fine, it means they're equal
      }
      
    }
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
###### END ######### ZOMBIE CODE! WILL NEED AD SOME POINT WHEN WE INTRODUCE IMPUTATION
  
  ## graph
  if( is.null(dim(structureGraph)) | !is.matrix(structureGraph) ){
    stop("The graph structure must be given as an adjacency matrix")
  }else{
    if( ncol(structureGraph) != length(blockList) ){
      stop("The graph structure must have the same number of blocks as blockList")
    }else{
     
      ## check no loops and at least one covariate-only group
      if( any( ((structureGraph!=0) + (t(structureGraph)!=0))>1 )  )
        stop("The graph should contain no loops nor diagonal elements")
      
      if( !any( rowSums(structureGraph) == 0 ) )
        stop("At least one set of varaible needs to be covariates (right-hand-side) only")
      
      ## else is fine       
    }
  }
  
  dir.create(outFilePath)
  
  dir.create("tmp/")
  write.table(blockLabels,"tmp/blockLabels.txt", row.names = FALSE, col.names = FALSE)
#   write.table(varType,"tmp/varType.txt", row.names = FALSE, col.names = FALSE)
  write.table(structureGraph,"tmp/structureGraph.txt", row.names = FALSE, col.names = FALSE)
  
  status = R2SSUR_internal(dataFile, "tmp/blockLabels.txt", "tmp/structureGraph.txt", outFilePath, nIter,  nChains, method, gammaSampler, gammaInit, usingGPrior)

  if(outFilePath != "tmp/")
    unlink("tmp",recursive = TRUE)
  
  return(status)
}