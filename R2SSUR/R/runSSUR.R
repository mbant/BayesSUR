#' R2SSUR -- Bayesian Sparse Seemingly Unrelated Regression
#' @title runSSUR
#' @description
#' Run a SUR Bayesian sampler
#' @name runSSUR
#' @param data either a matrix/dataframe or the path to (a plain text) data file with variables on the columns and observations on the rows 
#' @param blockList list of blocks in the model; each element of the list contains the (column) indices of the variables in each block, with respect to the data file
#' @param varType variable type for each column in the data file; coded as: 0 - continuous, 1- binary, 2 - categorical. Note that categorical variables cannot be imputed
#' @param structureGraph graph adjacency matrix representing the structure between the blocks. Edges represented as 2 indicate variables that will be always included in the regression model
#' @param outFilePath path to where the output files are to be written
#' @param nIter number of iterations for the MCMC procedure
#' @param burnin number of iterations (or fraction of iterations) to discard at the start of the chain; Default = 0
#' @param nChains number of parallel chains to run
#' @param method a string indicating the model type, either "SUR" for sparse and dense Seemingly Unrelated Regressions, or "HESS" for hierarchical regression with independent outcomes.
#' @param sparse,dense bool indicating wether the SUR model should have sparse or dense covariance ( not defined for HESS )
#' @param gammaPrior string indicating the gamma prior to use, either "hotspot" for the Hotspot prior of Bottolo (2011), "MRF" for the Markov Random Field prior or "hierarchical" for a simpler hierarchical prior
#' @param gammaSampler string indicating the type of sampler for gamma, either "bandit" for the Thompson sampling inspired samper or "MC3" for the usual $MC^3$ sampler
#' @param gammaInit gamma initialisation to either all-zeros ("0"), all ones ("1"), randomly ("R") or (default) MLE-informed ("MLE").
#' @param mrfG either a matrix or a path to the file containing the G matrix for the MRF prior on gamma (if necessary)
#' @param tmpFilder the path to a temporary folder where intermediate data files are stored (will be erased at the end of the chain); default to local tmpFolder
#'
#' @examples
#' \donttest{
#' data(example_data, package = "R2SSUR")
#' 
#' fit = R2SSUR::runSSUR(example_data[["data"]],outFilePath = "results/",
#'                       blockList = example_data[["blockList"]], structureGraph = example_data[["structureGraph"]],
#'                       nIter = 100, nChains = 2, method = "SUR", gammaPrior = "hotspot" )
#' 
#' ## check output
#' greyscale = grey((1000:0)/1000)
#' data(example_ground_truth, package = "R2SSUR")
#' 
#' est_gamma = as.matrix( read.table(fit$output$gamma) )
#' est_G = as.matrix( read.table(fit$output$G) )
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
                nIter=10, burnin=0, nChains=1, 
                method="SUR", sparse = NULL, dense=NULL,
                gammaPrior="",gammaSampler="bandit", gammaInit="MLE", mrfG=NULL,
                betaPrior="independent",
                tmpFolder="tmp/")
{
  
  if ( !( method %in% c("HESS","SUR") ) )
    stop("Method needs tobe either HESS or SUR")
  methodString = method
  if( methodString == "SUR" )
    methodString = paste( ifelse(sparse,"S","d") , methodString , sep="")
  
  if( ! dir.exists(tmpFolder) )
    dir.create(tmpFolder)
  
  # is the data given as matrix?
  if( ! (is.character(data) & length(data) == 1) )
  {
    if( is.null(dim(data)) )
      stop("Please provide either the path to a plain-text (.txt) file or a matrix for 'data'")
    
    ## If it's valid matrix, simply write it and re-assign the variable data to hold its path
    write.table(data,paste(sep="",tmpFolder,"data.txt"), row.names = FALSE, col.names = FALSE)
    data = paste(sep="",tmpFolder,"data.txt")
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
  
  dataString = head( strsplit( tail( strsplit(data,split = c("/"))[[1]] , 1 ) , ".txt" )[[1]] , 1 ) # magic stripping / and .txt from the data file name
  
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
    write.table(blockLabels,paste(sep="",tmpFolder,"blockLabels.txt"), row.names = FALSE, col.names = FALSE)
    blockList = paste(sep="",tmpFolder,"blockLabels.txt")
    
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
#   write.table(varType,paste(sep="",tmpFolder,"varType.txt"), row.names = FALSE, col.names = FALSE)
#   varType = paste(sep="",tmpFolder,"varType.txt") 
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
    write.table(structureGraph, paste(sep="",tmpFolder,"structureGraph.txt"), row.names = FALSE, col.names = FALSE)
    structureGraph = paste(sep="",tmpFolder,"structureGraph.txt")
  
  }# else assume it's already a valid path, should we do checks on this as well?

  # check how burnin was given
  if ( burnin < 0 ){
    stop("Burnin must be positive or 0")
  }else{ if ( burnin > nIter ){
    stop("Burnin might ont be greater then nIter")
  }else{ if ( burnin < 1 ){ # given as a fraction
    burnin = ceiling(nIter * burnin) # the zero case is taken into account here as well
  }}} # else assume is given as an absolute number

  
  if( is.null(dense) )
  {
    
    if( is.null(sparse))
      sparse = TRUE ## default

    dense = !sparse
    
  }else{
    
      if ( is.null(sparse) )
        sparse = !dense
      
      if( dense == sparse )
      {
        stop("Only sparse or only dense can be set to ",dense)
      }
  }
  
  # mrfG
  if( !(is.character(mrfG) & length(mrfG) == 1) )
  {
    if( is.null(mrfG) )
    {
      mrfG=""
    }else{
      write.table(mrfG,paste(sep="",tmpFolder,"mrfG.txt"), row.names = FALSE, col.names = FALSE)
      mrfG =  paste(sep="",tmpFolder,"mrfG.txt")
    }    
  }

  if( gammaPrior == "" )
	{
		if ( is.null(mrfG) )
		{
			cat( "Using default prior for Gamma - hotspot prior\n");
			mrfG=""
			gammaPrior = "hotspot";
		}
		else
		{
			cat( "No value for gammaPrior was specified, but mrfG was given - choosing MRF prior\n");
			gammaPrior = "MRF";
		}
		
	}
  
  ## Create the directory for the results
  dir.create(outFilePath)
  
  ## Create the return object
  ret = list( status=1, output = list() )
  class(ret) = "R2SSUR"
  
  ret$output["logP"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_logP_out.txt")
  ret$output["gamma"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_gamma_out.txt")
  
  if( gammaPrior %in% c("hierarchical","hotspot"))
    ret$output["pi"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_pi_out.txt")

  if( gammaPrior == "hotspot")
    ret$output["tail"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_hotspot_tail_p_out.txt")
  
  if( method == "SUR" )
  {
    if ( sparse )
      ret$output["G"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_G_out.txt")
    
    ret$output["beta"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_beta_out.txt")
    ret$output["sigmaRho"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_sigmaRho_out.txt")
    
  }
    
  
  ret$status = R2SSUR_internal(data, blockList, structureGraph, outFilePath, nIter, burnin, nChains, 
            method, sparse, gammaPrior, gammaSampler, gammaInit, mrfG, betaPrior)

  if(outFilePath != tmpFolder)
    unlink(tmpFolder,recursive = TRUE)
  
  return(ret)
}
