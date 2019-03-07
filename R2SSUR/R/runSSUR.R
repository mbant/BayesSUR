#' R2SSUR -- Bayesian Sparse Seemingly Unrelated Regression
#' @title runSSUR
#' @description
#' Run a SUR Bayesian sampler
#' @name runSSUR
#' @param data either a matrix/dataframe or the path to (a plain text) data file with variables on the columns and observations on the rows 
#' @param blockList list of blocks in the model each element of the list contains the (column) indices of the variables in each block, with respect to the data file
#' @param varType variable type for each column in the data file coded as: 0 - continuous, 1- binary, 2 - categorical. Note that categorical variables cannot be imputed
#' @param structureGraph graph adjacency matrix representing the structure between the blocks. Edges represented as 2 indicate variables that will be always included in the regression model
#' @param outFilePath path to where the output files are to be written
#' @param nIter number of iterations for the MCMC procedure
#' @param burnin number of iterations (or fraction of iterations) to discard at the start of the chain Default = 0
#' @param nChains number of parallel chains to run
#' @param method a string indicating the model type, either "SUR" for sparse and dense Seemingly Unrelated Regressions, or "HESS" for hierarchical regression with independent outcomes.
#' @param sparse,dense bool indicating wether the SUR model should have sparse or dense covariance ( not defined for HESS )
#' @param gammaPrior string indicating the gamma prior to use, either "hotspot" for the Hotspot prior of Bottolo (2011), "MRF" for the Markov Random Field prior or "hierarchical" for a simpler hierarchical prior
#' @param gammaSampler string indicating the type of sampler for gamma, either "bandit" for the Thompson sampling inspired samper or "MC3" for the usual $MC^3$ sampler
#' @param gammaInit gamma initialisation to either all-zeros ("0"), all ones ("1"), randomly ("R") or (default) MLE-informed ("MLE").
#' @param mrfG either a matrix or a path to the file containing the G matrix for the MRF prior on gamma (if necessary)
#' @param tmpFolder the path to a temporary folder where intermediate data files are stored (will be erased at the end of the chain) default to local tmpFolder
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
#' image(est_gamma,col=greyscale) image(example_ground_truth[["gamma"]],col=greyscale)
#' image((est_G+diag(s))[s:1,],col=greyscale) image(example_ground_truth[["G"]][s:1,],col=greyscale)
#' par(mfrow=c(1,1))
#' }
#' 
my_stop = function( msg , tmpFolder )
{
  unlink(tmpFolder,recursive = TRUE)
  stop(msg)
}

#' @export
runSSUR = function(data, Y, X, X_0=NULL,
                varType=NULL, structureGraph=NULL, outFilePath="", 
                nIter=10, burnin=0, nChains=1, 
                covariancePrior="HIW",
                gammaPrior="",gammaSampler="bandit", gammaInit="MLE", mrfG=NULL,
                betaPrior="independent",
                output_gamma = TRUE, output_beta = TRUE, output_G = TRUE, output_sigmaRho = TRUE,
                output_pi = TRUE, output_tail = TRUE, output_model_size = TRUE,
                tmpFolder="tmp/")
{
  
  # Create temporary directory
  if( ! dir.exists(tmpFolder) )
    dir.create(tmpFolder)
  

  ## Check the input: reasoning is that the user provides either
  # a data matrix or a data path-to-file
  #     - in this case Y, X (and X_0) need to be provided as vectors of indexes 
  # if the data matrix is not provided, the user will give 2/3 matrices for Y, X (and X_0)
  #     - in which case we write those in order into a new joint file
  # everything else throws an error

  # we'll check in reverse order, is data NULL?
  if ( is.null( data ) )
  {

    # Y,X (and if there X_0) need to be valid numeric matrices then
    # check Y and X have comfortable number of observations
    if( !is.numeric(Y) | is.null(dim(Y)) )
      my_stop("Y needs to be a valid matrix or data.frame with > 1 column ")
    
    nObservations = nrows(Y)
    if( !is.numeric(X) | is.null(dim(X)) | nrow(X) != nObservations )
      my_stop("X needs to be a valid matrix or data.frame with >= 1 column and the same number of rows of Y")
    
    if ( is.null ( X_0 ) )
      X_0 = matrix(NA,nrow=nObservations,ncol=0)
    
    if( !is.numeric(X_0) | is.null(dim(X_0)) | nrow(X_0) != nObservations )
      my_stop("if provided, X_0 needs to be a valid matrix or data.frame with >= 1 column and the same number of rows of Y")


    # Write the three down in a single data file
    write.table(cbind(Y,X,X_0),paste(sep="",tmpFolder,"data.txt"), row.names = FALSE, col.names = FALSE)
    data = paste(sep="",tmpFolder,"data.txt")

    blockLabels = c( rep(0,ncol(Y)) , rep(0,ncol(X)) , rep(0,ncol(X_0)) )

  }else{ # data is not null, so the user wants to use it to input the data

    # is the data given as matrix?
    ## If it's valid matrix, simply write it and re-assign the variable data to hold its path
    if( is.numeric(data) & !is.null(dim(data))  )
    {
      write.table(data,paste(sep="",tmpFolder,"data.txt"), row.names = FALSE, col.names = FALSE)
      data = paste(sep="",tmpFolder,"data.txt")
    }

    # is the data given as a string?
    if( is.character(data) & length(data) == 1 )
    {    
      if( substr(data,1,1) == "~" )
        data = path.expand(data)

      # now try and read from given data file
      if( !file.exists( data ) )
        my_stop("Input file doesn't exists!",tmpFolder)
    }

    ## at this point data contains the path to a file that exists
    # try and read one line to check dimensions
    dataHeader = read.table(data,header = FALSE,nrows = 1)
    nVariables = ncol(dataHeader)

    ## Y, X (and X_0) should be some fixed variables that needs to be included in the model
    if ( is.null(X_0) )
      X_0 = c()

    # be sure they can be converted to numeric
    Y = as.numeric(Y)
    X = as.numeric(X)
    X_0 = as.numeric(X_0)

    # be sure that they are vectors
    if ( !( is.vector(Y,"numeric") & is.vector(X,"numeric") & is.vector(X_0,"numeric") ) )
      my_stop("When the `data` argument is set, Y,X and X_0 need to be corresponding index vectors")

    # check thay do not overlap
    if ( length( c( intersect(Y,X) , intersect(Y,X_0) , intersect(X_0,X) ) ) != 0 )
      my_stop("Y,X and X_0 need to be distinct index vectors")

    # check if dimensions correspond -- higher dimensions gets an error
    if( length( c(Y,X,X_0) ) > nVariables )
      my_stop("When the `data` argument is set, Y,X and X_0 need to be corresponding index vectors")
    # equal dimensions are ok, but lower dimensions means some columns of the data will be disregarded ( set to -1  )

    # We can now init the blockList
    blockLabels = rep(NA,nVariables)
    blockLabels[Y] = 0
    blockLabels[X] = 1
    if ( length ( X_0 ) > 0 ) 
      blockLabels[X_0] = 2

    blockLabels[ is.na ( blockLabels ) ] = -1

  }

  ## Then init the structure graph
  # Consider that the indexes are written so that Y is 0 , X is 1 and (if there) X_0 is 2
  if ( length ( X_0 ) > 0 )
    structureGraph = structureGraph = matrix(c(0,0,0,1,0,0,2,0,0),3,3,byrow=TRUE)
  else structureGraph = structureGraph = matrix(c(0,0,1,0),2,2,byrow=TRUE)

  ## Finally write blockLabels and structureGraph to a file
  write.table(blockLabels,paste(sep="",tmpFolder,"blockLabels.txt"), row.names = FALSE, col.names = FALSE)
  blockList = paste(sep="",tmpFolder,"blockLabels.txt")
  write.table(structureGraph, paste(sep="",tmpFolder,"structureGraph.txt"), row.names = FALSE, col.names = FALSE)
  structureGraph = paste(sep="",tmpFolder,"structureGraph.txt")


  # cleanup file PATHS
  dataLength = nchar(data)
  if( dataLength == 0 )
    my_stop("Please provide a correct path to a plain-text (.txt) file", tmpFolder)
  
  outFilePathLength = nchar(outFilePath)
  if( outFilePathLength > 0 )
  {
    if( substr(outFilePath,outFilePathLength,outFilePathLength) != "/" )
      paste( outFilePath , "/" , sep="" )
  }

  # magicly strip '/' from the start and '.txt' from the end of the data file name
  dataString = head( strsplit( tail( strsplit(data,split = c("/"))[[1]] , 1 ) , ".txt" )[[1]] , 1 ) 

####### CARE ########## ZOMBIE CODE! WILL NEED AD SOME POINT WHEN WE INTRODUCE IMPUTATION
#   # default varType
#   if(is.null(varType)){
#     varType = rep( 0,length(dataHeader) )
#   }else{
    
#     if( length(varType) > length(dataHeader) ){
#       my_stop("more varible types provided than columns in data matrix",tmpFolder)
#     }else{
#       if( length(varType) < length(dataHeader) ){
        
#         if( length(varType) != sum(blockLabels!=-1) ){
#           my_stop("less varible types provided than used columns in data matrix",tmpFolder)
#         }
        
#       }# else is fine
#     }
#   }
#   write.table(varType,paste(sep="",tmpFolder,"varType.txt"), row.names = FALSE, col.names = FALSE)
#   varType = paste(sep="",tmpFolder,"varType.txt") 
###### END ######### ZOMBIE CODE! WILL NEED AD SOME POINT WHEN WE INTRODUCE IMPUTATION

  
  # check how burnin was given
  if ( burnin < 0 ){
    my_stop("Burnin must be positive or 0",tmpFolder)
  }else{ if ( burnin > nIter ){
    my_stop("Burnin might ont be greater then nIter",tmpFolder)
  }else{ if ( burnin < 1 ){ # given as a fraction
    burnin = ceiling(nIter * burnin) # the zero case is taken into account here as well
  }}} # else assume is given as an absolute number

  ###############################

  # method to use
  if ( covariancePrior == "sparse" | covariancePrior == "Sparse" | covariancePrior == "SPARSE" | covariancePrior == "HIW" | covariancePrior == "hiw" )
    covariancePrior = "HIW"
  else if ( covariancePrior == "dense" | covariancePrior == "Dense" | covariancePrior == "DENSE" | covariancePrior == "IW" | covariancePrior == "iw" ) 
    covariancePrior = "IW"
  else if ( covariancePrior == "independent" | covariancePrior == "Independent" | covariancePrior == "INDEPENDENT" | covariancePrior == "INDEP" | covariancePrior == "indep" | covariancePrior == "IG" | covariancePrior == "ig" ) 
    covariancePrior = "IG"
  else
    my_stop("Unknown covariancePrior argument: only sparse (HIW), dense(IW) or independent (IG) are available",tmpFolder)
  
  # mrfG and gammaPrior
  if( gammaPrior == "" )
	{
		if ( is.null(mrfG) )
		{
			cat( "Using default prior for Gamma - hotspot prior\n")
			mrfG=""
			gammaPrior = "hotspot"
		}
		else
		{
			cat( "No value for gammaPrior was specified, but mrfG was given - choosing MRF prior\n")
			gammaPrior = "MRF"
		}
	}else{

    if ( gammaPrior == "hotspot" | gammaPrior == "HOTSPOT" | gammaPrior == "hotspots" | gammaPrior == "HOTSPOTS" | gammaPrior == "hs" | gammaPrior == "HS" )
      gammaPrior = "hotspot"
    else if ( gammaPrior == "MRF" | gammaPrior == "mrf" | gammaPrior == "markov random field" | gammaPrior == "Markov Random Field" ) 
      gammaPrior = "MRF"
    else if ( gammaPrior == "hierarchical" | gammaPrior == "h" | gammaPrior == "H" ) 
      gammaPrior = "hierarchical"
    else
      my_stop("Unknown gammaPrior argument: only hotspot, MRF or hierarchical are available",tmpFolder)
  }

  # if mrfG is not a string
  if( !(is.character(mrfG) & length(mrfG) == 1) )
  {
    # if it's a matrix
    if( is.numeric(mrfG) & !is.null(dim(mrfG))  )
    {
      write.table(mrfG,paste(sep="",tmpFolder,"mrfG.txt"), row.names = FALSE, col.names = FALSE)
      mrfG = paste(sep="",tmpFolder,"mrfG.txt")
    }else if( is.null( mrfG )){
        mrfG=""
    }else
      my_stop("Unknown mrfG argument: check the help function for possibile values",tmpFolder)
  }


  if ( betaPrior == "independent" | betaPrior == "Independent" | betaPrior == "INDEPENDENT" | betaPrior == "indep" | betaPrior == "Indep" | betaPrior == "i" | betaPrior == "I" )
    betaPrior = "independent"
  else if ( betaPrior == "gprior" | betaPrior == "gPrior" | betaPrior == "g-prior" | betaPrior == "G-Prior" | betaPrior == "GPRIOR" ) 
    betaPrior = "g-prior"
  else
    my_stop("Unknown betaPrior method: only independent is available as of yet",tmpFolder) # g prior is accepted but will return an error later

  
  ## Create the directory for the results
  dir.create(outFilePath)
  
  ## Create the return object
  ret = list( status=1, input=list(), output = list() )
  class(ret) = "R2SSUR"
  
  # Copy the inputs
  ret$input["nIter"] = nIter
  ret$input["burnin"] = burnin
  ret$input["nChains"] = nChains 
  ret$input["covariancePrior"] = covariancePrior
  ret$input["gammaPrior"] = gammaPrior
  ret$input["gammaSampler"] = gammaSampler
  ret$input["gammaInit"] = gammaInit
  ret$input["mrfG"] = mrfG
  ret$input["betaPrior"] = betaPrior

  methodString = 
    switch( covariancePrior,
      "HIW" = "SSUR" ,
      "IW"  = "dSUR" ,
      "IG"  = "HESS" )

  # Prepare path to outputs
  ret$output["logP"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_logP_out.txt")

  if ( output_gamma )
    ret$output["gamma"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_gamma_out.txt")
  
  if( gammaPrior %in% c("hierarchical","hotspot") & output_pi )
    ret$output["pi"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_pi_out.txt")

  if( gammaPrior == "hotspot" & output_tail )
    ret$output["tail"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_hotspot_tail_p_out.txt")
  
  if ( output_beta )
    ret$output["beta"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_beta_out.txt")
  
  if ( covariancePrior == "HIW" & output_G )
    ret$output["G"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_G_out.txt")
    
  if ( covariancePrior %in% c("HIW","IW") & output_sigmaRho )
    ret$output["sigmaRho"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_sigmaRho_out.txt")

  if ( output_beta )
    ret$output["model_size"] = paste(sep="", outFilePath , dataString , "_",  methodString , "_model_size_out.txt")
  
  
  ret$status = R2SSUR_internal(data, blockList, structureGraph, outFilePath, nIter, burnin, nChains, 
            covariancePrior, gammaPrior, gammaSampler, gammaInit, mrfG, betaPrior,
            output_gamma, output_beta, output_G, output_sigmaRho, output_pi, output_tail, output_model_size)

  if(outFilePath != tmpFolder)
    unlink(tmpFolder,recursive = TRUE)
  
  return(ret)
}
