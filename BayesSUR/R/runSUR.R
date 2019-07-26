#' BayesSUR -- Bayesian Seemingly Unrelated Regression
#' @title runSUR
#' @description
#' Run a SUR Bayesian sampler
#' @name runSUR
#' @param data either a matrix/dataframe or the path to (a plain text) data file with variables on the columns and observations on the rows 
#' @param Y,X,X_0 vectors of indexes (with respect to the data matrix) for the outcomes, the covariates to select and the fixed covariates respectively if data is either a path to a file or a matrix;
#' if the 'data' argument is not provided, these needs to be matrices containing the data instead.
#' @param outFilePath path to where the output files are to be written
#' @param nIter number of iterations for the MCMC procedure
#' @param burnin number of iterations (or fraction of iterations) to discard at the start of the chain Default = 0
#' @param nChains number of parallel chains to run
#' @param covariancePrior string indicating the prior for the covariance $C$; it has to be either "HIW" for the hyper-inverse-Wishar (which will result in a sparse covariance matrix),
#' "IW" for the inverse-Wishart prior ( dense covariance ) or "IG" for independent inverse-Gamma on all the diagonal elements and 0 otherwise.
#' @param gammaPrior string indicating the gamma prior to use, either "hotspot" for the Hotspot prior of Bottolo (2011), "MRF" for the Markov Random Field prior or "hierarchical" for a simpler hierarchical prior
#' @param gammaSampler string indicating the type of sampler for gamma, either "bandit" for the Thompson sampling inspired samper or "MC3" for the usual $MC^3$ sampler
#' @param gammaInit gamma initialisation to either all-zeros ("0"), all ones ("1"), randomly ("R") or (default) MLE-informed ("MLE").
#' @param mrfG either a matrix or a path to the file containing the G matrix for the MRF prior on gamma (if necessary)
#' @param standardize Logical flag for X variable standardization. Default is standardize=TRUE. The coefficients are returned on the standardized scale.
#' @param standardize.response Standardization for the response variables. Default is standardize.response=TRUE.
#' @param hyperpar a list of named hypeparameters to use instead of the default values; valid names are mrf_d, mrf_e, a_sigma, b_sigma, a_tau, b_tau, nu, a_eta, b_eta, a_o, b_o, a_pi, b_pi, a_w and b_w; see Documentation for more details
#' @param output_* allow ( TRUE ) or suppress ( FALSE ) the outut for *; possible outputs are gamma, G, beta, sigmaRho, pi, tail (hotspot tail probability) or model_size
#' @param tmpFolder the path to a temporary folder where intermediate data files are stored (will be erased at the end of the chain) default to local tmpFolder
#'
#' @examples
#' \donttest{
#' 
#' data(example_eQTL, package = "BayesSUR")
#' hyperpar = list( a_w = 2 , b_w = 5 )
#' 
#' fit = BayesSUR::runSUR(example_eQTL[["data"]],outFilePath = "results/",
#'                      Y = example_eQTL[["blockList"]][[1]],
#'                      X = example_eQTL[["blockList"]][[2]],
#'                      nIter = 100, nChains = 2, gammaPrior = "hotspot",
#'                      hyperpar = hyperpar, tmpFolder="tmp/" )
#' 
#' # show the interaction of plots
#' plot(fit)

#' ## check output
#' # show the estimated beta, gamma and graph of responeses Gy
#' plotEstimator(fit);
#' 
#' }
#' 
#' @export
runSUR = function(data=NULL, Y, X, X_0=NULL,outFilePath="", 
                nIter=10, burnin=0, nChains=1, 
                covariancePrior="HIW", gammaPrior="",
                gammaSampler="bandit", gammaInit="MLE", mrfG=NULL,
                standardize = TRUE, standardize.response=TRUE,
                output_gamma = TRUE, output_beta = TRUE, output_G = TRUE, output_sigmaRho = TRUE,
                output_pi = TRUE, output_tail = TRUE, output_model_size = TRUE, output_Y = TRUE, output_X = TRUE,
                hyperpar=list(),tmpFolder="tmp/")
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
  
  outFilePathLength = nchar(outFilePath)
  if( outFilePathLength > 0 )
  {
    if( substr(outFilePath,outFilePathLength,outFilePathLength) != "/" )
      outFilePath = paste( outFilePath , "/" , sep="" )
    if( substr(outFilePath,1,1) != "/" )
      outFilePath = paste( getwd(), "/", outFilePath , sep="" )
    dir.create(outFilePath)
  }

  # we'll check in reverse order, is data NULL?
  if ( is.null( data ) )
  {

    # Y,X (and if there X_0) need to be valid numeric matrices then
    # check Y and X have comfortable number of observations
    if( !is.numeric(Y) | is.null(dim(Y)) )
      my_stop("Y needs to be a valid matrix or data.frame with > 1 column ")
    
    nObservations = nrow(Y)
    if( !is.numeric(X) | is.null(dim(X)) | nrow(X) != nObservations )
      my_stop("X needs to be a valid matrix or data.frame with >= 1 column and the same number of rows of Y")
    
    if ( is.null ( X_0 ) ){
      X_0 = matrix(NA,nrow=nObservations,ncol=0)
    }else{
      if( !is.numeric(X_0) | is.null(dim(X_0)) | nrow(X_0) != nObservations )
        my_stop("if provided, X_0 needs to be a valid matrix or data.frame with >= 1 column and the same number of rows of Y")
    }
      
    #if( !is.numeric(X_0) | is.null(dim(X_0)) | nrow(X_0) != nObservations )
    #  my_stop("if provided, X_0 needs to be a valid matrix or data.frame with >= 1 column and the same number of rows of Y")

    # Standarize the data
    if( standardize ){
      X = scale(X)[,]
      X_0 = scale(X_0)[,]
    } 
    if( standardize.response ) Y = scale(Y)[,]

    # Write the three down in a single data file
    write.table(cbind(Y,X,X_0),paste(sep="",tmpFolder,"data.txt"), row.names = FALSE, col.names = FALSE)
    data = paste(sep="",tmpFolder,"data.txt")

    blockLabels = c( rep(0,ncol(Y)) , rep(1,ncol(X)) , rep(2,ncol(X_0)) )

    # Write the data in a output file
    write.table(Y,paste(sep="",outFilePath,"data_Y.txt"), row.names = FALSE, col.names = TRUE)
    write.table(X,paste(sep="",outFilePath,"data_X.txt"), row.names = FALSE, col.names = TRUE)
    write.table(X_0,paste(sep="",outFilePath,"data_X0.txt"), row.names = FALSE, col.names = TRUE)

  }else{ # data is not null, so the user wants to use it to input the data

    # is the data given as matrix?
    ## If it's valid matrix, simply write it and re-assign the variable data to hold its path
    if( is.numeric(data) & !is.null(dim(data))  )
    {
      # Standarize the data
      if( standardize ){
        data[,X]= scale(data[,X])[,]
        data[,X_0] = scale(data[,X_0])[,]
      } 
      if( standardize.response ) data[,Y] = scale(data[,Y])[,]
      
      # Write the Y and X data in a output file
      write.table(data[,Y],paste(sep="",outFilePath,"data_Y.txt"), row.names = FALSE, col.names = TRUE)
      write.table(data[,X],paste(sep="",outFilePath,"data_X.txt"), row.names = FALSE, col.names = TRUE)
      write.table(data[,X_0],paste(sep="",outFilePath,"data_X0.txt"), row.names = FALSE, col.names = TRUE)
      
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
    dataHeader = read.table(data,header = FALSE,nrow = 1)
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
  if ( length ( X_0 ) > 0 ){
    structureGraph = structureGraph = matrix(c(0,0,0,1,0,0,2,0,0),3,3,byrow=TRUE)
  }else structureGraph = structureGraph = matrix(c(0,0,1,0),2,2,byrow=TRUE)

  ## Finally write blockLabels and structureGraph to a file
  write.table(blockLabels,paste(sep="",tmpFolder,"blockLabels.txt"), row.names = FALSE, col.names = FALSE)
  blockList = paste(sep="",tmpFolder,"blockLabels.txt")
  write.table(structureGraph, paste(sep="",tmpFolder,"structureGraph.txt"), row.names = FALSE, col.names = FALSE)
  structureGraph = paste(sep="",tmpFolder,"structureGraph.txt")


  # cleanup file PATHS
  dataLength = nchar(data)
  if( dataLength == 0 )
    my_stop("Please provide a correct path to a plain-text (.txt) file", tmpFolder)
  
  # magicly strip '/' from the start and '.txt' from the end of the data file name
  dataString = head( strsplit( tail( strsplit(data,split = c("/"))[[1]] , 1 ) , ".txt" )[[1]] , 1 ) 


  # check how burnin was given
  if ( burnin < 0 ){
    my_stop("Burnin must be positive or 0",tmpFolder)
  }else{ if ( burnin > nIter ){
    my_stop("Burnin might not be greater than nIter",tmpFolder)
  }else{ if ( burnin < 1 ){ # given as a fraction
    burnin = ceiling(nIter * burnin) # the zero case is taken into account here as well
  }}} # else assume is given as an absolute number

  ###############################

  # method to use
  if ( toupper(covariancePrior) %in% c("SPARSE", "HIW") ){
    covariancePrior = "HIW"
  }else if ( toupper(covariancePrior) %in% c("DENSE", "IW") ){
    covariancePrior = "IW"
  }else if ( toupper(covariancePrior) %in% c("INDEPENDENT", "INDEP", "IG") ){
    covariancePrior = "IG"
  }else
    my_stop("Unknown covariancePrior argument: only sparse (HIW), dense(IW) or independent (IG) are available",tmpFolder)
  
  # mrfG and gammaPrior
  if( gammaPrior == "" )
	{
		if ( is.null(mrfG) )
		{
			cat( "Using default prior for Gamma - hotspot prior\n")
			#mrfG=""
			gammaPrior = "hotspot"
		}
		else
		{
			cat( "No value for gammaPrior was specified, but mrfG was given - choosing MRF prior\n")
			gammaPrior = "MRF"
		}
    
	}else{

    if ( toupper(gammaPrior) %in% c("HOTSPOT", "HOTSPOTS", "HS") )
      gammaPrior = "hotspot"
    else if ( toupper(gammaPrior) %in% c("MRF", "MARKOV RANDOM FIELD") ) 
      gammaPrior = "MRF"
    else if ( toupper(gammaPrior) %in% c("HIERARCHICAL", "H") ) 
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
      # save a meaningless mrfG.txt file to pass the parameter to C++ 
      mrfG = matrix(c(0,0),ncol=2)
      write.table(mrfG, paste(sep="",tmpFolder,"mrfG.txt"), row.names = FALSE, col.names = FALSE)
      mrfG = paste(sep="",tmpFolder,"mrfG.txt")
    }else
      my_stop("Unknown mrfG argument: check the help function for possibile values",tmpFolder)
  }
  
  ## Set up the XML file for hyperparameters
  if("xml2" %in% rownames(installed.packages()) == FALSE)
    my_stop("Using non-default hyperparameters require an XML file to be written, please install the \"xml2\" package.",tmpFolder)
    
  xml  = xml2::as_xml_document(
    list( hyperparameters = list(
      lapply(hyperpar,function(x) list(x)) # every element in the list should be a list
    )))
  hyperParFile = paste(sep="",tmpFolder,"hyperpar.xml")
  xml2::write_xml(xml,file = hyperParFile)
  
  ## Create the directory for the results
  dir.create(outFilePath)
  
  ## Create the return object
  ret = list( status=1, input=list(), output = list() )
  class(ret) = "BayesSUR"
  
  # Copy the inputs
  ret$input["nIter"] = nIter
  ret$input["burnin"] = burnin
  ret$input["nChains"] = nChains 
  ret$input["covariancePrior"] = covariancePrior
  ret$input["gammaPrior"] = gammaPrior
  ret$input["gammaSampler"] = gammaSampler
  ret$input["gammaInit"] = gammaInit
  ret$input["mrfG"] = mrfG
  
  ret$input$hyperParameters = hyperpar

  methodString = 
    switch( covariancePrior,
      "HIW" = "SSUR" ,
      "IW"  = "dSUR" ,
      "IG"  = "HRR" )

  # Prepare path to outputs
  ret$output["outFilePath"] = outFilePath
  
  ret$output["logP"] = paste(sep="", dataString , "_",  methodString , "_logP_out.txt")

  if ( output_gamma )
    ret$output["gamma"] = paste(sep="", dataString , "_",  methodString , "_gamma_out.txt")
  
  if( gammaPrior %in% c("hierarchical","hotspot") & output_pi )
    ret$output["pi"] = paste(sep="", dataString , "_",  methodString , "_pi_out.txt")

  if( gammaPrior == "hotspot" & output_tail )
    ret$output["tail"] = paste(sep="", dataString , "_",  methodString , "_hotspot_tail_p_out.txt")
  
  if ( output_beta )
    ret$output["beta"] = paste(sep="", dataString , "_",  methodString , "_beta_out.txt")
  
  if ( covariancePrior == "HIW" & output_G )
    ret$output["G"] = paste(sep="", dataString , "_",  methodString , "_G_out.txt")
    
  if ( covariancePrior %in% c("HIW","IW") & output_sigmaRho )
    ret$output["sigmaRho"] = paste(sep="", dataString , "_",  methodString , "_sigmaRho_out.txt")

  if ( output_beta )
    ret$output["model_size"] = paste(sep="", dataString , "_",  methodString , "_model_size_out.txt")
  
  if ( output_Y )
    ret$output["Y"] = paste(sep="", "data_Y.txt")
  
  if ( output_X )
    ret$output["X"] = paste(sep="", "data_X.txt")

  betaPrior="independent"
  ret$status = BayesSUR_internal(data, mrfG, blockList, structureGraph, hyperParFile, outFilePath, 
            nIter, burnin, nChains, 
            covariancePrior, gammaPrior, gammaSampler, gammaInit, betaPrior,
            output_gamma, output_beta, output_G, output_sigmaRho, output_pi, output_tail, output_model_size)

  if(outFilePath != tmpFolder)
    unlink(tmpFolder,recursive = TRUE)
   
  class(summary) <- c(class(summary), "BayesSUR")
  class(plot) <- c(class(plot), "BayesSUR")
  return(ret)
}


my_stop = function( msg , tmpFolder )
{
  unlink(tmpFolder,recursive = TRUE)
  stop(msg)
}

