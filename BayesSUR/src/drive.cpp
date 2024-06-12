#include "drive.h"

#ifndef CCODE
using Rcpp::Rcout;
using Rcpp::Rcerr;
#else
#define Rcout std::cout
#define Rcerr std::cerr
#endif

#ifdef _OPENMP
extern omp_lock_t RNGlock; /*defined in global.h*/
#endif
/*extern std::vector<std::mt19937_64> rng;*/

#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]


using Utils::Chain_Data;

int drive_SUR( Chain_Data& chainData )
{
    
    // ****************************************
    // **********  INIT THE CHAIN *************
    // ****************************************
    Rcout << "Initialising the (SUR) MCMC Chain";
    
    ESS_Sampler<SUR_Chain> sampler( chainData.surData , chainData.nChains , 1.2 ,
                                   chainData.gamma_sampler_type, chainData.gamma_type, chainData.beta_type, chainData.covariance_type, 
                                   chainData.output_CPO, chainData.maxThreads, chainData.tick, chainData.burnin );
    
    Rcout << " ... ";
    
    // *****************************
    
    // extra step needed to read in the MRF prior G matrix if needed
    // ************************************
    /*    if ( chainData.gamma_type == Gamma_Type::mrf )
     {
     for( unsigned int i=0; i< chainData.nChains; ++i )
     {
     sampler[i]->mrfGInit( chainData.mrfG );
     sampler[i]->logPGamma();
     }
     }
     */
    // Init all parameters
    // *****************************
    sampler.setHyperParameters( chainData );
    
    Rcout << " ... ";
    
    // *****************************
    
    // set when the JT move should start
    unsigned int jtStartIteration = 0;
    if( chainData.covariance_type == Covariance_Type::HIW )
    {
        jtStartIteration = chainData.nIter/10;
        for( unsigned int i=0; i< chainData.nChains; ++i )
            sampler[i]->setJTStartIteration( jtStartIteration );
    }
    
    // Init gamma and beta for the main chain
    // *****************************
    sampler[0] -> gammaInit( chainData.gammaInit );
    if(!chainData.betaInit.is_empty())
        sampler[0] -> betaInit( chainData.betaInit );
    sampler[0] -> updateQuantities();
    sampler[0] -> logLikelihood();
    sampler[0] -> stepSigmaRhoAndBeta();
    if ( chainData.output_CPO )
        sampler[0] -> predLikelihood();
    
    
    // ****************************************
    Rcout << " DONE!\nDrafting the output files with the start of the chain ... ";
    
    // INIT THE FILE OUTPUT
    std::string outFilePrefix = chainData.outFilePath+chainData.filePrefix;
    
    // clear the content of previous files
    std::ofstream logPOutFile; logPOutFile.open(outFilePrefix+"logP_out.txt", std::ios::out | std::ios::trunc); logPOutFile.close();
    // openlogP file in append mode
    logPOutFile.open( outFilePrefix+"logP_out.txt" , std::ios_base::app); // note we don't close!
    // open avg files in trunc mode to cut previous content
    
    std::ofstream gammaOutFile;
    if ( chainData.output_gamma )
    {
        gammaOutFile.open( outFilePrefix+"gamma_out.txt" , std::ios_base::trunc);
        gammaOutFile.close();
    }
    
    std::ofstream gOutFile;
    if ( chainData.covariance_type == Covariance_Type::HIW && chainData.output_Gy )
    {
        gOutFile.open( outFilePrefix+"Gy_out.txt" , std::ios_base::trunc);
        gOutFile.close();
    }
    
    std::ofstream piOutFile;
    if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) && chainData.output_pi )
    {
        piOutFile.open( outFilePrefix+"pi_out.txt" , std::ios_base::trunc);
        piOutFile.close();
    }
    
    std::ofstream htpOutFile;
    if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
    {
        htpOutFile.open( outFilePrefix+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
        htpOutFile.close();
    }
    
    std::ofstream ModelSizeOutFile;
    if ( chainData.output_model_size )
    {
        ModelSizeOutFile.open( outFilePrefix+"model_size_out.txt", std::ios::out | std::ios::trunc); ModelSizeOutFile.close();
        ModelSizeOutFile.open( outFilePrefix+"model_size_out.txt" , std::ios_base::app); // note we don't close!
    }
    
    std::ofstream GVisitOutFile;
    std::ofstream ModelVisitGammaOutFile;
    std::ofstream ModelVisitGOutFile;
    if( chainData.output_model_visit )
    {
        GVisitOutFile.open( outFilePrefix+"Gy_visit.txt", std::ios::out | std::ios::trunc); GVisitOutFile.close();
        GVisitOutFile.open( outFilePrefix+"Gy_visit.txt" , std::ios_base::app); // note we don't close!
        
        ModelVisitGammaOutFile.open( outFilePrefix+"model_visit_gamma_out.txt", std::ios::out | std::ios::trunc); ModelVisitGammaOutFile.close();
        ModelVisitGammaOutFile.open( outFilePrefix+"model_visit_gamma_out.txt" , std::ios_base::app); // note we don't close!
        
        ModelVisitGOutFile.open( outFilePrefix+"model_visit_gy_out.txt", std::ios::out | std::ios::trunc); ModelVisitGOutFile.close();
        ModelVisitGOutFile.open( outFilePrefix+"model_visit_gy_out.txt" , std::ios_base::app); // note we don't close!
    }
    
    // Output to file the initial state (if burnin=0)
    arma::umat gamma_out; // out var for the gammas
    arma::umat g_out, tmpG; // out var for G and tmpG
    arma::urowvec g_visit; // out var for visted G (vectorized upper triangle)
    arma::mat beta_out, tmpB; // out var for the betas and standard deviation
    arma::mat betaSD_out = arma::zeros<arma::mat>(chainData.surData.nFixedPredictors+chainData.surData.nVSPredictors,chainData.surData.nOutcomes); // out var for the betas SD
    arma::mat sigmaRho_out; // out var for the sigmas and rhos
    arma::vec tmpVec; // temporary to store the pi parameter vector
    arma::vec pi_out;
    arma::vec hotspot_tail_prob_out;
    arma::mat cpo_out, predLik;
    arma::vec cposumy_out; // CPO with each element summerizing all response variables
    arma::mat lpd, waic_out, waic_frac_sum;

    if( chainData.burnin == 0 )
    {
        if ( chainData.output_gamma )
        {
            gamma_out = sampler[0] -> getGamma();
            gammaOutFile.open( outFilePrefix+"gamma_out.txt" , std::ios_base::trunc);
            gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out));
            gammaOutFile.close();
        }
        
        if ( chainData.covariance_type == Covariance_Type::HIW && chainData.output_Gy )
        {
            tmpG = arma::umat( sampler[0] -> getGAdjMat() );
            g_out = tmpG;
            gOutFile.open( outFilePrefix+"Gy_out.txt" , std::ios_base::trunc);
            gOutFile << ( arma::conv_to<arma::mat>::from(g_out) );   // this might be quite long...
            gOutFile.close();
        }
        
        if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) &&
            ( chainData.output_pi || chainData.output_tail ) )
        {
            tmpVec = sampler[0] -> getPi();
            if ( chainData.output_pi )
            {
                pi_out = tmpVec;
                piOutFile.open( outFilePrefix+"pi_out.txt" , std::ios_base::trunc);
                piOutFile << pi_out;
                piOutFile.close();
            }
            
            if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
            {
                tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
                hotspot_tail_prob_out = tmpVec;
                
                htpOutFile.open( outFilePrefix+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
                htpOutFile << hotspot_tail_prob_out;
                htpOutFile.close();
            }
        }
        
        if ( chainData.output_beta ){
            tmpB = sampler[0] -> getBeta();
            beta_out = tmpB;
            betaSD_out = arma::square( tmpB );
        }
        
        if ( chainData.output_sigmaRho )
            sigmaRho_out = sampler[0] -> getSigmaRho();
        
    }else{
        if ( chainData.output_gamma )
            gamma_out = sampler[0] -> getGamma();
        if ( chainData.covariance_type == Covariance_Type::HIW && chainData.output_Gy )
        {
            tmpG = arma::umat( sampler[0] -> getGAdjMat() );
            g_out = tmpG;
        }
            
        if ( chainData.output_beta ){
            tmpB = sampler[0] -> getBeta();
            beta_out = tmpB;
            betaSD_out = arma::square( tmpB );
        }
        if ( chainData.output_sigmaRho )
            sigmaRho_out = sampler[0] -> getSigmaRho();
        if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) &&
            ( chainData.output_pi || chainData.output_tail ) )
        {
            tmpVec = sampler[0] -> getPi();
            if ( chainData.output_pi )
                pi_out = tmpVec;
            if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
            {
                tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
                hotspot_tail_prob_out = tmpVec;
            }
        }
    }
    
    logPOutFile <<     sampler[0] -> getLogPTau() << " ";
    logPOutFile <<     sampler[0] -> getLogPEta() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPJT() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPSigmaRho() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPO() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPPi() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPGamma() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPW() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPBeta() <<  " ";
    logPOutFile <<     sampler[0] -> getLogLikelihood();
    logPOutFile <<     '\n';
    
    if ( chainData.output_model_size )
    {
        ModelSizeOutFile << sampler[0]->getModelSize() << " ";
        ModelSizeOutFile << '\n';
    }
    
    if ( chainData.output_model_visit )
    {
        if ( chainData.covariance_type == Covariance_Type::HIW && chainData.output_Gy )
        {
            g_visit.clear();
            for(unsigned int k=0; k < tmpG.n_cols-1; ++k)
            {
                g_visit = join_rows( g_visit, tmpG.submat(k,k+1, k,tmpG.n_cols-1) );
            }
            GVisitOutFile << g_visit << " ";
            GVisitOutFile << '\n';
        }
        
        ModelVisitGammaOutFile << arma::find((sampler[0] -> getGamma()) == 1).t() << " ";
        ModelVisitGammaOutFile << '\n';
        
        ModelVisitGOutFile << arma::find(arma::conv_to<arma::umat>::from(sampler[0] -> getGAdjMat()) == 1).t() << " ";
        ModelVisitGOutFile << '\n';
    }

    if ( chainData.output_CPO )
    {
        predLik =  sampler[0] -> getPredLikelihood() ;
        cpo_out = 1./predLik;
        cposumy_out = 1./arma::exp( arma::sum( arma::log(predLik), 1 ) );
        
        lpd = predLik;
        waic_out = arma::square( arma::log(predLik) );
        waic_frac_sum = arma::log(predLik);
    }
    
    
    // ########
    // ########
    // ######## Start
    // ########
    // ########
    
    Rcout << "DONE! \n\nStarting "<< chainData.nChains <<" (parallel) chain(s) for " << chainData.nIter << " iterations:" << '\n';
    
    //unsigned int tick = 1000; // how many iter for each print?
    
    for(unsigned int i=1; i < chainData.nIter ; ++i)
    {
        sampler.step();
       
        // #################### END LOCAL MOVES
        
        // ## Global moves
        // *** end Global move's section
        
        // UPDATE OUTPUT STATE
        if( i >= chainData.burnin )
        {
            if ( chainData.output_gamma )
                gamma_out += sampler[0] -> getGamma(); // the result of the whole procedure is now my new mcmc point, so add that up
            
            if ( chainData.covariance_type == Covariance_Type::HIW && chainData.output_Gy )
            {
                tmpG = arma::umat( sampler[0] -> getGAdjMat() );
                g_out += tmpG;
            }
            
            if ( chainData.output_beta ){
                tmpB = sampler[0] -> getBeta();
                beta_out += tmpB;
                betaSD_out += arma::square( tmpB );
            }
            
            if ( chainData.output_sigmaRho )
                sigmaRho_out += sampler[0] -> getSigmaRho();
            
            if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) &&
                ( chainData.output_pi || chainData.output_tail ) )
            {
                tmpVec = sampler[0] -> getPi();
                if ( chainData.output_pi )
                    pi_out += tmpVec;
                
                if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
                {
                    tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
                    hotspot_tail_prob_out += tmpVec;
                }
            }
            
            if( chainData.output_CPO )
            {
                predLik = sampler[0] -> getPredLikelihood();
                cpo_out += 1./predLik;
                cposumy_out += 1./arma::exp( arma::sum( arma::log(predLik), 1 ) );
                
                lpd += predLik;
                waic_out += arma::square( arma::log(predLik) );
                waic_frac_sum += arma::log(predLik);
            }
            
            // Nothing to update for model size
        }else{
            if ( chainData.covariance_type == Covariance_Type::HIW && chainData.output_Gy )
                tmpG = arma::umat( sampler[0] -> getGAdjMat() );
        }
        
        if ( chainData.output_model_visit )
        {
            ModelVisitGammaOutFile << arma::find((sampler[0] -> getGamma()) == 1).t() << " ";
            ModelVisitGammaOutFile << '\n';
            
            ModelVisitGOutFile << arma::find(arma::conv_to<arma::umat>::from(sampler[0] -> getGAdjMat()) == 1).t() << " ";
            ModelVisitGOutFile << '\n';
        }
        
        // Print something on how the chain is going
        if( (i+1) % chainData.tick == 0 )
        {
            Rcout << " Running iteration " << i+1 << " ... local Acc Rate: ~ gamma: " << Utils::round( sampler[0] -> getGammaAccRate() , 3 );
            Rcout << " -- JT: " << Utils::round( sampler[0] -> getJTAccRate() , 3 ) ;
            
            if( chainData.nChains > 1){
                Rcout << " -- Global: " << Utils::round( sampler.getGlobalAccRate() , 3 ) << '\n';
            }else{
                Rcout << '\n';
            }
            
#ifndef CCODE
            Rcpp::checkUserInterrupt(); // this checks for interrupts from R
#endif
            
            // Output to files every now and then
            
            if( (i >= chainData.burnin) && ( (i-chainData.burnin+1) % (chainData.tick*1) == 0 ) )
            {
                
                if ( chainData.output_gamma )
                {
                    gammaOutFile.open( outFilePrefix+"gamma_out.txt" , std::ios_base::trunc);
                    gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/(double)(i+1.0-chainData.burnin);
                    gammaOutFile.close();
                }
                
                if ( chainData.covariance_type == Covariance_Type::HIW && chainData.output_Gy )
                {
                    gOutFile.open( outFilePrefix+"Gy_out.txt" , std::ios_base::trunc);
                    gOutFile << ( arma::conv_to<arma::mat>::from(g_out) )/((double)(i-std::max(jtStartIteration,chainData.burnin))+1.0);   // this might be quite long...
                    gOutFile.close();
                }
                
                if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) &&
                    chainData.output_pi )
                {
                    piOutFile.open( outFilePrefix+"pi_out.txt" , std::ios_base::trunc);
                    piOutFile << pi_out/(double)(i+1.0-chainData.burnin);
                    piOutFile.close();
                }
                
                if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
                {
                    htpOutFile.open( outFilePrefix+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
                    htpOutFile << hotspot_tail_prob_out/(double)(i+1.0-chainData.burnin);
                    htpOutFile.close();
                }
            }
            
            //if( (i-chainData.burnin+1) % (tick*1) == 0 )
            if( (i+1) % (chainData.tick*1) == 0 )
            {
                logPOutFile <<     sampler[0] -> getLogPTau() << " ";
                logPOutFile <<     sampler[0] -> getLogPEta() <<  " ";
                logPOutFile <<     sampler[0] -> getLogPJT() <<  " ";
                logPOutFile <<     sampler[0] -> getLogPSigmaRho() <<  " ";
                logPOutFile <<     sampler[0] -> getLogPO() <<  " ";
                logPOutFile <<     sampler[0] -> getLogPPi() <<  " ";
                logPOutFile <<     sampler[0] -> getLogPGamma() <<  " ";
                logPOutFile <<     sampler[0] -> getLogPW() <<  " ";
                logPOutFile <<     sampler[0] -> getLogPBeta() <<  " ";
                logPOutFile <<     sampler[0] -> getLogLikelihood();
                logPOutFile <<     '\n';
                
                if ( chainData.output_model_size )
                {
                    if ( chainData.covariance_type == Covariance_Type::HIW && chainData.output_Gy )
                    {
                        //g_visit = arma::conv_to<arma::urowvec>::from( arma::trimatu(tmpG, 1) );
                        g_visit.clear();
                        for(unsigned int k=0; k < tmpG.n_cols-1; ++k)
                        {
                            g_visit = join_rows( g_visit, tmpG.submat(k,k+1, k,tmpG.n_cols-1) );
                        }
                        GVisitOutFile << g_visit << " ";
                        GVisitOutFile << '\n';
                    }
                    
                    ModelSizeOutFile << sampler[0]->getModelSize() << " ";
                    ModelSizeOutFile << '\n';
                }
            }
        }
        
    } // end MCMC
    
    // Print the end
    Rcout << " MCMC ends. " /* << " Final temperature ratio ~ " << temperatureRatio  */<< "  --- Saving results and exiting" << '\n';
    
    // ### Collect results and save them
    if ( chainData.output_gamma )
    {
        gammaOutFile.open( outFilePrefix+"gamma_out.txt" , std::ios_base::trunc);
        gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/(double)(chainData.nIter-chainData.burnin+1.);
        gammaOutFile.close();
    }
    
    if ( chainData.covariance_type == Covariance_Type::HIW && chainData.output_Gy )
    {
        gOutFile.open( outFilePrefix+"Gy_out.txt" , std::ios_base::trunc);
        gOutFile << ( arma::conv_to<arma::mat>::from(g_out) )/(double)(chainData.nIter-std::max(jtStartIteration,chainData.burnin)+1.);   // this might be quite long...
        gOutFile.close();
    }
    
    logPOutFile <<     sampler[0] -> getLogPTau() << " ";
    logPOutFile <<     sampler[0] -> getLogPEta() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPJT() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPSigmaRho() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPO() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPPi() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPGamma() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPW() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPBeta() <<  " ";
    logPOutFile <<     sampler[0] -> getLogLikelihood();
    logPOutFile <<     '\n';
    logPOutFile.close();
    
    if ( chainData.output_model_size )
    {
        ModelSizeOutFile << sampler[0]->getModelSize();
        ModelSizeOutFile << '\n';
        ModelSizeOutFile.close();
    }
    
    GVisitOutFile.close();
    
    if ( chainData.output_model_visit )
    {
        ModelVisitGammaOutFile.close();
        ModelVisitGOutFile.close();
        
    }
    
    // ----
    if ( chainData.output_beta )
    {
        beta_out = beta_out/(double)(chainData.nIter-chainData.burnin+1);
        beta_out.save(outFilePrefix+"beta_out.txt",arma::raw_ascii);
        
        betaSD_out = arma::sqrt( betaSD_out/(double)(chainData.nIter-chainData.burnin+1) - arma::square(beta_out) );
        betaSD_out.save(outFilePrefix+"betaSD_out.txt",arma::raw_ascii);
    }
    
    if ( chainData.output_sigmaRho )
    {
        sigmaRho_out = sigmaRho_out/(double)(chainData.nIter-chainData.burnin+1);
        sigmaRho_out.save(outFilePrefix+"sigmaRho_out.txt",arma::raw_ascii);
    }
    
    if ( chainData.output_CPO )
    {
        cpo_out = 1./( cpo_out/(double)(chainData.nIter-chainData.burnin+1) );
        cpo_out.save(outFilePrefix+"CPO_out.txt",arma::raw_ascii);
        
        cposumy_out = 1./( cposumy_out/(double)(chainData.nIter-chainData.burnin+1) );
        cposumy_out.save(outFilePrefix+"CPOsumy_out.txt",arma::raw_ascii);
        
        waic_out = arma::log( lpd/(double)(chainData.nIter-chainData.burnin+1) ) - ( waic_out/(double)(chainData.nIter-chainData.burnin+1) - arma::square(waic_frac_sum/(double)(chainData.nIter-chainData.burnin+1)) );
        waic_out.save(outFilePrefix+"WAIC_out.txt",arma::raw_ascii);
    }
    // -----
    
    // -----
    if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) && chainData.output_pi )
    {
        piOutFile.open( outFilePrefix+"pi_out.txt" , std::ios_base::trunc);
        piOutFile << pi_out/(double)(chainData.nIter-chainData.burnin+1);
        piOutFile.close();
    }
    
    if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
    {
        htpOutFile.open( outFilePrefix+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
        htpOutFile << hotspot_tail_prob_out/(double)(chainData.nIter-chainData.burnin+1);
        htpOutFile.close();
    }
    // -----
    
    Rcout << "Saved to :   "+outFilePrefix+"****_out.txt" << '\n';
    if( chainData.surData.nFixedPredictors > 0 )
        Rcout << "Final w0  : " << sampler[0] -> getW0() <<  '\n';
    Rcout << "Final w   : " << sampler[0] -> getW() <<  '\n';
    Rcout << "Final tau : " << sampler[0] -> getTau() << "    w/ proposal variance: " << sampler[0] -> getVarTauProposal() << '\n';
    if ( chainData.covariance_type == Covariance_Type::HIW )
        Rcout << "Final eta : " << sampler[0] -> getEta() <<  '\n';
    Rcout << "  -- Average Omega : " << arma::accu( sampler[0] -> getO() * sampler[0] -> getPi().t() )/((double)(sampler[0]->getP()*sampler[0]->getS())) <<  '\n';
    if( chainData.nChains > 1 )
        Rcout << "Final temperature ratio : " << sampler[1]->getTemperature() <<  '\n' << '\n' ;
    
    
    // Exit
    Rcout << "DONE, exiting! " << '\n' << '\n' ;
    return 0;
}

int drive_HRR( Chain_Data& chainData )
{
    
    // ****************************************
    // **********  INIT THE CHAIN *************
    // ****************************************
    Rcout << "Initialising the (HRR) MCMC Chain ";
    
    ESS_Sampler<HRR_Chain> sampler( chainData.surData , chainData.nChains , 1.2 ,
                                   chainData.gamma_sampler_type, chainData.gamma_type, chainData.beta_type, chainData.covariance_type, 
                                   chainData.output_CPO, chainData.maxThreads, chainData.tick, chainData.burnin );
    
    Rcout << " ... ";
    // *****************************
    
    // extra step needed to read in the MRF prior G matrix if needed
    // **********************************
    /*    if ( chainData.gamma_type == Gamma_Type::mrf )
     {
     for( unsigned int i=0; i< chainData.nChains; ++i)
     {
     sampler[i]->mrfGInit( chainData.mrfG );
     sampler[i]->logPGamma();
     }
     }
     */
    // Init all parameters
    // *****************************
    sampler.setHyperParameters( chainData );
    Rcout << " ... ";
    
    // Init gamma for the main chain
    // *****************************
    
    sampler[0] -> gammaInit( chainData.gammaInit ); // this updates gammaMask as well
    sampler[0] -> logLikelihood();
    if ( chainData.output_CPO )
        sampler[0] -> predLikelihood();
    
    // ****************************************
    
    Rcout << " DONE!\nDrafting the output files with the start of the chain ... ";
    
    // INIT THE FILE OUTPUT
    
    std::string outFilePrefix = chainData.outFilePath+chainData.filePrefix;
    
    // clear the content of previous files
    std::ofstream logPOutFile; logPOutFile.open(outFilePrefix+"logP_out.txt", std::ios::out | std::ios::trunc); logPOutFile.close();
    // openlogP file in append mode
    logPOutFile.open( outFilePrefix+"logP_out.txt" , std::ios_base::app); // note we don't close!
    // open avg files in trunc mode to cut previous content
    
    std::ofstream gammaOutFile;
    if ( chainData.output_gamma )
    {
        gammaOutFile.open( outFilePrefix+"gamma_out.txt" , std::ios_base::trunc);
        gammaOutFile.close();
    }
    
    std::ofstream piOutFile;
    if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) && chainData.output_pi )
    {
        piOutFile.open( outFilePrefix+"pi_out.txt" , std::ios_base::trunc);
        piOutFile.close();
    }
    
    std::ofstream htpOutFile;
    if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
    {
        htpOutFile.open( outFilePrefix+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
        htpOutFile.close();
    }
    
    logPOutFile <<     sampler[0] -> getLogPO() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPPi() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPGamma() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPW() <<  " ";
    logPOutFile <<     sampler[0] -> getLogLikelihood();
    logPOutFile <<     '\n';
    
    std::ofstream ModelSizeOutFile;
    if ( chainData.output_model_size )
    {
        ModelSizeOutFile.open(outFilePrefix+"model_size_out.txt", std::ios::out | std::ios::trunc); ModelSizeOutFile.close(); // clear out previous content
        ModelSizeOutFile.open( outFilePrefix+"model_size_out.txt" , std::ios_base::app); //note we don't close it
    }
    
    std::ofstream ModelVisitGammaOutFile;
    if ( chainData.output_model_visit )
    {
        ModelVisitGammaOutFile.open( outFilePrefix+"model_visit_gamma_out.txt", std::ios::out | std::ios::trunc); ModelVisitGammaOutFile.close();
        ModelVisitGammaOutFile.open( outFilePrefix+"model_visit_gamma_out.txt" , std::ios_base::app); // note we don't close!
    }
    
    // Output to file the initial state (if burnin=0)
    arma::umat gamma_out; // out var for the gammas
    arma::mat beta_out; // out var for the betas
    
    arma::vec tmpVec; // temporary to store the pi parameter vector
    arma::vec pi_out;
    arma::vec hotspot_tail_prob_out;
    
    arma::mat cpo_out, predLik;
    arma::vec cposumy_out;
    //arma::mat lpd;
    arma::mat waic_out, waic_frac_sum;
    
    if( chainData.burnin == 0 )
    {
        if ( chainData.output_gamma )
        {
            gamma_out = sampler[0] -> getGamma();
            gammaOutFile.open( outFilePrefix+"gamma_out.txt" , std::ios_base::trunc);
            gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out));
            gammaOutFile.close();
        }
        
        if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) &&
            ( chainData.output_pi || chainData.output_tail ) )
        {
            tmpVec = sampler[0] -> getPi();
            if ( chainData.output_pi )
            {
                pi_out = tmpVec;
                piOutFile.open( outFilePrefix+"pi_out.txt" , std::ios_base::trunc);
                piOutFile << pi_out;
                piOutFile.close();
            }
            
            if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
            {
                tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
                hotspot_tail_prob_out = tmpVec;
                
                htpOutFile.open( outFilePrefix+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
                htpOutFile << hotspot_tail_prob_out;
                htpOutFile.close();
            }
        }
        
    }else{
        if ( chainData.output_gamma )
            gamma_out = sampler[0] -> getGamma();
        
        if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) &&
            ( chainData.output_pi || chainData.output_tail ) )
        {
            tmpVec = sampler[0] -> getPi();
            if ( chainData.output_pi )
                pi_out = tmpVec;
            if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
            {
                tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
                hotspot_tail_prob_out = tmpVec;
            }
        }
    }
    
    if ( chainData.output_beta )
        beta_out = sampler[0] -> getBeta();
    
    if ( chainData.output_CPO )
    {
        predLik =  sampler[0] -> getPredLikelihood() ;
        cpo_out = predLik;
        cposumy_out = arma::exp( arma::sum( arma::log(predLik), 1 ) );
        
        //lpd = predLik;
        waic_out = arma::square( arma::log(predLik) );
        waic_frac_sum = arma::log(predLik);
    }
    
    logPOutFile <<     sampler[0] -> getLogPO() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPPi() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPGamma() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPW() <<  " ";
    logPOutFile <<     sampler[0] -> getLogLikelihood();
    logPOutFile <<     '\n';
    
    if ( chainData.output_model_size )
    {
        ModelSizeOutFile << sampler[0]->getModelSize() << " ";
        ModelSizeOutFile << '\n';
    }
    
    if ( chainData.output_model_visit )
    {
        ModelVisitGammaOutFile << arma::find((sampler[0] -> getGamma()) == 1).t() << " ";
        ModelVisitGammaOutFile << '\n';
    }
    
    // ########
    // ########
    // ######## Start
    // ########
    // ########
    
    Rcout << "DONE! \n\nStarting "<< chainData.nChains <<" (parallel) chain(s) for " << chainData.nIter << " iterations:" << '\n';
    
    //unsigned int tick = 1000; // how many iter for each print?
    
    for(unsigned int i=1; i < chainData.nIter ; ++i)
    {
        
        sampler.step();
        
        // #################### END LOCAL MOVES
        
        // ## Global moves
        // *** end Global move's section
        
        
        // UPDATE OUTPUT STATE
        if( i >= chainData.burnin )
        {
            if ( chainData.output_gamma )
                gamma_out += sampler[0] -> getGamma();; // the result of the whole procedure is now my new mcmc point, so add that up
            
            if ( chainData.output_beta )
                beta_out += sampler[0] -> getBeta();
            
            if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) &&
                ( chainData.output_pi || chainData.output_tail ) )
            {
                tmpVec = sampler[0] -> getPi();
                if ( chainData.output_pi )
                    pi_out += tmpVec;
                
                if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
                {
                    tmpVec.for_each( [](arma::vec::elem_type& val) { if(val>1.0) val = 1.0; else val=0.0; } );
                    hotspot_tail_prob_out += tmpVec;
                }
            }
            
            if( chainData.output_CPO )
            {
                predLik = sampler[0] -> getPredLikelihood();
                cpo_out += predLik;
                cposumy_out += arma::exp( arma::sum( arma::log(predLik), 1 ) );
                
                //lpd += predLik;
                waic_out += arma::square( arma::log(predLik) );
                waic_frac_sum += arma::log(predLik);
            }
            
            // Nothing to update for model size
        }
        
        if( chainData.output_model_visit )
        {
            ModelVisitGammaOutFile << arma::find((sampler[0] -> getGamma()) == 1).t() << " ";
            ModelVisitGammaOutFile << '\n';
        }
        
        // Print something on how the chain is going
        if( (i+1) % chainData.tick == 0 )
        {
            
            Rcout << " Running iteration " << i+1 << " ... local Acc Rate: ~ gamma: " << Utils::round( sampler[0] -> getGammaAccRate() , 3 );
            
            if( chainData.nChains > 1)
                Rcout << " -- Global: " << Utils::round( sampler.getGlobalAccRate() , 3 ) << '\n';
            else
                Rcout << '\n';
            
            
#ifndef CCODE
            Rcpp::checkUserInterrupt(); // this checks for interrupts from R ... or does it?
#endif
            
            // Output to files every now and then
            if( (i >= chainData.burnin) && ( (i-chainData.burnin+1) % (chainData.tick*1) == 0 ) )
            {
                
                if ( chainData.output_gamma )
                {
                    gammaOutFile.open( outFilePrefix+"gamma_out.txt" , std::ios_base::trunc);
                    gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/(double)(i+1.0-chainData.burnin);
                    gammaOutFile.close();
                }
                
                if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) &&
                    chainData.output_pi )
                {
                    piOutFile.open( outFilePrefix+"pi_out.txt" , std::ios_base::trunc);
                    piOutFile << pi_out/(double)(i+1.0-chainData.burnin);
                    piOutFile.close();
                }
                
                if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
                {
                    htpOutFile.open( outFilePrefix+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
                    htpOutFile << hotspot_tail_prob_out/(double)(i+1.0-chainData.burnin);
                    htpOutFile.close();
                }
                
            }
            
            if( (i+1) % (chainData.tick*1) == 0 )
            {
                logPOutFile <<     sampler[0] -> getLogPO() <<  " ";
                logPOutFile <<     sampler[0] -> getLogPPi() <<  " ";
                logPOutFile <<     sampler[0] -> getLogPGamma() <<  " ";
                logPOutFile <<     sampler[0] -> getLogPW() <<  " ";
                logPOutFile <<     sampler[0] -> getLogLikelihood();
                logPOutFile <<     '\n';
                if ( chainData.output_model_size )
                {
                    ModelSizeOutFile << sampler[0]->getModelSize() << " ";
                    ModelSizeOutFile << '\n';
                }
            }
        }
    } // end MCMC
    
    
    // Print the end
    Rcout << " MCMC ends. " /* << " Final temperature ratio ~ " << temperatureRatio  */<< "  --- Saving results and exiting" << '\n';
    
    // ### Collect results and save them
    if ( chainData.output_gamma )
    {
        gammaOutFile.open( outFilePrefix+"gamma_out.txt" , std::ios_base::trunc);
        gammaOutFile << (arma::conv_to<arma::mat>::from(gamma_out))/(double)(chainData.nIter-chainData.burnin+1.);
        gammaOutFile.close();
    }
    
    
    logPOutFile <<     sampler[0] -> getLogPO() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPPi() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPGamma() <<  " ";
    logPOutFile <<     sampler[0] -> getLogPW() <<  " ";
    logPOutFile <<     sampler[0] -> getLogLikelihood();
    logPOutFile <<     '\n';
    logPOutFile.close();
    
    
    if ( chainData.output_model_size )
    {
        ModelSizeOutFile << sampler[0]->getModelSize();
        ModelSizeOutFile.close();
    }
    
    if ( chainData.output_model_visit )
        ModelVisitGammaOutFile.close();
    
    // ----
    if ( chainData.output_beta )
    {
        beta_out = beta_out/(double)(chainData.nIter-chainData.burnin+1);
        beta_out.save(outFilePrefix+"beta_out.txt",arma::raw_ascii);
    }
    
    if ( chainData.output_CPO )
    {
        cpo_out = cpo_out/(double)(chainData.nIter-chainData.burnin+1);
        cpo_out.save(outFilePrefix+"CPO_out.txt",arma::raw_ascii);
        
        cposumy_out = cposumy_out/(double)(chainData.nIter-chainData.burnin+1);
        cposumy_out.save(outFilePrefix+"CPOsumy_out.txt",arma::raw_ascii);
        
        waic_out = arma::log( cpo_out ) - ( waic_out - arma::square(waic_frac_sum)/(double)(chainData.nIter-chainData.burnin+1) )/(double)(chainData.nIter-chainData.burnin);
        waic_out.save(outFilePrefix+"WAIC_out.txt",arma::raw_ascii);
    }
    
    // -----
    if ( ( chainData.gamma_type == Gamma_Type::hotspot || chainData.gamma_type == Gamma_Type::hierarchical ) && chainData.output_pi )
    {
        piOutFile.open( outFilePrefix+"pi_out.txt" , std::ios_base::trunc);
        piOutFile << pi_out/(double)(chainData.nIter-chainData.burnin+1);
        piOutFile.close();
    }
    
    if ( chainData.gamma_type == Gamma_Type::hotspot && chainData.output_tail )
    {
        htpOutFile.open( outFilePrefix+"hotspot_tail_p_out.txt" , std::ios_base::trunc);
        htpOutFile << hotspot_tail_prob_out/(double)(chainData.nIter-chainData.burnin+1);
        htpOutFile.close();
    }
    // -----
    Rcout << "Saved to :   "+outFilePrefix+"****_out.txt" << '\n';
    if( chainData.surData.nFixedPredictors > 0 )
        Rcout << "Final w0 : " << sampler[0] -> getW0() <<  '\n';
    Rcout << "Final w  : " << sampler[0] -> getW() << "       w/ proposal variance: " << sampler[0] -> getVarWProposal() << '\n';
    // Rcout << "Final o : " << sampler[0] -> getO().t() << "       w/ proposal variance: " << sampler[0] -> getVarOProposal() << '\n';
    // Rcout << "Final pi : " << sampler[0] -> getPi().t() << "       w/ proposal variance: " << sampler[0] -> getVarPiProposal() << '\n';
    Rcout << "  -- Average Omega : " << arma::accu( sampler[0] -> getO() * sampler[0] -> getPi().t() )/((double)(sampler[0]->getP()*sampler[0]->getS())) <<  '\n';
    if( chainData.nChains > 1 )
        Rcout << "Final temperature ratio : " << sampler[1]->getTemperature() <<  '\n' << '\n' ;
    
    
    // Exit
    Rcout << "DONE, exiting! " << '\n' << '\n' ;
    return 0;
    
}


// *******************************************************************************
// *******************************************************************************
// *******************************************************************************


int drive( const std::string& dataFile, const std::string& mrfGFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& hyperParFile, const std::string& outFilePath,
          unsigned int nIter, unsigned int burnin, unsigned int nChains,
          const std::string& covariancePrior,
          const std::string& gammaPrior, const std::string& gammaSampler, const std::string& gammaInit,
          const std::string& betaPrior, const int maxThreads, const int tick, 
          bool output_gamma, bool output_beta, bool output_Gy, bool output_sigmaRho, bool output_pi, bool output_tail, bool output_model_size, bool output_CPO, bool output_model_visit )
{
    
    Rcout << "BayesSUR -- Bayesian Seemingly Unrelated Regression Modelling" << '\n';

    #ifdef _OPENMP
    if( maxThreads == 1 ){
        omp_set_nested( 0 );
        omp_set_num_threads( 1 );
    } else {
            Rcout << "Using OpenMP: " << maxThreads << " threads \n";
            omp_init_lock(&RNGlock);  // init RNG lock for the parallel part
        
            omp_set_nested(1); // 1=enable, 0=disable nested parallelism (run chains in parallel + compute likelihoods in parallel at least wrt to outcomes + wrt to individuals)
            // MOST OF THE PARALLELISATION IMPROVEMENTS COME FROM OPENBLAS ANYWAY .. I WONDER IF ACCELERATING LA THOURGH GPU WOULD CHANGE THAT ..
            
            // int nThreads = std::min( omp_get_max_threads()-1, maxThreads );
            // omp_set_num_threads(  std::min( omp_get_max_threads()-1, maxThreads ) );
            omp_set_num_threads( maxThreads );
        }
    #endif
    
    // ###########################################################
    // ###########################################################
    // ## Read Arguments and Data
    // ###########################################################
    // ###########################################################
    
    // Declare all the data-related variables
    Chain_Data chainData; // this initialises the pointers and the strings to ""
    
    chainData.nChains = nChains;
    chainData.nIter = nIter;
    chainData.burnin = burnin;
    
    chainData.outFilePath = outFilePath;
    
    
    // ****************************************************
    // ***  DEFINE CHAIN / PARAMETER TYPES
    // ****************************************************
    
    // Covariance type
    // *****************************************************
    
    if ( covariancePrior == "HIW" )
        chainData.covariance_type = Covariance_Type::HIW;
    else if ( covariancePrior == "IW" )
        chainData.covariance_type = Covariance_Type::IW;
    else if ( covariancePrior == "IG" )
        chainData.covariance_type = Covariance_Type::IG;
    else
    {
        Rcout << "ERROR: Wrong type of Covariance prior given\n";
        return 1;
    }
    
    // Gamma Prior
    // ****************************************************
    if ( gammaPrior == "hotspot" )
        chainData.gamma_type = Gamma_Type::hotspot ;
    else if ( gammaPrior == "hierarchical" )
        chainData.gamma_type = Gamma_Type::hierarchical ;
    else if ( gammaPrior == "MRF" )
    {
        chainData.gamma_type = Gamma_Type::mrf ;
        
    }
    else
    {
        Rcout << "ERROR: Wrong type of Gamma prior given\n";
        return 1;
    }
    
    
    // Gamma Sampler
    // ****************************************************
    if ( gammaSampler == "bandit" )
        chainData.gamma_sampler_type = Gamma_Sampler_Type::bandit ;
    else if ( gammaSampler == "MC3" )
        chainData.gamma_sampler_type = Gamma_Sampler_Type::mc3 ;
    else
    {
        Rcout << "ERROR: Wrong type of Gamma Sampler given\n";
        return 1;
    }
    
    // Beta Prior
    // *****************************************************
    if ( betaPrior == "g-prior" )
        chainData.beta_type = Beta_Type::gprior ;
    else if ( betaPrior == "independent" )
        chainData.beta_type = Beta_Type::independent ;
    else if ( betaPrior == "reGroup" )
        chainData.beta_type = Beta_Type::reGroup ;
    else
    {
        Rcout << betaPrior << "\n";
        Rcout << "ERROR: Wrong type of Beta prior given\n";
        
        return 1;
    }
    
    // need this untill g-prior code is not finalised
    if( chainData.beta_type == Beta_Type::gprior )
    {
        Rcout << "ERROR: GPrior not implemented yet\n";
        return 1;
    }
    
    
    // Outputs
    // *****************************************************
    chainData.output_gamma = output_gamma;
    chainData.output_beta = output_beta;
    chainData.output_Gy = output_Gy;
    chainData.output_sigmaRho = output_sigmaRho;
    chainData.output_pi = output_pi;
    chainData.output_tail = output_tail;
    chainData.output_model_size = output_model_size;
    chainData.output_CPO = output_CPO;
    chainData.maxThreads = maxThreads;
    chainData.tick = tick;
    chainData.output_model_visit = output_model_visit;
    
    // ***********************************
    // ***********************************
    
    
    // read Data and format into usables
    Rcout << "Reading input files ... ";
    
    try
    {
        Utils::formatData(dataFile, mrfGFile, blockFile, structureGraphFile, chainData.surData );
    }
    catch(const std::exception& e)
    {
        Rcerr << e.what() << '\n';
        return 1;
    }
    
    try
    {
        Utils::readHyperPar(hyperParFile, chainData );
    }
    catch(const std::exception& e)
    {
        Rcerr << e.what() << '\n';
        return 1;
    }
    
    Rcout << "... successfull!" << '\n';
    
    // ############
    
    Rcout << "Clearing and initialising output files " << '\n';
    
    // Re-define dataFile so that I can use it in the output
    chainData.filePrefix = dataFile;
    std::size_t slash = chainData.filePrefix.find("/");  // remove the path from filePrefix
    while( slash != std::string::npos )
    {
        chainData.filePrefix.erase(chainData.filePrefix.begin(),chainData.filePrefix.begin()+slash+1);
        slash = chainData.filePrefix.find("/");
    }
    chainData.filePrefix.erase(chainData.filePrefix.begin()+chainData.filePrefix.find(".txt"),chainData.filePrefix.end());  // remove the .txt from filePrefix !
    
    // Update the "outFilePath" (filePrefix variable) with the method's name
    chainData.filePrefix += "_";
    
    switch ( chainData.covariance_type )
    {
        case Covariance_Type::HIW :
            chainData.filePrefix += "SSUR_";
            break;
            
        case Covariance_Type::IW :
            chainData.filePrefix += "dSUR_";
            break;
            
        case Covariance_Type::IG :
            chainData.filePrefix += "HRR_";
            break;
            
        default:
            throw Bad_Covariance_Type( chainData.covariance_type );
    }
    
    /*#ifdef _OPENMP
            // ENABLING NESTED PARALLELISM SEEMS TO SLOW DOWN CODE MORE THAN ANYTHING,
            // I SUSPECT THE THREAD MANAGING OVERHEAD IS GREATER THAN EXPECTED
            omp_set_nested(1); // 1=enable, 0=disable nested parallelism (run chains in parallel + compute likelihoods in parallel at least wrt to outcomes + wrt to individuals)
            // MOST OF THE PARALLELISATION IMPROVEMENTS COME FROM OPENBLAS ANYWAY .. I WONDER IF ACCELERATING LA THOURGH GPU WOULD CHANGE THAT ..
            
            int nThreads = std::min( omp_get_max_threads()-1, maxThreads );
            // omp_set_num_threads(  std::min( omp_get_max_threads()-1, maxThreads ) );
            omp_set_num_threads( nThreads );
    #endif*/
    
    // ###################################
    // Parameters Inits
    // ###################################
    
    if ( gammaInit == "R" )
    {
        // Random Init
        chainData.gammaInit = arma::umat(chainData.surData.nVSPredictors,chainData.surData.nOutcomes); // init empty
        for(unsigned int j=0; j<chainData.surData.nVSPredictors; ++j)
            for(unsigned int l=0; l< chainData.surData.nOutcomes; ++l)
                chainData.gammaInit(j,l) = randBernoulli( 0.1 );
        
    }else if( gammaInit == "1" ){
        // Static Init ***
        // ** 1
        chainData.gammaInit = arma::ones<arma::umat>(chainData.surData.nVSPredictors,chainData.surData.nOutcomes);
        
    }else if ( gammaInit == "0" ) {
        // ** 0
        chainData.gammaInit = arma::zeros<arma::umat>(chainData.surData.nVSPredictors,chainData.surData.nOutcomes);
        
    }else if ( gammaInit == "MLE" ) {
        // ** MLE
        arma::mat Q,R; arma::qr(Q,R, chainData.surData.data->cols( arma::join_vert( *chainData.surData.fixedPredictorsIdx , *chainData.surData.VSPredictorsIdx ) ) );
        
        chainData.betaInit = arma::solve(R,arma::trans(Q) * chainData.surData.data->cols( *chainData.surData.outcomesIdx ) );
        chainData.gammaInit = chainData.betaInit > 0.5*arma::stddev(arma::vectorise(chainData.betaInit));
        
        if( chainData.surData.nFixedPredictors > 0 )
            chainData.gammaInit.shed_rows( 0 , chainData.surData.nFixedPredictors-1 ); // shed the fixed preditors rows since we don't have gammas for those
        
    }else{
        // default case
        chainData.gammaInit = arma::zeros<arma::umat>(chainData.surData.nVSPredictors,chainData.surData.nOutcomes);
    }
    
    // ###################################
    // Samplers
    // ###################################
    
    int status;
    
    // TODO, I hate this, but I can't initialise/instanciate templated classes
    // at runtime so this seems fair (given that the different drive functions have their differences in output and stuff...)
    // still if there's a more elegant solution I'd like to find it
    
    try
    {
        switch ( chainData.covariance_type )
        {
            case Covariance_Type::HIW :
            case Covariance_Type::IW :
                status = drive_SUR(chainData);
                break;
                
            case Covariance_Type::IG :
                status = drive_HRR(chainData);
                break;
                
            default:
                status = 1;
                throw Bad_Covariance_Type( chainData.covariance_type );
        }
    }
    catch(...){ }
    
    return status;
}
