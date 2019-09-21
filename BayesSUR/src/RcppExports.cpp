// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// BayesSUR_internal
int BayesSUR_internal(const std::string& dataFile, const std::string& mrfGFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& hyperParFile, const std::string& outFilePath, unsigned int nIter, unsigned int burnin, unsigned int nChains, const std::string& covariancePrior, const std::string& gammaPrior, const std::string& gammaSampler, const std::string& gammaInit, const std::string& betaPrior, bool output_gamma, bool output_beta, bool output_G, bool output_sigmaRho, bool output_pi, bool output_tail, bool output_model_size, bool output_CPO);
RcppExport SEXP _BayesSUR_BayesSUR_internal(SEXP dataFileSEXP, SEXP mrfGFileSEXP, SEXP blockFileSEXP, SEXP structureGraphFileSEXP, SEXP hyperParFileSEXP, SEXP outFilePathSEXP, SEXP nIterSEXP, SEXP burninSEXP, SEXP nChainsSEXP, SEXP covariancePriorSEXP, SEXP gammaPriorSEXP, SEXP gammaSamplerSEXP, SEXP gammaInitSEXP, SEXP betaPriorSEXP, SEXP output_gammaSEXP, SEXP output_betaSEXP, SEXP output_GSEXP, SEXP output_sigmaRhoSEXP, SEXP output_piSEXP, SEXP output_tailSEXP, SEXP output_model_sizeSEXP, SEXP output_CPOSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::string& >::type dataFile(dataFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type mrfGFile(mrfGFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type blockFile(blockFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type structureGraphFile(structureGraphFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type hyperParFile(hyperParFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type outFilePath(outFilePathSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nChains(nChainsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type covariancePrior(covariancePriorSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type gammaPrior(gammaPriorSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type gammaSampler(gammaSamplerSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type gammaInit(gammaInitSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type betaPrior(betaPriorSEXP);
    Rcpp::traits::input_parameter< bool >::type output_gamma(output_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type output_beta(output_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type output_G(output_GSEXP);
    Rcpp::traits::input_parameter< bool >::type output_sigmaRho(output_sigmaRhoSEXP);
    Rcpp::traits::input_parameter< bool >::type output_pi(output_piSEXP);
    Rcpp::traits::input_parameter< bool >::type output_tail(output_tailSEXP);
    Rcpp::traits::input_parameter< bool >::type output_model_size(output_model_sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type output_CPO(output_CPOSEXP);
    rcpp_result_gen = Rcpp::wrap(BayesSUR_internal(dataFile, mrfGFile, blockFile, structureGraphFile, hyperParFile, outFilePath, nIter, burnin, nChains, covariancePrior, gammaPrior, gammaSampler, gammaInit, betaPrior, output_gamma, output_beta, output_G, output_sigmaRho, output_pi, output_tail, output_model_size, output_CPO));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesSUR_BayesSUR_internal", (DL_FUNC) &_BayesSUR_BayesSUR_internal, 22},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesSUR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
