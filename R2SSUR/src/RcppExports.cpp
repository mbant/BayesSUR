// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// R2SSUR_internal
int R2SSUR_internal(const std::string& dataFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& outFilePath, unsigned int nIter, unsigned int burnin, unsigned int nChains, const std::string& method, bool sparse, const std::string& gammaPrior, const std::string& gammaSampler, const std::string& gammaInit, const std::string& mrfGFile, const std::string& betaPrior);
RcppExport SEXP _R2SSUR_R2SSUR_internal(SEXP dataFileSEXP, SEXP blockFileSEXP, SEXP structureGraphFileSEXP, SEXP outFilePathSEXP, SEXP nIterSEXP, SEXP burninSEXP, SEXP nChainsSEXP, SEXP methodSEXP, SEXP sparseSEXP, SEXP gammaPriorSEXP, SEXP gammaSamplerSEXP, SEXP gammaInitSEXP, SEXP mrfGFileSEXP, SEXP betaPriorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::string& >::type dataFile(dataFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type blockFile(blockFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type structureGraphFile(structureGraphFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type outFilePath(outFilePathSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nChains(nChainsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type method(methodSEXP);
    Rcpp::traits::input_parameter< bool >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type gammaPrior(gammaPriorSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type gammaSampler(gammaSamplerSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type gammaInit(gammaInitSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type mrfGFile(mrfGFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type betaPrior(betaPriorSEXP);
    rcpp_result_gen = Rcpp::wrap(R2SSUR_internal(dataFile, blockFile, structureGraphFile, outFilePath, nIter, burnin, nChains, method, sparse, gammaPrior, gammaSampler, gammaInit, mrfGFile, betaPrior));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_R2SSUR_R2SSUR_internal", (DL_FUNC) &_R2SSUR_R2SSUR_internal, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_R2SSUR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
