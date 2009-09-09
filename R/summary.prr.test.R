summary.prr.test <-
function(object, digits = max(3, getOption("digits") - 3),...)
{
cat("\n\n    Permutation of Regressor Residuals Test:\n\n\n")
cat("Call: \n", deparse(object$call), "\n\n")
cat("number of observations used: ",  format(signif(object$nobs, digits)),"\n\n")
cat("null hypothesis: regression coefficient of covariate", object$var,"= 0")
cat("\nobserved Likelihood Ratio Test Statistics: ", format(signif(object$LRstat, digits)),"\n\n")
cat("    -----------------------------------------\n") 
cat("    Results based on chi-squared distribution \n")   
cat("    -----------------------------------------\n")  
cat("\nobserved p-value:", format(signif(object$p.value.obs, digits)),"\n\n")
cat("    ---------------------------------------------------\n") 
cat("    Results based on permutation of regressor residuals \n")   
cat("    ---------------------------------------------------\n")  
cat("\npermutation p-value for simulated p-values <= observed p-value:",  format(signif(object$p.value.perm$p0, digits)),"\n")
cat("\npermutation p-value for simulated p-values <= 1.005 observed p-value:",  format(signif(object$p.value.perm$p005, digits)),"\n")
cat("\npermutation p-value for simulated p-values <= 1.01 observed p-value:",  format(signif(object$p.value.perm$p01, digits)),"\n")
cat("\npermutation p-value for simulated p-values <= 1.02 observed p-value:",  format(signif(object$p.value.perm$p02, digits)),"\n")
cat("\npermutation p-value for simulated p-values <= 1.04 observed p-value:",  format(signif(object$p.value.perm$p04, digits)),"\n\n")
## Warnings for Over and Underdispersion 
if(object$fit1$family$family=="binomial" & object$estimated.Dispersion >1.5){
cat("*************************************************************************\n")
cat("WARNING: estimated dispersion is > 1.5, rather use family = quasibinomial\n")
cat("*************************************************************************\n")}
if(object$fit1$family$family=="binomial" & object$estimated.Dispersion<0.5){
cat("*************************************************************************\n")
cat("WARNING: estimated dispersion is < 0.5, rather use family = quasibinomial\n")
cat("*************************************************************************\n")}
if(object$fit1$family$family=="poisson" & object$estimated.Dispersion>1.5){
cat("************************************************************************\n")
cat("WARNING: estimated dispersion is > 1.5, rather use family = quasipoisson\n")
cat("************************************************************************\n")}
if(object$fit1$family$family=="poisson" & object$estimated.Dispersion<0.5){
cat("************************************************************************\n")
cat("WARNING: estimated dispersion is < 0.5, rather use family = quasipoisson\n")
cat("************************************************************************\n")}
    invisible(object)
}

