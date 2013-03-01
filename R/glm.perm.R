glm.perm <-
function(y, x, Family)
{   
    ## resample the data
    sel <- sample(1:nrow(x), replace = FALSE)
    x.resample <- cbind(x[, -which(colnames(x)=="resid"), drop = FALSE], x[, "resid", drop = FALSE][sel])
    ## calculate deviances and dispersion factors for resampled data
    if(Family$family == "Cox"){
      f1 <- try(coxph.fit(x=x.resample, y=y, strata=NULL, control=coxph.control(), method="efron", rownames=NULL))
      devi <- -2*f1$loglik[2]
      return(devi)  
    } else {  
      f1 <- glm.fit(x.resample, y, family=Family)
      df.r <- f1$df.residual
      if(f1$family$family %in% c("poisson", "binomial")){disp <- 1}
      if(all(f1$family$family != c("poisson", "binomial"))){
        if(df.r > 0){
              if(any(f1$weights == 0)){warning("observations with zero weight not used for calculating dispersion")}
           disp <- sum((f1$weights * f1$residuals^2)[f1$weights > 0])/df.r
          }
         if(df.r==0){disp <- NaN}
        }
      devi <- f1$deviance
      return(c(devi, disp))
    }
}


