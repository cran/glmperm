prr.test <-
function(formula, var, family=gaussian, data, nrep = 1000, seed=12345, weights,  subset, na.action,
    start = NULL, etastart, mustart, offset, control = glm.control(...), model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ...)
{
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data",  "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), glm.fit = 1, stop("invalid 'method': ",
        method))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    if(!(paste(var) %in% colnames(X))) stop("var not a covariate in the formular")
    X <- cbind(X, rep(NA,nrow(X)))
    colnames(X)[ncol(X)] <- "resid"            
    X[,"resid"] <- lm.fit(x = X[, -which(colnames(X) %in% c(paste(var),"resid")), drop = FALSE], y=X[, paste(var), drop = FALSE])$residuals
    ###
    fit1 <- glm.fit(x = X[, -which(colnames(X)==paste(var)), drop = FALSE], y = Y, weights = weights, start = start,
        etastart = etastart, mustart = mustart, offset = offset, family = family, control = control, intercept = attr(mt,
            "intercept") > 0)
    fit2 <- glm.fit(x = X[, -which(colnames(X) %in% c(paste(var),"resid")), drop = FALSE], y = Y, weights = weights, start = start,
        etastart = etastart, mustart = mustart, offset = offset, family = family, control = control, intercept = attr(mt,
            "intercept") > 0)
    ###  Dispersion factor for model fit1
    df.r <- fit1$df.residual
      if(df.r > 0){
            if(any(fit1$weights == 0)){warning("observations with zero weight not used for calculating dispersion")}
         dispersion <- sum((fit1$weights * fit1$residuals^2)[fit1$weights > 0])/df.r
        }
       if(df.r==0){ dispersion <- NaN}
    if(df.r==0){warning("dispersion is Na")}
    if(fit1$family$family=="binomial" & dispersion>1.5){warning("estimated dispersion is > 1.5, rather use family = quasibinomial")}
    if(fit1$family$family=="binomial" & dispersion<0.5){warning("estimated dispersion is < 0.5, rather use family = quasibinomial")}
    if(fit1$family$family=="poisson" & dispersion>1.5){warning("estimated dispersion is > 1.5, rather use family = quasipoisson")}
    if(fit1$family$family=="poisson" & dispersion<0.5){warning("estimated dispersion is < 0.5, rather use family = quasipoisson")}
    estimated.dispersion <- dispersion
    if(fit1$family$family %in% c("poisson", "binomial")){dispersion <- 1}
    ###
    p.value.obs <- 1 - pchisq(abs(fit1$deviance - fit2$deviance)/dispersion, 1)
    ### permutation part
    set.seed(seed)
    devi.disp <- matrix(0, ncol=2, nrow=nrep)
    options(warn = -1)
    oldtime <- proc.time()[1]
    for (i in 1:nrep){devi.disp[i,] <- glm.perm(Y, X[, -which(colnames(X)==paste(var)), drop = FALSE], Family=family)}
    print(c("execution time in minutes", round((proc.time()[1] - oldtime)/60, 2)))
    options(warn = 0)
    psim <- 1 - pchisq(abs(devi.disp[,1] - fit2$deviance)/devi.disp[,2], 1)
    ### output of prr.test by Potter 
    ret.val <- list(nobs = nrow(X), p0 = length(psim[psim <= p.value.obs])/nrep,
        p005 = length(psim[psim <= 1.005 * p.value.obs])/nrep,
        p01 = length(psim[psim <= 1.01 * p.value.obs])/nrep,
        p02 = length(psim[psim <= 1.02 * p.value.obs])/nrep,
        p04 = length(psim[psim <= 1.04 * p.value.obs])/nrep)
    names(ret.val$nobs) <- "number of observations used"
    names(ret.val$p0) <- "permutation p-value for simulated p-values <= observed p-value"
    names(ret.val$p005) <- "permutation p-value for simulated p-values <= 1.005 observed p-value"
    names(ret.val$p01) <- "permutation p-value for simulated p-values <= 1.01 observed p-value"
    names(ret.val$p02) <- "permutation p-value for simulated p-values <= 1.02 observed p-value"
    names(ret.val$p04) <- "permutation p-value for simulated p-values <= 1.04 observed p-value"
    ### new output
    if (model)
        fit1$model <- mf
    fit1$na.action <- attr(mf, "na.action")
    if (x)
        fit1$x <- X
    if (!y)
        fit1$y <- NULL
    out <- c(list(fit1 = c(fit1,list(terms = mt, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,mf))), fit2 = fit2, call = call, formula = formula, 
        seed=seed, fit1deviance=fit1$deviance, fit2deviance=fit2$deviance, Dispersion=dispersion, estimated.Dispersion = estimated.dispersion,
        LRstat= abs(fit1$deviance - fit2$deviance)/dispersion,  p.value.obs=p.value.obs, p.value.perm = ret.val, nobs= ret.val$nobs, var=var))
   class(out) <- "prr.test"
    return(out)
  }

