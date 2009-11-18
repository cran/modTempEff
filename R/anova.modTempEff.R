`anova.modTempEff`<-function(object, ..., dispersion = NULL, test = NULL){
#----
    anova.glmlist1<-function (object, ..., dispersion = NULL, test = NULL){
      responses <- as.character(lapply(object, function(x) {
        deparse(formula(x)[[2L]])
    }))
    sameresp <- responses == responses[1L]
    if (!all(sameresp)) {
        object <- object[sameresp]
        warning("models with response ", deparse(responses[!sameresp]),
            " removed because response differs from model 1")
    }
    ns <- sapply(object, function(x) length(x$residuals))
    if (any(ns != ns[1L]))
        stop("models were not all fitted to the same size of dataset")
    nmodels <- length(object)
    if (nmodels == 1)
        return(anova.glm(object[[1L]], dispersion = dispersion,
            test = test))
    resdf <- as.numeric(lapply(object, function(x) x$df.residual))
    resdev <- as.numeric(lapply(object, function(x) x$deviance))
    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA,
        -diff(resdev)))
    variables <- lapply(object, function(x) paste(deparse(formula(x)),
        collapse = "\n"))
    dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev",
        "Df", "Deviance"))
    title <- "Analysis of Deviance Table\n"
    topnote <- paste("Model ", format(1L:nmodels), ": ", variables,
        sep = "", collapse = "\n")
    if (!is.null(test)) {
        bigmodel <- object[[order(resdf)[1L]]]
        dispersion <- 1#summary(bigmodel, dispersion = dispersion)$dispersion
        df.dispersion <- if (dispersion == 1)
            Inf
        else min(resdf)
        if (test == "F" && df.dispersion == Inf) {
            fam <- bigmodel$family$family
            if (fam == "binomial" || fam == "poisson")
                warning(gettextf("using F test with a '%s' family is inappropriate",
                  fam), domain = NA, call. = FALSE)
            else warning("using F test with a fixed dispersion is inappropriate")
        }
        table <- stat.anova(table = table, test = test, scale = dispersion,
            df.scale = df.dispersion, n = length(bigmodel$residuals))
    }
    structure(table, heading = c(title, topnote), class = c("anova",
        "data.frame"))
  }
#start anova.modTempEff
    dotargs <- list(...)
    named <- if (is.null(names(dotargs)))
        rep(FALSE, length(dotargs))
    else (names(dotargs) != "")
    if (any(named))
        warning("The following arguments to anova.modTempEff(..) are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.modTempEff <- unlist(lapply(dotargs, function(x) inherits(x, "modTempEff")))
    dotargs <- dotargs[is.modTempEff]
    if (length(dotargs) <= 0) stop("this anova method for single object is not allowed")
    return(anova.glmlist1(c(list(object), dotargs), dispersion = dispersion,
            test = test))
    }
    

