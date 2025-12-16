#' @title Degrees of Freedom Calculations

#' @description
#' 
#' Methods for calculating degrees of freedom:
#' 
#' \code{residual_aov}: DFs are calculated from `df.residual`
#' \code{asymptotic_aov}: DFs are set to \code{Inf} for all model terms
#' \code{nlme_aov}: DFs are calculated based on the algorithm given in Mixed-Effects Models in S and S-Plus and employed in the \code{nlme} package
#' \code{containment_aov}: DFs are calculated using the containment algorithm
#' \code{satterthwaite_aov}: DFs are calculated using the Satterthwaite approximation, as implemented in \code{glmmTMB::dof_satt}
#' 
#' @param model a \code{\link{glmmTMB}} model object
#' @param type type of test, \code{II, III, 2, 3}. Roman numerals are equivalent
#'     to the corresponding Arabic numerals.
#'
#' @return a \code{data.frame}
#' @export
residual_aov = function(model = model, type = type){
  y_name<- names(model$modelInfo$respCol)
  dc <- dataClasses(model)
  TMBaov <- suppressPackageStartupMessages(car::Anova(model, type=type))
  
  # Pull the DFs associated with each term
  basic_aov_dfs <- base_aov_dfs(model)
  
  # Setup the final output
  if(row.names(TMBaov)[1] == "(Intercept)") {
    df_output <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ]
  } else {
    df_output <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ][-1, ]
  }
  
  df_output$vartype <- NULL
  df_output$denDf <- df.residual(model)
  
  chisq <- as.vector(TMBaov$Chisq)
  nDF <- as.vector(TMBaov$Df)
  Fval <- chisq/nDF
  dDF <- df_output$denDf
  Pval <- pf(Fval, nDF, dDF, lower.tail = FALSE)
  
  aod <- data.frame(numDF = nDF, denDF = dDF, Fvalue = round(Fval, 2), pvalue = round(Pval, 4))
  row.names(aod) <- df_output$terms
  class(aod) <- c("bdf", "residual", "data.frame")
  
  if (type == 3) {
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type III F-tests)", "\n\nResponse: ", y_name)
  } else if (type == 2){
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type II F-tests)", "\n\nResponse: ", y_name)
  }
  
  return(aod)
}

#' @rdname residual_aov
#' @export
asymptotic_aov = function(model = model, type = type){
  y_name<- names(model$modelInfo$respCol)
  dc <- dataClasses(model)
  TMBaov <- suppressPackageStartupMessages(car::Anova(model, type=type))
  
  # Pull the DFs associated with each term
  basic_aov_dfs <- base_aov_dfs(model)
  
  # Setup the final output
  if(row.names(TMBaov)[1] == "(Intercept)") {
    df_output <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ]
  } else {
    df_output <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ][-1, ]
  }
  
  df_output$vartype <- NULL
  df_output$denDf <- Inf
  
  chisq <- as.vector(TMBaov$Chisq)
  nDF <- as.vector(TMBaov$Df)
  Fval <- chisq/nDF
  dDF <- df_output$denDf
  Pval <- pf(Fval, nDF, dDF, lower.tail = FALSE)
  
  aod <- data.frame(numDF = nDF, denDF = dDF, Fvalue = round(Fval, 2), pvalue = round(Pval, 4))
  row.names(aod) <- df_output$terms
  class(aod) <- c("bdf", "asymptotic", "data.frame")
  
  if (type == 3) {
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type III F-tests)", "\n\nResponse: ", y_name)
  } else if (type == 2){
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type II F-tests)", "\n\nResponse: ", y_name)
  }
  
  return(aod)
}
