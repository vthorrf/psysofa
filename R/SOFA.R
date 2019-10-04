SOFA <- function(x,y,method="bayes") {
  if(method == "bayes") {
    source("BayesSPSE.R")
    result <- BayesSPSE(x,y)
    Results <- list("disposition"=result$eff,"full"=result,
                    "pred"=result$pred,"dic"=result$dic)
  } else if (method == "gam") {
    require(semsfa)
    dati <- data.frame(x,y)
    result <- semsfa::semsfa(y~s(x),dati,sem.method="gam",ineffDecrease=FALSE)
    Results <- list("disposition"=semsfa::efficiencies.semsfa(result)$efficiencies,
                    "full"=result,"pred"=fitted(result),"bic"=result$bic)
  } else if (method == "kernel") {
    require(semsfa)
    dati <- data.frame(x,y)
    result <- semsfa::semsfa(y~x,dati,sem.method="kernel",ineffDecrease=FALSE)
    Results <- list("disposition"=semsfa::efficiencies.semsfa(result)$efficiencies,
                    "full"=result,"pred"=fitted(result),"bic"=result$bic)
  } else if (method == "loess") {
    require(semsfa)
    dati <- data.frame(x,y)
    result <- semsfa::semsfa(y~x,dati,sem.method="loess",ineffDecrease=FALSE)
    Results <- list("disposition"=semsfa::efficiencies.semsfa(result)$efficiencies,
                    "full"=result,"pred"=fitted(result),"bic"=result$bic)
  } else { warning("Unknow method. Please refer to SOFA's documentation.") }
  return(Results)
}
