#' after the Actual Impurity Reduction Importance is computed with a Random Forest, pvalues from different probes or SNPs belonging to the same gene can be aggregated in order to obtain a single pvalue for that gene. Correlation between probes can also be taken into account.
#'
#' @param rf a ranger object with "importance="impurity_corrected"
#' @param genenames a vector of the name of the gene to which each probe or SNP belongs, it has to be of size dim(x)[1]
#' @param x design matrix used by the random forest
#' @param method one of Tippett, Fisher, Kost, EBM
#' @param adjust "no" / "yes" depending if correlation has to be taken into account
#' @import ranger
#' @import EmpiricalBrownsMethod

#' @export
#' @examples
#' n <- 250
#' x=replicate(50, runif(n))
#' dat <- data.frame(y = factor(rbinom(n, 1, .5)), x)
#' library(ranger)
#' rf <- ranger(y ~ ., dat, importance = "impurity_corrected",num.trees=100)
#' genenames=colnames(x)=rep(c("G1","G2"),50/2)
#' poolVIM(rf,genenames,x,method="Fisher",adjust="no")


#'
poolVIM=function (rf,genenames,x,method="Tippett",adjust){

  imp=rf$variable.importance
  null=c(imp[imp<0],imp[imp=0],-1*imp[imp<0])
  p=  sapply(1:length(imp), function(i) {

    pnorm(gaussianize(null,imp[i]),lower.tail = FALSE)
  })

  #p[which(p==0)]=1/(length(null)+1)

  p[which(p==0)]=.Machine$double.eps
  names(p)=genenames

  genenames=names(table(genenames))


  if (method == "Tippett") {
  res=sapply(1:length(genenames), function(i) {
    tippett(p[which(names(p)==genenames[i])],adjust=adjust,R=cor(t(x[,which(names(p)==genenames[i])])))
  })

}

  if (method == "Fisher") {

  res=sapply(1:length(genenames), function(i) {
    fisher(p[which(names(p)==genenames[i])],adjust=adjust, R=cor(t(x[,which(names(p)==genenames[i])])))
  })
}



if (method == "EBM") {

res=sapply(1:length(genenames), function(i) {


    empiricalBrownsMethod(t(x[,which(colnames(x)==genenames[i])]),as.vector(p[which(names(p)==genenames[i])]))


})
}

  if (method == "Kost") {

    res=sapply(1:length(genenames), function(i) {


      kostsMethod(t(x[,which(colnames(x)==genenames[i])]),as.vector(p[which(names(p)==genenames[i])]))


    })
  }
names(res)=genenames
return(res)
}
