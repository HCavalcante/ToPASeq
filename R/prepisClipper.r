clipperSingle<-function(expr, classes, g, nperm, alphaV, b, permute, method){

out<-list(unlist(pathQ(expr, classes, g, nperm=nperm, alphaV=alphaV, b=nperm, permute=TRUE)),
  easyClip(expr, classes, g, method ="mean"))


 return(out)
}  
        

CLIPPER<-function (pathways, expr, classes, method, testCliques=FALSE, nperms, alphaV=0.05, b=NULL, permute=TRUE){
    #if (!require(clipper))        stop("library clipper is missing")

    classes<-factor(as.numeric(factor(classes)))
message("Analysing pathway:\n")
      out<-catchErr(pathways, function(p) {
      cat(p[[2]],"\n")
clipperSingle(expr, classes, p[[1]], nperms, alphaV, b, permute, method)
})
cliq.test<-list()
if (testCliques) {
 cliq.test<-lapply(pathways, function(p) {
 if (method=="mean") cliq<-cliqueMeanTest(expr, classes, p[[1]], nperms, alphaV, b, permute=permute)
 if (method=="var") cliq<-cliqueVarianceTest(expr, classes, p[[1]], nperms, alphaV, b, permute=permute)
 return(cliq)
 })
 
  }
if (length(out[[1]])>0) {
res<-data.frame(t(sapply(out[[1]],function(x) x[[1]])))
res$mean.q.value<-p.adjust(res$alphaMean,"fdr")
res$var.q.value<-p.adjust(res$alphaVar,"fdr")

paths<-lapply(out[[1]],function(x) x[[2]])

out[[1]]<-list(res, paths, cliq.test)
}
return(out)
}

    
