#' meta_path.RData
#' @param qvec This is the Metabolite_pathway_table from MetaTOpathway function. This table includes the metabolites ids and the its corssponding hmdb ids
#' @keywords PDS
#' @import dplyr pathifier
#' @export
#' @examples
#' \donttest{
#' @@ -20,22 +23,64 @@
#' }

lilikoi.PDSfun<-function(qvec){
  Metadata<-as.data.frame(Metadata)
  Metadata$label<-as.factor(Metadata$label)
  phe=Metadata$label %>% as.numeric %>% -1
  newData1=qvec %>% filter(pathway!='NA')%>% dplyr::select(Query,HMDB)
  newData=Metadata[,t(newData1['Query'])]
  colnames(newData)=t(newData1['HMDB'])
  newData=t(newData)
  newData <- matrix(as.numeric(newData),
                    nrow = nrow(newData),
                    ncol = ncol(newData),
                    dimnames = dimnames(newData))
  PDS<-PDSchanged(as.matrix(newData), row.names(newData), metabolites.list,pathway.list,normal=as.logical(phe),attempts = 5, min_exp=0, min_std=0)
  qpdmat <- matrix(as.data.frame(PDS$scores), nrow=length(names(PDS$scores)), byrow=TRUE)
  #qpdmat <<- data.frame(lapply(PDS$scores,function (x) {as.numeric(unlist(x))}),check.names=F)
  colnames(qpdmat) <- colnames(newData)
  rownames(qpdmat) <- names(PDS$scores)
  mode(qpdmat) <- "numeric"
  return(PDS)
}
