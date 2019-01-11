#' Generalized Mass Spectrum missing peaks imputation with Two-Step Lasso as default algorithm
#'
#' GMS.Lasso recovers the abundance of missing peaks via either TS.Lasso or the minimum abundance per compound.  
#' @param input_data Raw abundance matrix with missing value, with features in rows and samples in columns.
#' @param alpha Weights for L1 penalty in Elastic Net. The default and suggested value is alpha=1, which is for Lasso. 
#' @param nfolds The number of folds used in parameter (lambda) tuning. 
#' @param log.scale Whether the input_data needs log scale transform.The default is log.scale=T, assuming input_data is the 
#' raw abundance matrix. If input_data is log abundance matrix, log.scale=F.  
#' @param TS.Lasso Whether to use TS.Lasso or the minimum per compound for imputation. 
#' @return 
#'    \item{imputed.final}{The imputed abundance matrix at the scale of input_data.}
#' 
#'    
#' @examples 
#' data('tcga.bc')
#' # tcga.bc contains mass specturm abundance of 150 metabolites for 30 breast cancer 
#' # tumor and normal tissue samples with missing values.
#' 
#' imputed.compound.min=GMS.Lasso(tcga.bc,log.scale=TRUE,TS.Lasso=FALSE)
#' # Impute raw abundance matrix tcga.bc with compound minimum
#' 
#' imputed.tslasso=GMS.Lasso(tcga.bc,log.scale=TRUE,TS.Lasso=TRUE)
#' # Impute raw abundance matrix tcga.bc with TS.Lasso
#' 
#' @export 

GMS.Lasso<-function(input_data,alpha=1,nfolds=10,log.scale=TRUE,TS.Lasso=TRUE){
  row.missing=rowSums(is.na(input_data))  
  if(sum(row.missing==0)<5){
    cat('TS.Lasso requires at least 5 compounds without missing peaks. Switched to TS.Lasso=FALSE ')
    TS.Lasso=FALSE
  }
  if(TS.Lasso){
    imputed.final=TS.Lasso(input_data,alpha,nfolds,log.scale)
  }else{
    cat('Impute by the minimum per compound')
    imputed.final=apply(input_data,1,function(x){
      y=ifelse(is.na(x),min(x,na.rm = T),x)
      return(y)
    })
    
  }
  return(imputed.final)
  
}