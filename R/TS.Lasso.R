#' Two-Step Lasso for missing peaks imputation
#'
#' TS.Lasso recovers the abundance of various types of missing peaks.
#' @param input_data Raw abundance matrix with missing value, with features in rows and samples in columns.
#' @param alpha Weights for L1 penalty in Elastic Net. The default and suggested value is alpha=1, which is for Lasso.
#' @param nfolds The number of folds used in parameter (lambda) tuning.
#' @param log.scale Whether the input_data needs log scale transform.The default is log.scale=T, assuming input_data is the
#' raw abundance matrix. If input_data is log abundance matrix, set log.scale=F.
#'
#' @return
#'    \item{imputed.final}{The imputed abundance matrix at the scale of input_data.}
#'
#' @import utils
#' @import glmnet
#'
#'
#' @examples
#' data('tcga.bc')
#' # tcga.bc contains mass specturm abundance of 150 metabolites for 30 breast cancer
#' # tumor and normal tissue samples with missing values.
#'
#' imputed=TS.Lasso(tcga.bc,log.scale=TRUE)
#' # Impute raw abundance matrix tcga.bc
#'
#' @export
#'


TS.Lasso<-function(input_data,alpha=1,nfolds=10,log.scale=TRUE){
  if(!requireNamespace('glmnet')) {
    install.packages('glmnet')
    if(!requireNamespace('glmnet',character.only = TRUE)) stop("Package 'glmnet' not found")}

  row.missing=rowSums(is.na(input_data))
  if(sum(row.missing==0)<5){
    stop('Need at least 5 compounds without missing peaks')
  }


  if(sum(row.missing<0.5*ncol(input_data) & row.missing>0 )==0){
    stop('No missing value detected after removing rows with more than 50% missing')
  }

  cat('Filtering: remove rows with more than 50% missing')
  input_data=input_data[row.missing<0.5*ncol(input_data),]

  if(log.scale){
      input_data=log2(input_data)
    }

    if(is.null(rownames(input_data))){
      rownames(input_data)=paste('Var',1:dim(input_data)[1],sep='')
    }
    if(is.null(colnames(input_data))){
      colnames(input_data)=paste('sample',1:dim(input_data)[2],sep='')
    }

      miss=input_data[rowSums(is.na(input_data))>0,]
      nomiss=input_data[rowSums(is.na(input_data))==0,]
      imputed=miss
    for(i in 1:dim(miss)[1]){
      cvfit = cv.glmnet(t(nomiss[,!is.na(miss[i,])]), t(miss[i,][!is.na(miss[i,])]),alpha=alpha,nfolds = nfolds,family='gaussian')
      imputed[i,][is.na(miss[i,])]=predict(cvfit, newx = t(nomiss[,is.na(miss[i,])]),s = "lambda.min")
    }

    imputed.second=imputed
    for(i in 1:dim(miss)[1]){
      imputed.all=rbind(nomiss,imputed[-i,])
      cvfit = cv.glmnet(t(imputed.all[,!is.na(miss[i,])]), t(miss[i,][!is.na(miss[i,])]),alpha=alpha,nfolds = nfolds,family='gaussian')
      imputed.second[i,][is.na(miss[i,])]=predict(cvfit, newx = t(imputed.all[,is.na(miss[i,])]),s = "lambda.min")
    }

    imputed.final=rbind(nomiss,imputed.second)
    imputed.final=imputed.final[rownames(input_data),]
    if(log.scale){
      imputed.final=2^imputed.final
    }
    return(imputed.final)
  }
