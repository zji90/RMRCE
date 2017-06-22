#' RMRCE_cv
#' 
#' Cross validation for choosing RMRCE tuning parameters
#'
#' This function performs a five-fold cross validation to choose the optimal tuning parameters lambda and alpha which yield the smallest testing error.
#'
#' @param X A numeric matrix of explanatory variables.
#' @param Y A numeric vector of response.
#' @param lambda A numeric vector of tuning parameter lambda for which the testing error is to be computed.
#' @param alpha A numeric vector of tuning parameter alpha for which the testing error is to be computed.
#' @return A data.frame which has three columns: lambda, alpha and corresponding testing error. Lambda and alpha that yield the smallest testing error are selected as the optimal tuning parameters.
#' @export
#' @author Fang Han, Hongkai Ji, Zhicheng Ji, Honglang Wang <zji4@@jhu.edu>
#' @examples
#' Y <- rnorm(100)
#' X <- cbind(Y+rnorm(100,sd=0.1),-0.5*Y+rnorm(100,sd=0.1),rnorm(100,sd=0.1))
#' cv <- RMRCE_cv(X,Y,c(0.001,0.01),c(1,5))
#' optlambda <- cv[which.min(cv[,3]),1]
#' optalpha <- cv[which.min(cv[,3]),2]
#' fit <- RMRCE(X,Y,optlambda,optalpha)

RMRCE_cv=function(X,Y,lambdavec,alphavec) {
      n=nrow(X)
      p=ncol(X)
      nk=round(n/5)
      folder=c(rep(nk,4),n-4*nk)
      para <- expand.grid(lambdavec,alphavec)
      error <- apply(para,1,function(i) {
            count=1
            testerror=0
            for(k in 1:5) {
                  Xk=X[-(count:(count+folder[k]-1)),]
                  Yk=Y[-(count:(count+folder[k]-1))]
                  betahatk=RMRCE(Xk,Yk,i[1],i[2])
                  X_t=X[count:(count+folder[k]-1),]
                  Y_t=Y[count:(count+folder[k]-1)]
                  Y_Y=matrix(rep(Y_t,folder[k]),nrow=folder[k])-t(matrix(rep(Y_t,folder[k]),nrow=folder[k]))
                  signYY=sign(Y_Y)
                  x=X_t%*%as.matrix(betahatk,ncol=1)
                  testerror=testerror+sum(signYY*sign(matrix(rep(x,folder[k]),nrow=folder[k])-t(matrix(rep(x,folder[k]),nrow=folder[k]))))/(folder[k]*(folder[k]-1))
                  count=count+folder[k] 
            }
            testerror/5      
      })
      data.frame(lambda=lambdavec,alpha=alphavec,testerror=error)
}
