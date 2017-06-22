#' RMRCE
#' 
#' Regularized Maximum Rank Correlation Estimator
#'
#' This function fits RMRCE (Regularized Maximum Rank Correlation Estimator). It is recommended that the tuning parameters be tuned using cross validation (see RMRCE_cv). RMRCE may need some time to finish with large datasets.
#'
#' @param X A numeric matrix of explanatory variables.
#' @param Y A numeric vector of response.
#' @param lambda Tuning parameter lambda that can be tuned using cross validation. See RMRCE_cv.
#' @param alpha Tuning parameter alpha that can be tuned using cross validation. See RMRCE_cv. Default value is 5.
#' @return A numeric vector of regression coefficients.
#' @export
#' @author Fang Han, Hongkai Ji, Zhicheng Ji, Honglang Wang <zji4@@jhu.edu>
#' @examples
#' Y <- rnorm(100)
#' X <- cbind(Y+rnorm(100,sd=0.1),-0.5*Y+rnorm(100,sd=0.1),rnorm(100,sd=0.1))
#' fit <- RMRCE(X,Y,0.01,5)

RMRCE=function(X,Y,lambda,alpha=5) {
      uni_opt_obj=function(theta,index,lambda,alpha,beta,X,tmpYsign)
      {
            beta[index]=theta
            x=alpha*X%*%as.matrix(beta,ncol=1)
            tmpx <- matrix(rep(x,n),nrow=n)
            core=tmpYsign*(tmpx-t(tmpx))
            sum(pnorm(core))/(n*(n-1))+lambda*abs(theta)
      }
      n=nrow(X)
      p=ncol(X)
      tmpY <- matrix(rep(Y,n),nrow=n)
      Y_Y=tmpY-t(tmpY)
      signYY=sign(Y_Y)
      tmpYsign <- -sign(tmpY-t(tmpY))
      T=apply(X,2,function(x) sum(signYY*sign(matrix(rep(x,n),nrow=n)-t(matrix(rep(x,n),nrow=n))))/(n*(n-1)))
      j1=which.max(abs(T))
      betahat1=T/T[j1]
      betahat0=betahat1
      for(j in setdiff(1:p,j1)){
            result=optimize(uni_opt_obj,interval=c(-2,2),index=j,lambda=lambda,alpha=alpha,beta=betahat1,X=X,tmpYsign=tmpYsign)
            betahat1[j]=result[[1]]
      }
      count=1
      
      while(sqrt(sum((betahat1-betahat0)*(betahat1-betahat0)))>1e-5 & (count<=500))
      {
            count=count+1
            
            betahat0=betahat1
            for(j in setdiff(1:p,j1)){
                  result=optimize(uni_opt_obj,interval=c(-2,2),index=j,lambda=lambda,alpha=alpha,beta=betahat1,X=X,tmpYsign=tmpYsign)
                  betahat1[j]=result[[1]]
            }
      }
      betahat1
}
