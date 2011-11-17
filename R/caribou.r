###########################################################################
#   Description:
#	The function abundance applies the methodology found in Rivest, L.-P.,  
#	Couturier, S. & Crépeau, H. (1998) for estimating caribou abundance 
#   using postcalving aggregations detected by radio telemetry
# 
#   Usage:
#	abundance(mat, n, model=c("H","I","T"), B)
# 
#   Arguments:
#
# 	mat : matrix containning in the first column the number of radio-collared animals 
#           in the detected (photographed) groups and in the second column the
#           corresponding size of the detected  groups 
#  	
#  	n : total number of active collars during the census 
#
#  	model: model to determine the probability that a group with collared animals is detected pi;
#	    H: homogeneity model pi = r , provides xi>1
#	    I: independence model pi =  1-(r^xi)	
# 	    T: threshold model pi = 1 if Xi>=B
#					    r if  1<=Xi<B
#
#  	B : bound for the threshold model 	
#
#   Values:
#
#	mp : the number of detected groups having radio-collared animals
#		
# 	xt : the total number of radio-collared animals found in the detected groups 
#
#	gnt : the total number of animals counted in the detected groups
#
# 	rr :  estimated parameter related to the probability of detection
#
#	se_rr: estimated standard error of rr
#    
#	mat_pi: matrix containing xi,gni and pi the detection probability i
#
#	T.hat : estimator for the total number of animals in a herd
# 
#	se_T.hat : estimated standard error of T.hat
#
#	loglikelihood : the maximum value of the loglikelihood function for the detected model
#
#	randomness_test: vector with the statistic and the p-value of a score test for the randomness assumption
#		 available only for the homogeneity model, the independance model and the threshold model for B=2 or 3
# 
#   call : function call
#
#
#   Source : Rivest, L.-P., Couturier, S. & Crépeau, H. (1998). Statistical Methods for estimating 
#            caribou abundance using postcalving aggregations detected by radio telemetry. Biometrics, 54 (3), 865-876
#
#
#  program realized by Hélène Crépeau 
#  Last update: 21-02-2011
#
############################################################################

abundance<-function(mat, n, model=c("H","I","T"), B)
{ 
  ### Validation des arguments
  mat <- as.matrix(mat)
  if (dim(mat)[2]!=2) stop("'mat' must have two columns")
  #if (any((mat%%1)!=0)||any(mat<0)) stop("'mat' must contain non-negative integers")
  if (!is.numeric(mat)||any(mat<0)) stop("'mat' must contain non-negative numerics")
  if (any(mat[,2]<mat[,1])) stop("each number in the second column of 'mat' must be greater or equal to the number in the first column on the same row") 
  
  valid.one(n,"numeric")
  if ((n%%1)!=0||n<=0) stop("'n' must be a positive integer")
  
  model <- model[1]
  if(!model%in%c("H","I","T")) stop("'model' can only take one of these values: 'H', 'I' or 'T'") 
  
  if(!missing(B)) {
    valid.one(B,"numeric")
    if ((B%%1)!=0||B<2) stop("'B' must be an integer greater or equal to 2")
  } else { if(model=="T") stop("Argument 'B' must be given with model='T'") }
  ############ 
  
  
# xi : vector of the number of radio-collared in the detected groups
  xi=mat[,1]
# gni : vector of the size of detected  groups 
  gni=mat[,2]
# mp : the number of detected groups having radio-collared animals
  mp=length(xi) 
# xt : the total number of radio-collared animals found in the detected groups 
  xt=sum(xi)
# gnt : the total number of animals counted in the detected groups
  gnt=sum(gni)
  
# 1. Estimation of r and its variance (Section 4   Rivest et al. 1998) 
  
  if (model=="H"){
    rr=sum(xi)/n
    pi<-rep(rr,mp)
    pip=rep(1,mp)
  }
  
  if (model=="I") {
    fr = function(r){sum(xi/(1-r^xi))- n}
    # Find the zero of a function <=> find the minimum in absolute value 
    abs_fr = function(r){abs(fr(r))}
    rr = optimize(abs_fr ,interval=c(0,1))$minimum
    pi=1-rr^xi
    pip=-xi*rr^(xi-1)
  }
  
  if (model=="T") {
    rr=((n-sum(xi[xi>(B-1)]))^(-1))*sum(xi[xi<B])
    pi=rr*(xi<B) + 1*(xi>=B)
    pip=1*(xi<B) + 0*(xi>=B)
  }
  
# variance estimation
  f1=sum((xi*pip/(pi^2)))^-2
  f2=sum(xi^2*(1-pi)/(pi^2))
  v2rr<-f1*f2
  ser=sqrt(v2rr)
  
# 2. Herd size estimation T.hat (equation 4  Rivest et al. 1998)
  
  min_t=sum(gni)
  max_t=n*max(gni)
  ft = function(t){t - sum(gni/(pi*(1-(1-gni/t)^n)))}
  abs_ft = function(t){abs(ft(t))}
  T.hat = optimize(abs_ft ,interval=c(min_t,max_t))$minimum
  T.hat
  
  
# 3. Variance estimation of T.hat (equation 6  Rivest et al. 1998)
  
  pi.i=1-(1-gni/T.hat)^n
  T0.n=sum(pip*gni/((pi^2)*pi.i))
  X0.d=sum(pip*xi/(pi^2))
  
# Detection probability for each group 
  mat_pi=cbind(xi,gni,pi,pi.i)
  mat_pi=mat_pi[order(mat_pi[,"gni"]),]
  
  terme1=1-(n/(T.hat^2))*sum((gni^2*(1-gni/T.hat)^(n-1))/(pi*pi.i^2))
  terme2=sum((gni^2*(1-pi.i))/(pi*pi.i^2))
  terme4=sum(((1-pi)/pi^2)*(gni/pi.i-xi*T0.n/X0.d)^2)
  
# computation of the third term (equation 1 et 2 Rivest et al. 1998)
  gamij=pi.ij=matrix(0,mp,mp) 
  for (i in 1:mp){
    for (j in 1:mp){
      pi.ij[i,j]=1-(1-gni[i]/T.hat)^n-(1-gni[j]/T.hat)^n+(1-(gni[i]+gni[j])/T.hat)^n
      gamij[i,j]<-pi.ij[i,j]-pi.i[i]*pi.i[j]
    } 
  }
  term3=matrix(0,mp,mp)
  for (i in 1:mp-1){
    for (j in (i+1):mp)
      term3[i,j]=(gni[i]*gni[j]*gamij[i,j])/(pi[i]*pi[j]*pi.i[i]*pi.i[j]*pi.ij[i,j])
  }
  terme3=2*sum(term3)
  
# Variance and SE of the herd size estimation
  
  var_T.hat=(terme2+terme3 +terme4)/(terme1^2)
  se_T.hat=sqrt(var_T.hat)
  
  
  
# 4. loglikelihood estimation (vraissemblance)  ( section 5.2 Rivest et al. 1998)
  
  if (model=="H"){
    vraissemblance=function(t){
      lambda=n*gni/t
      sum(xi*log(lambda)-log(exp(lambda)-1))}
  }
  
  if (model=="I") {
    vraissemblance = function(t){
      lambda=n*gni/t ;
      sum(xi*log(lambda) + log(1 - rr^xi)-log(exp(lambda) - exp(rr*lambda)))}
  }
  
  if (model=="T") {
    vraissemblance = function(t){
      rr=((n-sum(xi[xi>(B-1)]))^-1)*sum(xi[xi<B])
      lambda=n*gni/t
      somBm1=rep(0,mp)
      for(i in 1:(B-1)){
        somBm1=somBm1+(lambda^i)/factorial(i)
      }
      sum(xi*log(lambda)+log(1-(1-rr)*pip)-log(exp(lambda)-1-(1-rr)*somBm1))
    }
  }
  
  lik_T.hat=vraissemblance(T.hat)
  
# 5. Score test for the randomness assumption (section 5.1 Rivest et al. 1998)
  
#for the Homogeneity and the independence Model 
  
  if (model=="H" || model=="I") {
    if (model=="H") {r=0}
    if (model=="I") {r=rr}
    
# estimate theta 
    ftheta = function(t)
    {sum(xi) - 
          sum(exp(t+log(n) + log(gni))
                  + ((1-r)*exp(t + log(n) + log(gni))/
                    (exp((1-r)*exp(t + log(n) + log(gni))) - 1 )))
    };
    min_theta=-log(max(T.hat-2*se_T.hat,n));
    max_theta=-log(T.hat+2*se_T.hat);
    abs_ftheta = function(t){abs(ftheta(t))}
    theta.hat = optimize(abs_ftheta ,interval=c(min_theta,max_theta))$minimum
    theta.hat
    
# define the function g(tdi)
    g_tdi=expression(r*exp(tdi)+log(exp((1-r)*exp(tdi))-1))
    
# Define the derivative of the function
    gp_tdi=D(g_tdi,"tdi")
    gpp_tdi=D(gp_tdi,"tdi")
    gppp_tdi=D(gpp_tdi,"tdi")
    gpppp_tdi=D(gppp_tdi,"tdi")
    
    tdi<-theta.hat+log(n)+log(gni)
    
    vv2<-sum(eval(gpppp_tdi)+2*(eval(gpp_tdi))^2)-(sum(eval(gppp_tdi)))^2/(sum(eval(gpp_tdi)))
    
# Z statistic value  which is compared to normal distribution
    
    zobs<-sum((xi-eval(gp_tdi))^2-eval(gpp_tdi))/sqrt(vv2)
    pvalue=1-pnorm(zobs)
  }
  if (model=="T") {
    if (B==2)	{  
# estimate theta 
      r=rr
      ftheta = function(t)
      {sum(xi) - 
            sum((exp(exp(t+log(n)+log(gni))) * exp(t+log(n)+log(gni)) - (1 - r) * exp(t+log(n)+log(gni)))/(exp(exp(t+log(n)+log(gni))) - 
                      1 - (1 - r) * exp(t+log(n)+log(gni))))
      } 
      min_theta=-log(max(T.hat-2*se_T.hat,n));
      max_theta=-log(T.hat+2*se_T.hat);
      abs_ftheta = function(t){abs(ftheta(t))}
#      theta.hat = optimize(abs_ftheta ,interval=c(min_theta,max_theta))$minimum
      theta.hat = optimize(abs_ftheta ,interval=c(min_theta,max_theta))$minimum
      theta.hat
      
# define the function g(tdi)
      g_tdi=expression(log(exp(exp(tdi))-1-(1-r)*exp(tdi)))
      
# Define the derivative of the function
      gp_tdi=D(g_tdi,"tdi")
      gpp_tdi=D(gp_tdi,"tdi")
      gppp_tdi=D(gpp_tdi,"tdi")
      gpppp_tdi=D(gppp_tdi,"tdi")
      
      tdi<-theta.hat+log(n)+log(gni)
      
      vv2<-sum(eval(gpppp_tdi)+2*(eval(gpp_tdi))^2)-(sum(eval(gppp_tdi)))^2/(sum(eval(gpp_tdi)))
      
# Z statistic value  which is compared to normal distribution
      
      zobs<-sum((xi-eval(gp_tdi))^2-eval(gpp_tdi))/sqrt(vv2)
      pvalue=1-pnorm(zobs)
    }
    if (B==3) 	{
# estimate theta 
      r=rr
      ftheta = function(t)
      {sum(xi) - 
            sum((exp(exp(t+log(n)+log(gni)))*exp(t+log(n)+log(gni))-
                      (1-r)*(exp(t+log(n)+log(gni))+2*(exp(t+log(n)+log(gni))*
                          exp(t+log(n)+log(gni)))/2))/(exp(exp(t+log(n)+log(gni)))
                      -1-(1 - r)*(exp(t+log(n)+log(gni))+(exp(t+log(n)+log(gni))^2)/2)))
      } 
      min_theta=-log(max(T.hat-2*se_T.hat,n));
      max_theta=-log(T.hat+2*se_T.hat);
      abs_ftheta = function(t){abs(ftheta(t))}
      theta.hat = optimize(abs_ftheta ,interval=c(min_theta,max_theta))$minimum
      theta.hat
      
# define the function g(tdi)
      g_tdi=expression(log(exp(exp(tdi))-1-(1-r)*(exp(tdi)+(exp(tdi)^2)/2)))
      
# Define the derivative of the function
      gp_tdi=D(g_tdi,"tdi")
      gpp_tdi=D(gp_tdi,"tdi")
      gppp_tdi=D(gpp_tdi,"tdi")
      gpppp_tdi=D(gppp_tdi,"tdi")
      
      tdi<-theta.hat+log(n)+log(gni)
      
      vv2<-sum(eval(gpppp_tdi)+2*(eval(gpp_tdi))^2)-(sum(eval(gppp_tdi)))^2/(sum(eval(gpp_tdi)))
      
# Z statistic value  which is compared to normal distribution
      
      zobs<-sum((xi-eval(gp_tdi))^2-eval(gpp_tdi))/sqrt(vv2)
      pvalue=1-pnorm(zobs)
    }
    if (B>3)	{
      zobs=NA
      pvalue=NA
    }
  }		
  score=c(z_obs=zobs,pvalue=pvalue)
  
  out=list(mp=mp,xt=xt,gnt=gnt,rr=rr,se_rr=ser,mat_pi=mat_pi,T.hat=T.hat,se_T.hat=se_T.hat,loglikelihood=lik_T.hat,randomness_test=score,call=match.call())
  class(out)<- "abundance"
  return(out)
  
}

"print.abundance" <- function(x, ...) {
  cat("Given parameters:")
  cat("\nn: Total number of active collars during the census =",eval(x$call$n))
  nmodel <- if(is.null(x$call$model)) "Homogeneity" else { if(x$call$model=="H") "Homogeneity" 
            else { if(x$call$model=="I") "Independence" else "Threshold" }}
  cat("\nmodel:",nmodel,"model") 
  if(x$call$model=="T") cat(" with a bound of",x$call$B)
  
  cat("\n\nDescriptive statistics:")
  cat("\nmp: Number of detected groups having radio-collared animals =",x$mp)
  cat("\nxt: Total number of radio-collared animals in detected groups =",x$xt)
  cat("\ngnt: Total number of animals counted in detected groups =",x$gnt)
  
  cat("\n\nEstimation:\n")
  tableau2 <- cbind(estimate=round(x$rr,6),se=round(x$se_rr,6))
  rownames(tableau2) <- c("Parameter related to probability of detection:")
  print.default(tableau2, print.gap = 2, quote = FALSE, right=TRUE, ...)
  
#     cat("\n")
  tableau3 <- cbind(estimate=round(x$T.hat,0),se=round(x$se_T.hat,1))
  rownames(tableau3) <- c("Total number of animals in a herd:")
  print.default(tableau3, print.gap = 2, quote = FALSE, right=TRUE, ...)
  
  invisible(x)
}


###########################################################################
#   Description:
#	The function petersen estimates the total population size based
#   on the Lincoln-Petersen Index as applied to radio-telemetry data by White and Garrott (1990).
#	We use the Petersen estimator with Chapman's bias correction and the bias corrected 
#	standard error estimator of of Seber and Wittes.
# 
#	WHITE, G.C., and GARROTT, R.A. 1990. Analysis of wildlife radio-tracking data. San Diego: Academic Press.
#	
#
#   Usage:
#	petersen(mat,M,S)
# 
#   Arguments:
#
# 	mat : matrix containing in the first column the number of radio-collared 
#           in the detected (photographed) groups and in the second column the
#           corresponding size of the detected  groups 
#  	
#  	M : total number of active collars during the census 
#
#  	S : the minimum size that define well aggregated group
#
#
#
#   Values:
#
#  	M : total number of active collars during the survey
#
#  	S : the minimum size that define well aggregated group
#
#	G : Number of well define aggregated groups
#
#  	R : total number of radio-collared animal observed in the well aggregated groups  

#  	C : total number of animals observed in the well aggregated groups containing at least one radio-collared animal during the survey

#	LP_T.hat : modified lincoln-petersen estimator for the total number of animals in a herd
# 
#	ser_LP : estimated standard error of LP_T.hat

#   mat_aggre : A matrix containing a subset of the input matrix mat: the data for the well aggregated groups only.

#   call : function call
#
#
#  program realized by Hélène Crépeau
#  Last update: 21-02-2011
#
############################################################################

petersen=function(mat,M,S)
{ 
  ### Validation des arguments
  mat <- as.matrix(mat)
  if (dim(mat)[2]!=2) stop("'mat' must have two columns")
#  if (any((mat%%1)!=0)||any(mat<0)) stop("'mat' must contain non-negative integers")
  if (!is.numeric(mat)||any(mat<0)) stop("'mat' must contain non-negative numerics")
  if (any(mat[,2]<mat[,1])) stop("each number in the second column of 'mat' must be greater or equal to the number in the first column on the same row") 
  
  valid.one(M,"numeric")
  
  valid.one(S,"numeric")
  ############
  
# gni : vector of the size of detected  groups 
  gni=mat[,2]
  
  mat_aggre=mat[gni>S,]
  G=dim(mat_aggre)[1]
  R=sum(mat_aggre[,1])
  C=sum(mat_aggre[,2])
  
  if(M<R) warning("")
  
  LP_T.hat=((M+1)*(C+1))/(R+1)-1
  var_LP=((M+1)*(C+1)*(M-R)*(C-R))/((R+1)^2*(R+2))
  se_LP=sqrt(var_LP)
  
  out<-list(G=G,R=R,C=C,LP_T.hat=LP_T.hat,se_LP=se_LP,mat_aggre=mat_aggre,call=match.call())
  class(out)<- "petersen"    
  return(out)
  
}


"print.petersen" <- function(x, ...) {
  cat("Given parameters:")
  cat("\nM: Total number of active collars during the census =",eval(x$call$M))
  cat("\nS: Minimum size of well define aggregated groups =",eval(x$call$S)) 
  
  cat("\n\nDescriptive statistics on well aggregated groups:")
  cat("\nG: Number of well define aggregated groups =",x$G) 
  cat("\nR: Total number of radio-collared animal observed in the well aggregated groups =",x$R)
  cat("\nC: Total number of animals observed in the well aggregated groups =",x$C)
  
  cat("\n\nEstimation:\n")
  tableau3 <- cbind(estimate=round(x$LP_T.hat,0),se=round(x$se_LP,1))
  rownames(tableau3) <- c("Total number of animals in a herd:")
  print.default(tableau3, print.gap = 2, quote = FALSE, right=TRUE, ...)
  
  cat("\nThis function calculates a Lincoln-Petersen estimator with Chapman's bias correction and the bias corrected standard error estimator of Seber and Wittes\n")
  
  invisible(x)
}

"valid.one" <- function(x,type)
{
  if(!eval(call(paste("is.",type,sep=""),x))||length(x)!=1) {
    error <- paste("'",deparse(substitute(x)),"' must be a length-one object of type ",type,sep="")
    stop(error)
  } 
}
