
library(MASS)
library(matrixStats)

N=3  #size of the first sample
M=4#size of second sample
sigma0_2=1 #standard deviation of the second sample 
lnsigma=seq(from=-1,to=1,by=0.01)
sigma0_1=10^(lnsigma)
N_simulations=1000000 

#constants as defined by Aspin
lambda_i=c(1/N,1/M)
lambda_ir=matrix(data=NA,nrow=5,ncol=2)
for (i in 1:2)
  for (r in 1:5)
    lambda_ir[r,i]=lambda_i[i]^r
f_iu=matrix(data=NA,nrow=5,ncol=2)
  for (u in 1:5)
  {    f_iu[r,1]=1/(N-1)^u
  f_iu[r,2]=1/(M-1)^u
  }
# \V_{ru} as defined by Aspin
VRU<-function(r,u,s12,s22) (lambda_ir[r,1]*(s12^r)+lambda_ir[r,2]*(s22^r))/(lambda_i[1]*s12+lambda_i[2]*s22)^r

mu=0

dekninga3=matrix(data=rep(0,length(sigma0_1)*length(mu)), nrow=length(mu), ncol=length(sigma0_1)) #coverage of Aspin test, 3 summands
dekninga2=matrix(data=rep(0,length(sigma0_1)*length(mu)), nrow=length(mu), ncol=length(sigma0_1)) #coverage of Aspin test, 2 summands
for (i in 1:length(mu)) #mu could be longer
{
  for (j in 1:length(sigma0_1))
  { 
    #Simulate samples X og Y
    X=matrix(rnorm(N*N_simulations,mean=mu[i],sd=sigma0_1[j]),
             nrow=N,              # number of rows 
             ncol=N_simulations)  
    Y=matrix(rnorm(M*N_simulations,mean=0,sd=sigma0_2),
             nrow=M,              # number of rows 
             ncol=N_simulations)   

    s12=colVars(X) #not to calculate it many times
    s22=colVars(Y) #not to calculate it many times

    U=colMeans(X)-colMeans(Y) #not to calculate it many times

#3.3. Aspin
    dzeta=1.644854 # it is qnorm(0.95)
    
    #summands in notation of Aspin
    h0=dzeta*sqrt(lambda_i[1]*s12+lambda_i[2]*s22)
    h1=0.25*(1+dzeta^2)*VRU(2,1,s12,s22)*h0
    h2=(-0.5*(1+dzeta^2)*VRU(2,2,s12,s22)+1/3*(3+5*dzeta^2+dzeta^4)*VRU(3,2,s12,s22)-
         1/32*(15+32*dzeta^2+9*dzeta^4)*(VRU(2,1,s12,s22))^2)*h0
    h3=h0*((1+dzeta^2)*VRU(2,3,s12,s22)-2*(3+5*dzeta^2+dzeta^4)*VRU(3,3,s12,s22)+1/8*(15+32*dzeta^2+9*dzeta^4)*VRU(2,2,s12,s22)*VRU(2,1,s12,s22)+
      1/8*(75+173*dzeta^2+63*dzeta^4+5*dzeta^6)*VRU(4,3,s12,s22)-1/12*(105+298*dzeta^2+140*dzeta^4+15*dzeta^6)*VRU(3,2,s12,s22)*VRU(2,1,s12,s22)+
        1/384*(945+3169*dzeta^2+1811*dzeta^4+243*dzeta^6)*((VRU(2,1,s12,s22))^3 ) )
    h=h0+h1+h2+h3
    forcast_aspin3=sign(abs(U)-h)
    forcast_aspin2=sign(abs(U)-h0-h1-h2)
    dekninga3[i,j]=sum(forcast_aspin3>0)/N_simulations
    dekninga2[i,j]=sum(forcast_aspin2>0)/N_simulations    
  }
}

plot(lnsigma,dekninga3[1,],type='l',col='red',ylim=c(0.00,0.35), ylab="", xlab=expression(paste("log"[10],sigma[X],"/",sigma[Y])))
par ( new=TRUE)
plot(lnsigma,dekninga2[1,],type='l',ylim=c(0.00,0.35), ylab="", xlab="")
par ( new=TRUE)
plot(lnsigma,rep(0.05,length(sigma0_1)),type='l',col='darkgrey',lty='dashed',ylim=c(0.00,0.35), ylab="", xlab="")
legend("topleft",legend=c("2 ledd","3 ledd"),col=c('black', 'red'), title="Dekning for Aspin test",pch=1)
