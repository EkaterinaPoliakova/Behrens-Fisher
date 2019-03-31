
library(MASS)
library(matrixStats)
partitition=100  #of quantiles here

teta=pi/partitition*0.5*(1:partitition) # the angle in Behrens-Fisher distribution
freedom1=2
freedom2=3
#The quantiles for Behrens-Fisher distribution for 100 values of teta

quant0025=rep(0,partitition)
quant0975=rep(0,partitition)

# Step 2. Read quantiles for the Behrens-Fisher distribution from a file

quant0025=scan(file = "df 2 3.txt", what = quant0025, nmax = 100, skip = 10) #we read quantiles that were calculated before
quant0975=scan(file = "df 2 3.txt", what = quant0975, nmax = 100, skip = 40)
partitition_teta_to_use=9901 #we will actually use 99*100+1 values of teta
#there will be 99 more values of teta between each two neighbour tetas from before


q0025=rep(0,partitition_teta_to_use)

q0975=rep(0,partitition_teta_to_use)

#We find the new quantiles with linear interpolation:
for (i in 1:partitition_teta_to_use)
{ distance_to_lower_int=(i+99)/100-floor((i+99)/100)       
q0025[i]=quant0025[floor((i+99)/100)]*(1-distance_to_lower_int)+quant0025[ceiling((i+99)/100)]*distance_to_lower_int
q0975[i]=quant0975[floor((i+99)/100)]*(1-distance_to_lower_int)+quant0975[ceiling((i+99)/100)]*distance_to_lower_int
}
#Quantiles of Behrens-Fisher distribution are ready for use!
###########################################################################################################
#Calculating the table
BF=0 #This variable is of 7.5 Gb or should be deleted, othervise the workspace in R crashes very often
N=3  #size of the first sample
M=4#size of second sample
sigma0_2=1 #standard deviation of the second sample 
lnsigma=seq(from=-3,to=3,by=0.1)
sigma0_1=10^(lnsigma)
N_simulations=4000000 


#mu=c(0,0.1,0.2,0.3,0.4,0.5,0.7,1,2,4,10,20)
mu=0
dekning=matrix(data=rep(0,length(sigma0_1)*length(mu)), nrow=length(mu), ncol=length(sigma0_1))
dekning1=matrix(data=rep(0,length(sigma0_1)*length(mu)), nrow=length(mu), ncol=length(sigma0_1))

for (i in 1:length(mu))
{
  for (j in 1:length(sigma0_1))
  { 
    #Step 3.1 Simulate samples X og Y
    X=matrix(rnorm(N*N_simulations,mean=mu[i],sd=sigma0_1[j]),
             nrow=N,              # number of rows 
             ncol=N_simulations)  
    Y=matrix(rnorm(M*N_simulations,mean=0,sd=sigma0_2),
             nrow=M,              # number of rows 
             ncol=N_simulations)   
    #Step 3.2 Do Welch-Satterwaite test
    Teta=colVars(X)/colVars(Y)
    frihetsgrad=((Teta+N/M)^2)/((Teta^2)/(N-1)+((N/M)^2)/(M-1))  
    U=colMeans(X)-colMeans(Y)
    nevner=sqrt((colVars(X)/N+colVars(Y)/M))
    forcast_welch=sign(abs(U)-nevner*qt((1-0.025), df=frihetsgrad))
    dekning[i,j]=dekning[i,j]+sum(forcast_welch>0)/N_simulations
   
    #Step 3.2 Do Behrens-Fisher test
    #compute the angle for Behrens-Fisher distribution 
    Theta=atan(sqrt(colVars(X)*M/colVars(Y)/N))*2/pi+0.000101  #we avoid using "round"because it looks like ii
    dekning_bf=sum((U/nevner-q0025[Theta*partitition_teta_to_use])<0)+sum((U/nevner-q0975[Theta*partitition_teta_to_use])>0)
    dekning1[i,j]=dekning1[i,j]+dekning_bf/N_simulations
    #do the Behrens-Fisher tests
  }
}

write((dekning1), file="power BF df1=2, df2=3.txt", ncol=length(mu), sep="\t")
write((dekning), file="power W df1=2, df2=3.txt", ncol=length(mu),sep="\t")

plot(lnsigma,dekning1[1,],type='l',col='red',ylim=c(0.00,0.085), ylab="", xlab=expression(paste("log"[10],sigma[X],"/",sigma[Y])))
par ( new=TRUE)
plot(lnsigma,dekning[1,],type='l',ylim=c(0.00,0.085), ylab="", xlab="")
par ( new=TRUE)
plot(lnsigma,rep(0.05,length(sigma0_1)),type='l',col='grey',lty='dashed',ylim=c(0.00,0.085), ylab="", xlab="")
legend("bottomleft",legend=c("Welch-Tatterthwaite test","Behrens-Fisher test"),col=c('black', 'red'), title="Dekning",pch=1)
