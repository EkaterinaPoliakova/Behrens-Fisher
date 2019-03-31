#First we compute quantiles of Behrens-Fisher distribution
n_simulation=10000000  # to get critical values exact enough
# when there were fewer, there was no correspondance with the table
partitition=100  #of quantiles here

teta=pi/partitition*0.5*(1:partitition) # the angle in Behrens-Fisher distribution
###############
###############
freedom1=2     #here can be any degrees of freedom
freedom2=3
###############
################
#The quantiles for Behrens-Fisher distribution for 100 values of teta
quant0005=rep(0,partitition) #this is a 0.005 - quantile.
quant005=rep(0,partitition)
quant0025=rep(0,partitition)
quant095=rep(0,partitition)
quant0995=rep(0,partitition)
quant0975=rep(0,partitition)

num_runs=10  #it is impossible with more simulations because the large vectors
#can't be stored in the workspace and the "for"-cycles are too slow
#so we estimate the quantiles 10 times and find the average
#at some simulations there was 20 but it takes long time

for (one_more_counter in 1:num_runs)
{
  T1=rt(n_simulation,freedom1) #we generate two student distributed vectors
  T2=rt(n_simulation,freedom2)
  BF= T2 %*% t(cos(teta)) - T1 %*% t(sin(teta)) #a Behrens- Fisher distributed variable
  #unfortunately the size for BF is limited but working with matrises increases the  speed
  
  for (i in 1:partitition)  #calculate quantiles. It works slowly
  {
    quant0005[i]=quantile(BF[,i],0.005)/num_runs+quant0005[i]
    quant005[i]=quantile(BF[,i],0.05)/num_runs+quant005[i]
    quant0025[i]=quantile(BF[,i],0.025)/num_runs+quant0025[i]
    quant095[i]=quantile(BF[,i],0.95)/num_runs+quant095[i]
    quant0995[i]=quantile(BF[,i],0.995)/num_runs+quant0995[i]
    quant0975[i]=quantile(BF[,i],0.975)/num_runs+quant0975[i]
  }
  
}


#Quantiles of Behrens-Fisher distribution are ready for use!
###################################################################################

write(quant005, file=paste("df ",freedom1," ",freedom2,".txt"), append=TRUE,  ncolumns=10, sep="\t")
write(quant0025, file=paste("df ",freedom1," ",freedom2,".txt"), append=TRUE,  ncolumns=10, sep="\t")
write(quant0005, file=paste("df ",freedom1," ",freedom2,".txt"), append=TRUE,  ncolumns=10, sep="\t")
write(quant095, file=paste("df ",freedom1," ",freedom2,".txt"), append=TRUE,  ncolumns=10, sep="\t")
write(quant0975, file=paste("df ",freedom1," ",freedom2,".txt"), append=TRUE,  ncolumns=10, sep="\t")
write(quant0995, file=paste("df ",freedom1," ",freedom2,".txt"), append=TRUE,  ncolumns=10, sep="\t")