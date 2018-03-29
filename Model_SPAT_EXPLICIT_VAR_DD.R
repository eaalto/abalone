

library(logitnorm)
#library(gplots)

###•••••••••••••••••••••  Spatial Arrangement  •••••••••••••••••••••••••••••••••••••••••••••••

dispersal_percentage=99                               # Define percentage of larve that are retained between +- dispersal distance
dispersal_quantile=(1-dispersal_percentage/100)/2

length_block=100                                       # Lunghezza di un blocco in m
length_coast=15000                                     # Lunghezza della linea di costa in m
nblock=round(length_coast/length_block)                # Number of blocks
Areablock=length_block*500/10000                       # Area of each block in ha

AreaTot_ha=Areablock*nblock                            # Area totale in ha


riempitoro=function(j,devst){pnorm(length_block*(j+1/2),0,devst)-pnorm(length_block*(j-1/2),0,devst)}
dev_stand_disp=function(y){uniroot(function(x)qnorm(dispersal_quantile,mean=0,sd=x)+y,lower=0,upper=8000)$root}



###•••••••••••••••••••••  Demographic Settings  •••••••••••••••••••••••••••••••••••••••••••••••

MLS=155  # set Minimum Landing Size
limits=c(5,30,55,80,105,130,MLS,180,205,230) # Set Class Limits [25 mm width adjusted to 155 mm]
                                                 # Actually, the first size class should include 1 year individuals
                                                 # Initial size of larvae in Bardos et al. 2006 is 0.3 mm
                                                 # Viana et al. measures 2000 postlarvae (unknown age): 909 +- 130 um
                                                 # Leighton et al.1981 uses postlarvae (unknown age) of 1-1.5 mm
                                                 # Leighton 1974: postlarvae of 30-40 days are 1.7-2 mm SL
                                                 # In collectors at Isla Natividad, larvae at settlment were ~500 um
                                                 # Using Bardos(2005) and starting with larvae all of 500um, >90% is between 3.5 and 27 mm
                                                 # We could simplify and say that recruitment goes all to 1st size class
                                                                                                  
num_class=length(limits)-1                       # Calculate number of size class

mean_size=rep(0,num_class)                       # Mean class length
for(i in 1:num_class) 
{mean_size[i]=(limits[i+1]+limits[i])/2}          

num_fished_class=sum(mean_size>MLS)              # Calculate number of fished class


W <- function(L){2.24*10^(-5)*(L^(3.36))}        # Weight [g] at Length relationship da Shepherd 1998 [with shell]
                                                 # to be multiplied by 0.4 to obtain weight [g] withouth shell [From Natividad Data]

SM=round(0.5/(1+exp(-(mean_size-135.99)/30.20)),1)        # Size at sexual maturity, Rossetto et al.(2013), accounting for 1:1 Sex Ratio


load("Fitting_Bardos")                             # load growth parameters

estimate5_3

source("trans_probabilities.R")                  # load R script for growth prob.
P=Trans_matrix(limits)    # Growth matrix

#rm(H)
#rm(lambda)  

#load('P')

H=function(harvest_rate) {N <- matrix(0,num_class,num_class);  # Harvest matrix (survival from harvest)
		diag(N)=ifelse(mean_size>=MLS,1-harvest_rate,1);
		return(N)}


Reclut=function(tot_settlers) # tot_settlers is the number of arrived settlers    
{recruits=c(0.01*tot_settlers*(exp(-tot_settlers/10000000)),rep(0,num_class-1)); return(recruits)} 



mdet=function(peso){exp(0.6346126-0.317388*log(peso))} ### ALL NAT DATA

###••••••••••••••••••••• INITIALIZATION•••••••••••••••••••••••••••••••••••••••••••••••
##### INITITALIZE

Abb_tot=matrix(0,nrow=nreplic,ncol=tmax)                      # Initialize output: total abundance (ind)
Abb_SSBtot=matrix(0,nrow=nreplic,ncol=tmax)                   # Initialize output: total catch (ind)
Catt_tot_w=matrix(0,nrow=nreplic,ncol=tmax)                   # Initialize output: total catch in weight (tons)
S=matrix(rep(0),num_class,num_class)                 # Initialize survival matrix 

Dens_SEXPL=array(0,c(num_class,tmax,nblock))         # Initialize Dens_SEXPL 
C_SEXPL=array(0,c(num_fished_class,tmax,nblock))     # Initialize Catch_SEXPL 
egg_produced=rep(0,nblock)                           # Initialize produced eggs
egg_arrived=rep(0,nblock)                            # Initialize arrived eggs
Disp=matrix(0,nrow=nblock,ncol=nblock)               # Initialize dispersal matrix

t_presim=150                                            # Pre-simulation 
Dens_pre=array(0,c(num_class,t_presim,nreplic))         # Pre-simulation density
Catch_pre=array(0,c(num_fished_class,t_presim,nreplic))         # Pre-simulation catch


                     
Dens_pre[,1,]=rep(50,num_class)                         # arbitrary initial condition

larval_survival=0.00309

egg_per_ind =3772*W(mean_size)*SM

diag(S)=exp(-mdet(W(mean_size)))
Eff_pre=c(rep(0,100),rep(0.3,50))
for (i in 1:nreplic)
{for (tp in 1:(t_presim-1))      
	{
		
		Dens_pre[,tp+1,i]=(P%*%S%*%H(Eff_pre[tp]))%*%Dens_pre[,tp,i]+Reclut(sum(larval_survival*egg_per_ind*Dens_pre[,tp,i]))
		
	Catch_pre[,tp,i]=((diag(1,num_class)-H(Eff_pre[tp]))%*%Dens_pre[,tp,i])[(num_class-num_fished_class+1):num_class]
	Catch_pre[,t_presim,i]=((diag(1,num_class)-H(Eff_pre[t_presim]))%*%Dens_pre[,t_presim,i])[(num_class-num_fished_class+1):num_class]
	
	} }

#save(Dens_pre,file=paste(save_path,'/','Dens_pre',sep=''))
#save(Catch_pre,file=paste(save_path,'/','Catch_pre',sep=''))


#Dens_above100mm=apply(Dens_pre[5:9,,],c(2,3),sum)
#mean(Dens_above100mm[150,])

InVect=apply(Dens_pre,c(1,2),mean)[,t_presim]
#UnfishedVect=apply(Dens_pre,c(1,2),mean)[,100]



###••••••••••••••••••••• SIMULATION MODEL•••••••••••••••••••••••••••••••••••••••••••••••


dd_matrix=matrix(rgamma(tmax*nreplic,shape=3,rate=0.006),nrow=tmax,ncol=nreplic)



model=function(reservesizeb,percentage,harvest_rate)   # function that returns matrix of Abb_tot, Catt_tot, Catt_totw                                                          
                                                       # given the spatial arrangement
                                                       # arguments are:
                                                       # I.   reservesizeb (size of individual reserves [block])
                                                       # II.  protection level (from 1 to 100), 
                                                       # III. harvest rate (from 0 to 1)
{

nblockprot=round(percentage*nblock/100,0)     # Calculation of spatial effort
numres=round(nblockprot/reservesizeb,0)
nblockfish=nblock-nblockprot

x=floor(nblockfish/numres)
y=ceiling(nblockfish/numres)

nx=nblockfish-x*numres
ny=numres-nx


if(is.na(nx)|is.na(ny)) {Eff<-rep(harvest_rate,nblock)} else
if (nx==ny) {Eff<-rep(c(rep(0,reservesizeb),rep(harvest_rate,x),rep(0,reservesizeb),rep(harvest_rate,y)),ny)} else 
if (nx>ny){Eff<-c(rep(c(rep(0,reservesizeb),rep(harvest_rate,x),rep(0,reservesizeb),rep(harvest_rate,y)),ny),
	rep(c(rep(0,reservesizeb),rep(harvest_rate,y)),(nx-ny)))} else 
if (nx<ny) {Eff<-c(rep(c(rep(0,reservesizeb),rep(harvest_rate,x),rep(0,reservesizeb),rep(harvest_rate,y)),nx),rep(c(rep(0,reservesizeb),rep(harvest_rate,x)),(ny-nx)))}

                                
for(rep in 1:nreplic)  # Start loop on replicates
{
Dens_SEXPL[,1,]=InVect

for(t in 1:(tmax-1))    # Start loop on time
{
	devst=dev_stand_disp(dd_matrix[t,rep])
	                                                                
for(j in 1:(nblock/2))
{Disp[,j]=riempitoro(c((1-j):(nblock-j)),devst)+rev(riempitoro(c(j:(nblock+(j-1))),devst))}
Disp[1:nblock,(nblock/2+1):nblock]=Disp[nblock:1,(nblock/2):1]



egg_produced= apply(sweep(Dens_SEXPL[,t,],egg_per_ind,MARGIN=1,'*'),2,sum)    # vectorized 
egg_arrived=apply(sweep(Disp,MARGIN=2,egg_produced,'*'),1,sum)   # vectorized



for(block in 1:nblock)
{
Dens_SEXPL[,t+1,block]=((P%*%S%*%H(Eff[block])))%*%Dens_SEXPL[ ,t,block]+Reclut(larval_survival*egg_arrived[block])
C_SEXPL[,t,block]=((diag(1,num_class)-H(Eff[block]))%*%Dens_SEXPL[,t,block])[(num_class-num_fished_class+1):num_class]
C_SEXPL[,tmax,block]=((diag(1,num_class)-H(Eff[block]))%*%Dens_SEXPL[,tmax,block])[(num_class-num_fished_class+1):num_class]

}


}                   # End loop on time

C_SEXPL_W=sweep(C_SEXPL,W(mean_size)[(num_class-num_fished_class+1):num_class]*0.4/1000,MARGIN=1,'*')    # Weight of catch [kg]


Dens_per_block=apply(Dens_SEXPL,c(3,2),sum)
Abb_SSB=apply(Dens_SEXPL[,,]*SM,c(3,2),sum)
Catt_per_block_w=apply(C_SEXPL_W,c(3,2),sum)  

Abb_tot[rep,]=apply(Dens_per_block,2,sum)*Areablock
Abb_SSBtot[rep,]=apply(Abb_SSB,2,sum)*Areablock
Catt_tot_w[rep,]=apply(Catt_per_block_w,2,sum)*Areablock
}                   # End loop on replicates

return(list(Abb_tot=Abb_tot,Catt_tot_w=Catt_tot_w,Abb_SSB=Abb_SSBtot))}









