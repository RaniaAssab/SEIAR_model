#Model Representation
#      ___________________________
#      |                          |
#      |                          |   
# ---->Sr---->Er---->Ir---->Ar---->Rr
#      |      |     |      |      |
#      |      |     |      |      |
#      |      |     |      |      |
#
#------------------------------------
#     ___________________________
#     |                          |
#     |                          |   
#    Ss---->Es---->Is---->As---->Rs


#Requires and libraries needed

require(stats)
require(graphics)
library(gdata)
require(ggplpot2)

#cmd line for utf-8 encoding : eval(parse("testTime.R",encoding="UTF-8"))

#####################################################
#########IMPORTATION OF THE DATA TO MODELISE#########
#####################################################

#Will be developed later
#datas<-read.csv("../../TrainingR/datas/datas.csv")
# => Recovery of the total number of people


#####################################################
#####################################################

#function to calculate one step of stochastic seir
SEIAR.onestep <- function (x, params,a,b,c,iso_code) {
	Sr <- x[2]
	Er <- x[3]
	Ir <- x[4]
	Ar <- x[5]
	Rr <- x[6]
	Sp <- x[7]
	Ep <- x[8]
	Ip <- x[9]
	Ap <- x[10]
	Rp <- x[11]
	Nr <- Sr+Er+Ir+Ar+Rr
	incP<-0
	incR<-0
	with(
		as.list(params),
		{
			#List of relative rates :
			e1 <- mu*Nr
			e2 <- e1+mu*Sr
			if (iso_code==1){
				e3 <- e2+Sr*((1-b)*betaRR*(Ir+alpha*Er+alpha*Ar)+(1-a)*betaPR*((1-eta)*Ip+alpha*Ep+alpha*Ap))
			}
			else if (iso_code==2){
				e3 <- e2+Sr*((1-b)*betaRR*(0*Ir+alpha*Er+alpha*Ar)+(1-a)*betaPR*((1-eta)*Ip+alpha*Ep+alpha*Ap))
			}
			e4 <- e3+mu*Er
			e5 <- e4+Er/dincub
			e6 <- e5+mu*Ir
			e7 <- e6+mu*Ar
			e8 <- e7+Ir/dsymp
			e9 <- e8+Ar/ro
			e10 <- e9+mu*Rr
			e11 <- e10+Rr/theta
			
			if (iso_code==1){
				e12 <- e11+Sp*((1-c)*betaPP*((1-eta)*Ip+alpha*Ep+alpha*Ap)+(1-a)*betaPR*(Ir+alpha*Er+alpha*Ar))
			}
			else if (iso_code==2){
				e12 <- e11+Sp*((1-c)*betaPP*((1-eta)*Ip+alpha*Ep+alpha*Ap)+(1-a)*betaPR*(0.4*Ir+alpha*Er+alpha*Ar))
			}

			
			e13 <- e12+Ep/dincub
			e14 <- e13+Ip/dsymp
			e15 <- e14+Ap/ro
			e16 <- e15+Rp/theta
			
			#calculate total rate :
			total.rate <- e16
			#inter event time :
			tau <- rexp(n=1,rate=total.rate)
			#uniform random deviate :
			U <- runif(1)
			#death of a recovered event id 'default'
			
			new.seiar <- c(Sr,Er,Ir,Ar,Rr,Sp+1,Ep,Ip,Ap,Rp-1)

			ct <- 0

			if (U<=(e1)/total.rate) { new.seiar<-c(Sr+1,Er,Ir,Ar,Rr,Sp,Ep,Ip,Ap,Rp)} 
			else if (U<=(e2)/total.rate) {new.seiar<-c(Sr-1,Er,Ir,Ar,Rr,Sp,Ep,Ip,Ap,Rp)}
			else if (U<=(e3)/total.rate) {new.seiar <- c(Sr-1,Er+1,Ir,Ar,Rr,Sp,Ep,Ip,Ap,Rp)}
			else if (U<=(e4)/total.rate) {new.seiar <- c(Sr,Er-1,Ir,Ar,Rr,Sp,Ep,Ip,Ap,Rp)}
			else if (U<=(e5)/total.rate) {
			                                new.seiar <- c(Sr,Er-1,Ir+1,Ar,Rr,Sp,Ep,Ip,Ap,Rp)
			                                incR<-1}
			else if (U<=(e6)/total.rate) {new.seiar <- c(Sr,Er,Ir-1,Ar,Rr,Sp,Ep,Ip,Ap,Rp)}
			else if (U<=(e7)/total.rate) {new.seiar <- c(Sr,Er,Ir,Ar-1,Rr,Sp,Ep,Ip,Ap,Rp)}
			else if (U<=(e8)/total.rate) {new.seiar <- c(Sr,Er,Ir-1,Ar+1,Rr,Sp,Ep,Ip,Ap,Rp)}
			else if (U<=(e9)/total.rate) {new.seiar <- c(Sr,Er,Ir,Ar-1,Rr+1,Sp,Ep,Ip,Ap,Rp)}
			else if (U<=(e10)/total.rate) {new.seiar <- c(Sr,Er,Ir,Ar,Rr-1,Sp,Ep,Ip,Ap,Rp)}
			else if (U<=(e11)/total.rate) {new.seiar <- c(Sr+1,Er,Ir,Ar,Rr-1,Sp,Ep,Ip,Ap,Rp)}
			
			else if (U<=(e12)/total.rate) { new.seiar <- c(Sr,Er,Ir,Ar,Rr,Sp-1,Ep+1,Ip,Ap,Rp)} 
			else if (U<=(e13)/total.rate) {
			                                new.seiar <- c(Sr,Er,Ir,Ar,Rr,Sp,Ep-1,Ip+1,Ap,Rp)
			                                incP<-1}
			else if (U<=(e14)/total.rate) {new.seiar <- c(Sr,Er,Ir,Ar,Rr,Sp,Ep,Ip-1,Ap+1,Rp)}
			else if (U<=(e15)/total.rate) {new.seiar <- c(Sr,Er,Ir,Ar,Rr,Sp,Ep,Ip,Ap-1,Rp+1)}

			ct <- ct+tau
			#store results
			c(tau,new.seiar,ct,incP,incR)
		}
		
		)
}


#function to simulate stochastic SEIR
SEIAR.model <- function(x, params, nstep, T,hpr,hrr,hpp,iso_code) {
	#set up array to store results
	output <- array(dim=c((T+2),14))
	#names variables :
	colnames(output) <- c("time",
						"Sr","Er","Ir","Ar","Rr",
						"Sp","Ep","Ip","Ap","Rp",
						"cum.time",
						"Incidence_p","Incidence_r")
	#first record of output is initial conditions
	output[1,] <- x
	var <- output[1,12] #var corresponde a cum.time
	n<-0 #jour simule en cours
	ligne<-1 #ligne dans le fichier de sortie
	inciR<-0 #incidence (par jour) chez les residents
	inciP<-0 #incidence (par jour) chez le personnel
	k<-0
	while (n<=T){ 
		x <- SEIAR.onestep(x,params,hpr,hrr,hpp,iso_code)
		var <- var+x[12]
		inciP<-inciP+x[13]
		inciR<-inciR+x[14]
		if(floor(var)>n){ #si au moins 1 jour entier s'est ecoule depuis la dernière fois qu'il s'est passe qq chose
			k<-floor(var)-n #nb entier de jours ecoules depuis la dernière fois qu'il s'est passe qq chose
			if (k>1){ #si des jours ont ete sautes
		  		for (j in 1:(k-1)){#on commence par noter ces jours où il ne s'est rien passe (dans la limite de T+2 lignes)
		    		ligne<-ligne+1
		    		if (ligne<(T+3))
 		    			{output[ligne,]<-output[ligne-1,]
 		    			output[ligne,12]<-n+j
 		    			output[ligne,13]<-0
 		    			output[ligne,14]<-0}
 		    		j<-j+1
		    	}
		  	}
  			ligne<-ligne+1 #puis on note ce qu'il s'est passe (dans la limite de T+2 lignes)
      		if (ligne<(T+3)){
      			output[ligne,]<-x
				output[ligne,12]<-floor(var)
				output[ligne,13]<-inciP
				output[ligne,14]<-inciR
			}
			n<-floor(var)
			inciP<-0
			inciR<-0
		}
	}

	#return output
	output
}


data_recuperation<-function(data){
	##############################################################
	####################################
	# Incidence
	####################################
	##############################################################

	Incidence.vector<-array(dim=c(nSimulations,(T+1)))
	for (i in 1:nSimulations){
	Incidence.vector[i,]<-data[[i]]$Incidence_r[1:(T+1)]
	}


	####################################
	# Moyenne
	####################################
	
	#mean
	CasCumR1<-c(nSimulations)
	for (i in 1:nSimulations){
		CasCumR1[i]<-sum(Incidence.vector[i,])
	}
	val1<-mean(CasCumR1)
	
	Mean.vector<-c()
	for (m in 1:T+1){
		Mean.vector[m]<-mean(Incidence.vector[,m])
	}
	

	

	####################################
	# Ecart type
	####################################

	Sd.vector<-c()
	#sd
	val2<-sd(CasCumR1)
	for (s in 1:length(Mean.vector)){
	Sd.vector[s]<-sd(sum(Incidence.vector[,m]))
	}
	Sd.vector[1]<-sd(Incidence.vector[,1])
	Sd.vector[2]<-sd(Incidence.vector[,2])

	####################################
	# Intervalle de confiance
	####################################

	IC_min.vector<-c()
	IC_max.vector<-c()
	for (ic in 1:length(Mean.vector)){
	IC_min.vector[ic]<-Mean.vector[ic]-1.96*Sd.vector[ic]/sqrt(nSimulations)
	IC_max.vector[ic]<-Mean.vector[ic]+1.96*Sd.vector[ic]/sqrt(nSimulations)
	}
	IC_min.vector[1]<-Mean.vector[1]-1.96*Sd.vector[1]/sqrt(nSimulations)
	IC_max.vector[1]<-Mean.vector[1]+1.96*Sd.vector[1]/sqrt(nSimulations)

	IC_min.vector[2]<-Mean.vector[2]-1.96*Sd.vector[2]/sqrt(nSimulations)
	IC_max.vector[2]<-Mean.vector[2]+1.96*Sd.vector[2]/sqrt(nSimulations)

	####################################
	# Temps
	####################################

	time_interval<-c(0:T)


	####################################
	# Tableau recapitulatif
	####################################

	Recap.data<-data.frame(Time_interval=time_interval,
					Residents=Mean.vector,
					IC_min=IC_min.vector,
					IC_max=IC_max.vector,
					rownames=c(0:T)
					)
	Recap.data
}

treatment<-function(data){
	vec<-c()
	#Incidence
	Incidence.vector<-array(dim=c(nSimulations,(T+1)))
	for (i in 1:nSimulations){
		Incidence.vector[i,]<-data[[i]]$Incidence_r[1:(T+1)]
	}
	#mean
	CasCumR<-c(nSimulations)
	for (i in 1:nSimulations){
		CasCumR[i]<-sum(Incidence.vector[i,])
	}
	val1<-mean(CasCumR)
	#sd
	val2<-sd(CasCumR)
	#Confidence interval
	Sub_IC_min.vector<-c()
	Sub_IC_max.vector<-c()
	for (ic in 1:nSimulations){
		Sub_IC_min.vector[ic]<-val1-1.96*val2/sqrt(nSimulations)
		Sub_IC_max.vector[ic]<-val1+1.96*val2/sqrt(nSimulations)
	}
	val3<-mean(Sub_IC_min.vector)
	val4<-mean(Sub_IC_max.vector)
	
	vec[1]<-val1
	vec[2]<-val2
	vec[3]<-val3
	vec[4]<-val4
	
	vec
}


casCumules1 <- function(T,h){
	adh_interval<-c(0:10)
	Mean.vector<-c()
	Sd.vector<-c()
	IC_min.vector<-c()
	IC_max.vector<-c()
	case<-1
	if (h==hpr){
		hpr<-0
		hP.vec <- c(1:10)
		while(hpr<=1){
			params$hpr<-hpr
			print(params$hpr)
			print("data generation")
			data <- data.generate(nSimulations,hpr,hrr,hpp,iso_code)
			print("treatment")
			vec_treat<-treatment(data)
			Mean.vector[case]<-vec_treat[1]
			Sd.vector[case]<-vec_treat[2]
			IC_min.vector[case]<-vec_treat[3]
			IC_max.vector[case]<-vec_treat[4]
			print("mean")
			print(Mean.vector)
			print("sd")
			print(Sd.vector)
			print("IC_min")
			print(IC_min.vector)
			print("IC_max")
			print(IC_max.vector)
			hP.vec[case]<-hpr
			case<-case+1
			hpr<-hpr+0.1		
			
		}
	}
	else if (h==hrr){
		hrr<-0
		hP.vec <- c(1:10)
		while(hrr<=1){
			params$hrr<-hrr
			print(params$hrr)
			data <- data.generate(nSimulations,hpr,hrr,hpp,iso_code)
			vec_treat<-treatment(data)
			Mean.vector[case]<-vec_treat[1]
			Sd.vector[case]<-vec_treat[2]
			IC_min.vector[case]<-vec_treat[3]
			IC_max.vector[case]<-vec_treat[4]
			print("mean")
			print(Mean.vector)
			print("sd")
			print(Sd.vector)
			print("IC_min")
			print(IC_min.vector)
			print("IC_max")
			print(IC_max.vector)
			hP.vec[case]<-hrr
			case<-case+1
			hrr<-hrr+0.1		
			
		}
	}
	else if (h==hpp){
		hpp<-0
		hP.vec <- c(1:10)
		while(hpp<=1){
			params$hpp<-hpp
			print(params$hpp)
			data <- data.generate(nSimulations,hpr,hrr,hpp,iso_code)
			vec_treat<-treatment(data)
			Mean.vector[case]<-vec_treat[1]
			Sd.vector[case]<-vec_treat[2]
			IC_min.vector[case]<-vec_treat[3]
			IC_max.vector[case]<-vec_treat[4]
			print("mean")
			print(Mean.vector)
			print("sd")
			print(Sd.vector)
			print("IC_min")
			print(IC_min.vector)
			print("IC_max")
			print(IC_max.vector)
			hP.vec[case]<-hpp
			case<-case+1
			hpp<-hpp+0.1		
			
		}
	}
	
	CC.matrix<-data.frame(adh=c(0,10,20,30,40,50,60,70,80,90,100),nb_Cas_r=Mean.vector, IC_min_r=IC_min.vector, IC_max_r=IC_max.vector)
	CC.matrix
}

createFile <- function(mat){
	write.table(mat, "~/out.csv", sep = "|", quote = FALSE, col.names = NA)
}

maximum234 <- function(data){
	max.y <- 0
	for( i in 1:length(data[[1]]$Ir)){
		if (is.na(data[[1]]$Ir[i]) == FALSE){
			max.y <- max(max.y,data[[1]]$Ir[i])
		}
	}
	max.y
}


data.generate <- function(nSimulations,hpr,hrr,hpp,iso_code){
	datas <- vector(mode='list',length=nSimulations)
	#simulate nSimulations times :
	for (k in 1:nSimulations){
		datas[[k]] <- as.data.frame(SEIAR.model(xstart,params,nstep,T,hpr,hrr,hpp,iso_code))
	}
	datas
}



##################
#######MAIN#######
##################


#set seed pour conserver les memes valeurs aleatoires
set.seed(38499583)

####RESIDENTS####
pop.size_res <- 100
#pop.size_res <- 89
#Initial number of exposed
Ezero_res <- 0
#Initial number of infected
Izero_res <- 1
#Initial number of asymptomatic infected
Azero_res <- 0
#initial realistic number of susceptibles
Szero_res <- pop.size_res-Ezero_res-Izero_res-Azero_res


####STAFF####
pop.size_per <- 63
#Initial number of exposed
Ezero_per <- 0
#Initial number of infected
Izero_per <- 0
#Initial number of asymptomatic infected
Azero_per <- 0
#initial realistic number of susceptibles
Szero_per <- pop.size_per-Ezero_per-Izero_per-Azero_per

#number of event to simulate
nstep <- 5000
#Time interval
T<-100

 
xstart <- c(time=0,
			Sr=Szero_res,
			Er=Ezero_res,
			Ir=Izero_res,
			Ar=Azero_res,
			Rr=pop.size_res-Szero_res-Izero_res-Ezero_res-Azero_res,
			Sp=Szero_per,
			Ep=Ezero_per,
			Ip=Izero_per,
			Ap=Azero_per,
			Rp=pop.size_per-Szero_per-Izero_per-Ezero_per-Azero_per,
			ct=0, #cum.time
			ip=0, #incidence moyenne
			ir=0
			)
nSimulations<-4000
print("Pending...")
 	
#parameters :
pTrans<-0.07 #0.069
params=list(mu=0.02/30,
    		betaRR=0.025*pTrans,betaPP=0.045*pTrans,betaPR=0.1*pTrans,
		alpha=10/100,
		dincub=1,
		dsymp=2,
		ro=10,
		theta=5*365,
		eta=68/100
		)
hpp<-0.00
hrr<-0.10
hpr<-0.15
data1 <- data.generate(nSimulations,hpr,hrr,hpp,iso_code=1)
data2 <- data.generate(nSimulations,hpr,hrr,hpp,iso_code=2)

max.time <- T


Recup_dat_Inc1<-data_recuperation(data1)
#Recup_dat_Inc1<-treatment(data1)
#Recup_dat_Inc2<-data_recuperation(data2)
Recup_dat_Inc2<-treatment(data2)

####################################
# Display
####################################

qplot(x=Time_interval, color=Residents,
      y=Residents, ymax=IC_max, ymin=IC_min,
      geom = c("path", "pointrange"), data=Recup_dat_Inc1)


##############################################################
####################################
# Impact hygiene des mains
####################################
##############################################################

####################################
# Recuperation informations HM
####################################

CC.m1<-casCumules1(T,h=hpr)
CC.m2<-casCumules1(T,h=hrr)
CC.m3<-casCumules1(T,h=hpp)

####################################
# Creation vecteurs mean et IC
####################################

hpr.ICmin<-CC.m1$IC_min_r
hpr.mean<-CC.m1$nb_Cas_r
hpr.ICmax<-CC.m1$IC_max_r

hrr.ICmin<-CC.m2$IC_min_r
hrr.mean<-CC.m2$nb_Cas_r
hrr.ICmax<-CC.m2$IC_max_r

hpp.ICmin<-CC.m3$IC_min_r
hpp.mean<-CC.m3$nb_Cas_r
hpp.ICmax<-CC.m3$IC_max_r


####################################
# Creation tableau mean et IC
####################################

matr<-data.frame(hpr_ICmin=hpr.ICmin,
			Residents_hpr=hpr.mean,
			hpr_ICmax=hpr.ICmax,
			hrr_ICmin=hrr.ICmin,
			Residents_hrr=hrr.mean,
			hrr_ICmax=hrr.ICmax,
			hpp_ICmin=hpp.ICmin,
			Residents_hpp=hpp.mean,
			hpp_ICmax=hpp.ICmax,
			Adherence_HH=CC.m1$adh)
				
createFile(matr)
####################################
# Display
####################################


#hpr
qplot(x=Adherence_HH, color=Residents_hpr,
      y=Residents_hpr, ymax=hpr_ICmax, ymin=hpr_ICmin, ylim=c(0,30),
      geom = c("path", "pointrange"), data=matr)

#hrr
qplot(x=Adherence_HH, color=Residents_hrr,
      y=Residents_hrr, ymax=hrr_ICmax, ymin=hrr_ICmin, ylim=c(0,30),
      geom = c("path", "pointrange"), data=matr)

#hpp
qplot(x=Adherence_HH, color=Residents_hpp,
      y=Residents_hpp, ymax=hpp_ICmax, ymin=hpp_ICmin, ylim=c(0,30),
      geom = c("path", "pointrange"), data=matr)


