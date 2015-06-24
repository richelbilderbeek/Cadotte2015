##Script to transform and randomize trees and calculate log-lokelihood values.
##from paper in review at Functional Ecology
##Marc Cadotte 2014

library(picante)
library(geiger)

##read in tree
tr<-read.tree("EXP1_MLtree_codes.txt")

##remove unused tips
todrop<-c("SOPT","SYNA","ASUM","LICY")
tr<-drop.tip(tr,todrop)

##scale tree to 160 million years	
tr<-rescale(tr,"depth",160)

##planted data
dat<-read.csv("planted.csv")

##community matrix
species<-dat[2:15] #RJCB: original was `species<-dat[5:18]`
spp.mat<-species
rownames(spp.mat)<-as.character(dat$PLOT)
spp.mat<-as.matrix(spp.mat)

##plots with > 1 species
p2<-apply(spp.mat,1,sum,na.rm=TRUE)
spp.mat2<-spp.mat[p2>1,]
spp.mat2[is.na(spp.mat2)]<-0

##biomass
datt<-read.csv("anal_dat.csv")
dat2<-datt[,c(14,4)]
rownames(dat2)<-datt$Plot


##logLik for richness
lm.r<-glm(datt$biomass~datt$Real.rich)
LLR<-logLik(lm.r) #-436.1247



##tree rescale, delta 0.01-100
tr.list<-NULL
delta<-seq(0,100,0.1)
for (i in 1:length(delta)){
	tr.list[[i]]<-rescale(tr,"delta",delta[i])
}

datB<-dat2

#run linear models
for (i in 1:length(delta)){
	
	tmp<-pd(spp.mat2,tr.list[[i]])
	datB[,i+2]<-tmp$PD[match(rownames(datB),rownames(tmp))]
}

datB[is.na(datB)]<-0
out<-data.frame(Delta=delta,Slope=0,loglikes=0)

for (i in 1:length(delta)){
	
	tmp.lm<-glm(datB$biomass~datB[,i+2])
	out$Slope[i]<-coefficients(tmp.lm)[2]
	out$loglikes[i]<-logLik(tmp.lm)
	
}

#best model close to observed tree, delta = 1.1
out$Delta[out$loglikes==max(out$loglikes)]
ll_delta<-out$loglikes[out$loglikes==max(out$loglikes)]

par(mfrow=c(1,2))
plot(out$Delta,out$loglikes,type="l",xlab="Delta",ylab="Loglikelihood")
abline(v=1,lty="dashed") #line for observed tree

#line for richness
abline(h=LLR,lty="dotted")

###max ll = -434.4512, at delta = 1.1


##Kappa transformation, K = 0 - 2.1
tr.listK<-NULL
kappa<-seq(0,2.1,0.01)
for (i in 1:length(kappa)){
	tr.listK[[i]]<-rescale(tr,"kappa",kappa[i])
}

datK<-dat2

for (i in 1:length(kappa)){
	
	tmp<-pd(spp.mat2,tr.listK[[i]])
	datK[,i+2]<-tmp$PD[match(rownames(datK),rownames(tmp))]
}

datK[is.na(datK)]<-0

outK<-data.frame(Kappa=kappa,Slope=0,loglikes=0)

for (i in 1:length(kappa)){
	
	tmp.lm<-glm(datK$biomass~datK[,i+2])
	outK$Slope[i]<-coefficients(tmp.lm)[2]
	outK$loglikes[i]<-logLik(tmp.lm)
	
}

#best model K=0.41
outK$Kappa[outK$loglikes==max(outK$loglikes)]
ll_kappa<-outK$loglikes[outK$loglikes==max(outK$loglikes)]

##observed LL
obs.LL<-outK$loglikes[outK$Kappa==1] #-434.4566

plot(outK$Kappa,outK$loglikes,type="l",xlab="Kappa",ylab="Loglikelihood")
abline(v=1,lty="dashed")

#line for richness
abline(h=LLR,lty="dotted")

#####three randomizations
#first randomization: shuffle edge lengths
nrand<-1000

datR<-dat2

for (i in 1:nrand){
	tr.tmp<-tr
	tr.tmp$edge.length<-sample(tr$edge.length,length(tr$edge.length))
	tmp<-pd(spp.mat2,tr.tmp)
	datR[,i+2]<-tmp$PD[match(rownames(datR),rownames(tmp))]
}

datR[is.na(datR)]<-0

outR<-data.frame(Run=1:nrand,Slope=0,loglikes=0)

for (i in 1:nrand){
	
	tmp.lm<-glm(datR$biomass~datR[,i+2])
	outR$Slope[i]<-coefficients(tmp.lm)[2]
	outR$loglikes[i]<-logLik(tmp.lm)
	
}

mean(outR$loglikes) #-435.8876
1-prob<-rank(c(obs.LL,outR$loglikes))[1]/(nrand+1) #prob = obs. PD, P = 0.271
1-prob<-rank(c(ll_delta,outR$loglikes))[1]/(nrand+1) # prob = max delta, P = 0.255
1-prob<-rank(c(ll_kappa,outR$loglikes))[1]/(nrand+1) #prob = max Kappa, P = 0.185


##compare with swap tip.label
datT<-dat2

for (i in 1:nrand){
	tr.tmp<-tr
	tr.tmp$tip.label<-sample(tr$tip.label,length(tr$tip.label))
	tmp<-pd(spp.mat2,tr.tmp)
	datT[,i+2]<-tmp$PD[match(rownames(datT),rownames(tmp))]
}

datT[is.na(datT)]<-0

outT<-data.frame(Run=1:nrand,Slope=0,loglikes=0)

for (i in 1:nrand){
	
	tmp.lm<-glm(datT$biomass~datT[,i+2])
	outT$Slope[i]<-coefficients(tmp.lm)[2]
	outT$loglikes[i]<-logLik(tmp.lm)
	
}


mean(outT$loglikes) #-437.9607
prob<-rank(c(obs.LL,outT$loglikes))[1]/(nrand+1) #prob = obs. PD = 0.002
prob<-rank(c(ll_delta,outT$loglikes))[1]/(nrand+1) #prob = optimal delta P < 0.001
prob<-rank(c(ll_kappa,outT$loglikes))[1]/(nrand+1) #prob = optimal Kappa < 0.001

##randomize topology and tip labels, root to tip =160

datRT<-dat2

for (i in 1:nrand){
	tr.tmp<-rcoal(length(tr$tip.label))
	tr.tmp$tip.label<-sample(tr$tip.label)
	tr.tmp<-rescale(tr.tmp,"depth",160)
	tmp<-pd(spp.mat2,tr.tmp)
	datRT[,i+2]<-tmp$PD[match(rownames(datRT),rownames(tmp))]
}

datRT[is.na(datRT)]<-0

outRT<-data.frame(Run=1:nrand,Slope=0,loglikes=0)

for (i in 1:nrand){
	
	tmp.lm<-glm(datRT$biomass~datRT[,i+2])
	outRT$Slope[i]<-coefficients(tmp.lm)[2]
	outRT$loglikes[i]<-logLik(tmp.lm)
	
}
mean(outRT$loglikes) #-438.7721
prob<-rank(c(obs.LL,outRT$loglikes))[1]/(nrand+1) #prob = obs. PD = 0.005
prob<-rank(c(ll_delta,outRT$loglikes))[1]/(nrand+1) #prob = optimal Delta = <0.001
prob<-rank(c(ll_kappa,outRT$loglikes))[1]/(nrand+1) #prob = optimal Kappa = <0.001


##plot 3 histograms
par(mfrow=c(3,1))
hist(outR$loglikes,main=NULL,xlab="Loglikelihood",col="grey95")
#abline(v=obs.LL,lwd=2,col="violetred")
abline(v=ll_delta,lwd=4,col="wheat2")
abline(v=ll_kappa,lwd=4,col="skyblue2")

abline(v=LLR,lty="dotted",lwd=2)


hist(outT$loglikes,main=NULL,xlab="Loglikelihood",col="grey95", xlim=c(-442,-432))
abline(v=ll_delta,lwd=4,col="wheat2")
abline(v=ll_kappa,lwd=4,col="skyblue2")

abline(v=LLR,lty="dotted",lwd=2)

hist(outRT$loglikes,main=NULL,xlab="Loglikelihood",col="grey95")
abline(v=ll_delta,lwd=4,col="wheat2")
abline(v=ll_kappa,lwd=4,col="skyblue2")

abline(v=LLR,lty="dotted",lwd=2)