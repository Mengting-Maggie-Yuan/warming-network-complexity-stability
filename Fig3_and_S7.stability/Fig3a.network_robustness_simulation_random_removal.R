# contributor: Linwei Wu

setwd("/Users/maggieyuan/Documents/!annual_network/GitHub/Fig3_and_S7.stability/")

##read otu table
otutab<-read.table("input_file/OTUtable_AllOTUs_Y14_W.txt",header = T,row.names=1,sep="\t")
otutab[is.na(otutab)]<-0
##keep 12 otus
counts<-rowSums(otutab>0)
otutab<-otutab[counts>=12,]

comm<-t(otutab)
sp.ra<-colMeans(comm)/30000   #relative abundance of each species

###### there are two choices to get the correlation matrix #######
###### choice 1 (slow): calculate correlation matrix from OTU table
cormatrix=matrix(0,ncol(comm),ncol(comm))

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i]>0,comm[k,i],ifelse(comm[k,j]>0,0.01,NA))
    })
    speciesj<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j]>0,comm[k,j],ifelse(comm[k,i]>0,0.01,NA))
    })
    corij<-cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormatrix[i,j]<-cormatrix[j,i]<-corij
    
  }}

###### choice 2: directely read in correlation matrix downloaded from MENAP. MENAP downloaded correlation matrix is an upper triangle. Need to make it into a symetric matrix.
cormatrix.input  <- matrix(0, 2567, 2567)
cormatrix.input[row(cormatrix.input) >= col(cormatrix.input)] <- scan("input_files/CorrelationMatrix_Y14_W.txt")
cormatrix <- t(cormatrix.input)
for (i in 1:nrow(cormatrix)){
  for (j in 1:ncol(cormatrix)){
    if (i>j){cormatrix[i,j]<-cormatrix[j,i]}
  }
}
###### end of the two choices of correlation matrix ########

row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm) # if processed using MENAP, OTU order should match in the original OTU table and the correlation matrix downloaded from MENAP.

cormatrix2<-cormatrix*(abs(cormatrix)>=0.80)  #only keep links above the cutoff point
cormatrix2[is.na(cormatrix2)]<-0
diag(cormatrix2)<-0    #no links for self-self    
sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.

network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched

## robustness simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
 
  sp.meanInteration<-colMeans(net.stength)

  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 year=rep(2014,40),treat=rep("Warmed",40))

currentdat<-dat1

write.csv(currentdat,"random_removal_result_Y14_W.csv")


##plot
library(ggplot2)

currentdat$year

ggplot(currentdat[currentdat$weighted=="weighted",], aes(x=Proportion.removed, y=remain.mean, group=treat, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()+
  facet_wrap(~year, ncol=3)
  

ggplot(currentdat[currentdat$weighted=="unweighted",], aes(x=Proportion.removed, y=remain.mean, group=treat, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()+
  facet_wrap(~year, ncol=3)
