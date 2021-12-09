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
#input network matrix, number of removed keystone species, keystonespecies list,  and ra of all species
#return the proportion of species remained

#get the keystone species list
node.attri<-read.csv("input_files/NodeAttribute_Y14_W.txt",header = T, sep="\t")
module.hub<-as.character(node.attri$Name[node.attri$Zi > 2.5 & node.attri$Pi <= 0.62])

#consider cascade effects: removed species will further influence the remaining nodes

rand.remov2.once<-function(netRaw, rm.num, keystonelist, sp.ra, abundance.weighted=T){
  rm.num2<-ifelse(rm.num > length(keystonelist), length(keystonelist), rm.num)
  id.rm<-sample(keystonelist, rm.num2)
  net.Raw=netRaw #don't want change netRaw
  
  net.new=net.Raw[!names(sp.ra) %in% id.rm, !names(sp.ra) %in% id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[!names(sp.ra) %in% id.rm]
    
    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }
    
    sp.meanInteration<-colMeans(net.stength)
    
    
    while ( length(sp.meanInteration)>1 & min(sp.meanInteration) <=0){
      id.remain<- which(sp.meanInteration>0) 
      net.new=net.new[id.remain,id.remain]
      sp.ra.new=sp.ra.new[id.remain]
      
      if (abundance.weighted){
        net.stength= net.new*sp.ra.new
      } else {
        net.stength= net.new
      }
      
      if (length(net.stength)>1){
        sp.meanInteration<-colMeans(net.stength)
      } else{
        sp.meanInteration<-0
      }
      
    }
    
    remain.percent<-length(sp.ra.new)/length(sp.ra)
    
    remain.percent}
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=length(module.hub)),
                 year=rep(2014,2*length(module.hub)),treat=rep("warmed",2*length(module.hub)))

currentdat = dat1
write.csv(currentdat,"targeted_deletion_results_Y14_W.csv")

##plot
library(ggplot2)

currentdat$Number.hub.removed

ggplot(currentdat[currentdat$weighted=="weighted",], aes(x=Number.hub.removed, y=remain.mean, group=treat, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("red","blue"))+
  xlab("Number of module hubs removed")+
  ylab("Proportion of species remained")+
  theme_light()+
  facet_wrap(~year, ncol=3,scales="free_x")


ggplot(currentdat[currentdat$weighted=="unweighted",], aes(x=Number.hub.removed, y=remain.mean, group=treat, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Number of module hubs removed")+
  ylab("Proportion of species remained")+
  theme_light()+
  facet_wrap(~year, ncol=3,scales="free_x")
