
##read otu table
setwd("D:/Dropbox/warming diversity/warming network")
otutab<-read.csv("file:///D:/Dropbox/warming diversity/warming network/node/2014allwarm24_networkOTU_16s_ITS_UPARSE.txt",header = T,row.names=1,sep="\t")
comm<-t(otutab)
sp.ra<-colMeans(comm)/40000   #relative abundance of each species

cormatrix<-cor(comm,use="pairwise.complete.obs")

cormatrix2<-cormatrix*(abs(cormatrix)>0.80)  #only keep links above the cutoff point
cormatrix2[is.na(cormatrix2)]<-0
diag(cormatrix2)<-0    #no links for self-self    
sum(abs(cormatrix2)>0)/2  #this should be the number of links. It's 2117, not 1499. Do we use additional methods to remove links?

sum(colSums(abs(cormatrix2))>0)  #511 species have at least one linkage with others.

network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched

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

#rm.p.list=seq(0.05,0.2,by=0.05)
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
                 year=rep(2009,40),treat=rep("unwarmed",40))

#currentdat<-dat1
currentdat<-rbind(dat1,currentdat)

write.csv(currentdat,"node/simuresult/simuresult.csv")



##plot
library(ggplot2)

#weighted result
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


ggplot(currentdat[currentdat$weighted=="weighted" & currentdat$year>2009,], aes(x=Proportion.removed, y=remain.mean, group=year, color=year)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()+
  facet_wrap(~treat)


ggplot(currentdat[currentdat$weighted=="unweighted" & currentdat$year>2009,], aes(x=Proportion.removed, y=remain.mean, group=year, color=year)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()+
  facet_wrap(~treat)
