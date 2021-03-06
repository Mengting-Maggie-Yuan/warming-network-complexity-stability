# contributor: Linwei Wu

setwd("/Users/maggieyuan/Documents/!annual_network/GitHub/Fig3_and_S7.stability/") # working directory

otutab = read.table("input_file/OTUtable_NetworkedOTUs_AllSamples.txt", sep="\t", header=T, row.names=1)
treatused = read.table("input_file/SampleMap_AllSamples.txt", header=T, sep="\t", row.names=1)

###calculate persistence for each plot
##here persistence= proportion of species persist in  certain number of years
comm<-t(otutab)
sum(row.names(comm)==row.names(treatused))
fake09.N<-comm[1:24,];row.names(fake09.N)<-gsub("S","N",row.names(fake09.N))
comm2<-rbind(fake09.N,comm)

fake09.treat<-treatused[1:24,];row.names(fake09.treat)<-gsub("S","N",row.names(fake09.treat))
fake09.treat$Plot_full_name<-gsub("S","N",fake09.treat$Plot_full_name)
fake09.treat$Plot_ID<-gsub("S","N",fake09.treat$Plot_ID)
treat2<-rbind(fake09.treat,treatused)
treat2$Plot_ID=factor(treat2$Plot_ID)
treat2$Warming[1:48]<-treat2$Warming[49:96][match(treat2$Plot_ID[1:48], treat2$Plot_ID[49:96])] # organize warming treatment for Year 09.

test1<-split(data.frame(comm2),treat2$Plot_ID)

years=c("Y09","Y10","Y11","Y12","Y13","Y14")
l <- rep(list(0:1), 6)
allcombines<-expand.grid(l)
allcombines<-allcombines[c(4,7,13,25,49, 8,15,29,57, 16,31,61, 32,63, 64),] # only select adjacent years
allcombines<-t(allcombines)

test<-sapply(test1, function(x){
  test3<-sapply(1:ncol(allcombines), 
                function(i){
    coms<-allcombines[,i]*x
    total.sp<-sum(colSums(coms>0)>0)
    persist.sp<-sum(colSums(coms>0)==sum(allcombines[,i]))
    prop=persist.sp/total.sp
    
    names(prop)=paste0(years[allcombines[,i]>0],collapse="")
    prop
  })
  test3
})

persist.result<-t(test)
output=data.frame(treat2[match(rownames(persist.result),treat2$Plot_ID),1:5,drop=FALSE],persist.result,stringsAsFactors = FALSE)
names(output)<-paste(c(rep("",5),
                       rep("Zeta2",5),
                       rep("Zeta3",4),
                       rep("Zeta4",3),
                       rep("Zeta5",2),
                       rep("Zeta6",1)), names(output), sep="")
head(output)
#save.file(output,filename = "MultiOrderPersistence.csv")

## plot Figure S7fghij - multiorder persistence. Adapted from the concept of Zeta diveristy in Hui and McGeoch American Naturalist 2014.
op.mo = output

op.long <- gather(op.mo, year, op, Zeta2Y09Y10:Zeta6Y09Y10Y11Y12Y13Y14, factor_key=TRUE)
op.long$EndYear = as.numeric(gsub(".*Y", "", op.long$year))

get_text <- function(y, x){
  lm = lm(y~x) 
  lm_p = round(summary(lm)$coefficients[2,4],3)
  lm_r = round(summary(lm)$adj.r.squared,3)
  lm_s = round(summary(lm)$coefficients[2,1],3)	
  txt = paste(lm_s, " (", lm_r, ", ", lm_p, ")", sep="")	
  return(txt)
}

plot_op_full <- function(xc, yc, xw, yw){
  #  yl=c(min(yc, yw), max(yc, yw))
  
  plot(yw~xw, ylim=c(0,0.3), xlim=c(10,14), ylab="Node persistence", xlab="Year", pch=19, cex.axis=0.7, cex.lab=0.7, tck=-0.03, col="#e7211f")
  abline(lm(yw~xw), col="#e7211f")
  
  points(yc~xc, col="#214da0")
  abline(lm(yc~xc), col="#214da0")
  
  txt1 = get_text(yw,xw)
  txt2 = get_text(yc,xc)
  mtext(txt1, side=3, line=0.8, cex=0.5, col="#e7211f")
  mtext(txt2, side=3, line=0.1, cex=0.5, col="#214da0")
}

plot_op_zeta6 <- function(xc, yc, xw, yw){
  # yl=c(min(yc, yw), max(yc, yw))
  
  plot(yw~xw, ylim=c(0,0.3), xlim=c(10,14), ylab="Node persistence", xlab="Year", pch=19, cex.axis=0.7, cex.lab=0.7, tck=-0.03, col="#e7211f")
  points(yc~xc, col="#214da0")
  
  wt = wilcox.test(yw, yc)
  mtext(paste(wt$statistic, wt$p.value, sep=","), cex=0.5)
}

par(mfrow=c(2,3))
plot_op_full(xc=op.long$EndYear[intersect(which(op.long$Warming=="N"), grep("Zeta2", op.long$year))], 
             yc=op.long$op[intersect(which(op.long$Warming=="N"), grep("Zeta2", op.long$year))], 
             xw=op.long$EndYear[intersect(which(op.long$Warming=="W"), grep("Zeta2", op.long$year))], 
             yw=op.long$op[intersect(which(op.long$Warming=="W"), grep("Zeta2", op.long$year))])

plot_op_full(xc=op.long$EndYear[intersect(which(op.long$Warming=="N"), grep("Zeta3", op.long$year))], 
             yc=op.long$op[intersect(which(op.long$Warming=="N"), grep("Zeta3", op.long$year))], 
             xw=op.long$EndYear[intersect(which(op.long$Warming=="W"), grep("Zeta3", op.long$year))], 
             yw=op.long$op[intersect(which(op.long$Warming=="W"), grep("Zeta3", op.long$year))])

plot_op_full(xc=op.long$EndYear[intersect(which(op.long$Warming=="N"), grep("Zeta4", op.long$year))], 
             yc=op.long$op[intersect(which(op.long$Warming=="N"), grep("Zeta4", op.long$year))], 
             xw=op.long$EndYear[intersect(which(op.long$Warming=="W"), grep("Zeta4", op.long$year))], 
             yw=op.long$op[intersect(which(op.long$Warming=="W"), grep("Zeta4", op.long$year))])

plot_op_full(xc=op.long$EndYear[intersect(which(op.long$Warming=="N"), grep("Zeta5", op.long$year))], 
             yc=op.long$op[intersect(which(op.long$Warming=="N"), grep("Zeta5", op.long$year))], 
             xw=op.long$EndYear[intersect(which(op.long$Warming=="W"), grep("Zeta5", op.long$year))], 
             yw=op.long$op[intersect(which(op.long$Warming=="W"), grep("Zeta5", op.long$year))])

plot_op_zeta6(xc=op.long$EndYear[intersect(which(op.long$Warming=="N"), grep("Zeta6", op.long$year))], 
              yc=op.long$op[intersect(which(op.long$Warming=="N"), grep("Zeta6", op.long$year))], 
              xw=op.long$EndYear[intersect(which(op.long$Warming=="W"), grep("Zeta6", op.long$year))], 
              yw=op.long$op[intersect(which(op.long$Warming=="W"), grep("Zeta6", op.long$year))])


# Observed persistence of adjacent two years

colors = c(rep("#214da0",5), rep("#e7211f", 5))
syms = c(1,1,1,1,1, 19,19,19,19,19)

plot_op <- function(yw, yc, ebw, ebc){
  x=c(1,2,3,4,5)
  yl=c(min(c(yw-ebw,yc-ebc)), max(c(yw+ebw,yc+ebc)))
  
  plot(yw~x, ylim=yl, ylab="Observed Node Persistence", xlab="Year of warming", pch=19, cex.axis=0.7, cex.lab=0.7, tck=-0.03, col="red")
  abline(lm(yw~x), col="red")
  arrows(x, yw-ebw, x, yw+ebw, length=0.05, angle=90, code=3, col="red")
  
  points(yc~x, col="blue")
  abline(lm(yc~x), lty=2, col="blue")
  arrows(x, yc-ebc, x, yc+ebc, length=0.05, angle=90, code=3, col="blue")
  
  txt1 = get_text(yw,x)
  txt2 = get_text(yc,x)
  mtext(txt1, side=3, line=0.8, cex=0.5, col="red")
  mtext(txt2, side=3, line=0.1, cex=0.5, col="blue")
}


zeta2 = op.long[grep("Zeta2", op.long$year),]
zeta2_means = aggregate(zeta2$op, list(zeta2$Warming, zeta2$EndYear), FUN=mean)
zeta2_sds = aggregate(zeta2$op, list(zeta2$Warming, zeta2$EndYear), FUN=sd)

plot_op(yw=zeta2_means$x[which(zeta2_means$Group.1=="W")], 
        yc=zeta2_means$x[which(zeta2_means$Group.1=="N")], 
        ebw=zeta2_sds$x[which(zeta2_sds$Group.1=="W")], 
        ebc=zeta2_sds$x[which(zeta2_sds$Group.1=="N")])

# Figure 3d. combine above data with compositional stablity
