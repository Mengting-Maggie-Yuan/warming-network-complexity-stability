# contributor: Daliang Ning

source("/Users/maggieyuan/Documents/!annual_network/GitHub/Fig3_and_S7.stability/iieggrtools.r")

wd="/Users/maggieyuan/Documents/!annual_network/GitHub/Fig3_and_S7.stability/"

com.file="input_file/OTUtable_NetworkedOTUs_AllSamples.txt"
treat.file="input_file/SampleMap_AllSamples.txt"

#####################
wd=iwd(wd)

comm=t(lazyopen(com.file))
treat=lazyopen(treat.file)
sum(is.na(comm)) # check NA
comm[is.na(comm)]=0# if have, should change to zero

sampc=match.name(rn.list = list(comm=comm,treat=treat))
comm=sampc$comm
treat=sampc$treat

head(treat) # check which column is plot ID, which column is year.
plot.lev=unique(treat$Plot_ID)
year.lev=sort(unique(treat$Year))

zeta.lev=2:length(year.lev)

year.windows=lapply(1:length(zeta.lev),
                  function(i)
                  {
                    zetai=zeta.lev[i]
                    lapply(1:(length(year.lev)-zetai+1),function(j){year.lev[j:(j+zetai-1)]})
                  })
names(year.windows)=zeta.lev
year.windows

# function of community stability among a group of samples
comstab<-function(subcom){((nrow(subcom)*sum(apply(subcom,2,min)))/sum(subcom))^0.5}

#######

stabl=lapply(1:length(year.windows),
             function(i)
             {
               stabi=t(sapply(1:length(plot.lev),
                              function(j)
                              {
                                plotj=plot.lev[j]
                                sapply(1:length(year.windows[[i]]),
                                       function(k)
                                       {
                                         yearwdk=year.windows[[i]][[k]]
                                         sampijk=rownames(treat)[which((treat$Plot_ID==plotj) & (treat$Year %in% yearwdk))]
                                         outijk=NA
                                         if(length(sampijk) < length(yearwdk))
                                         {
                                           warning("plot ",plotj," has missing year in year window ",paste(yearwdk,collapse = ","))
                                         }else if(length(sampijk) > length(yearwdk)){
                                           warning("plot ",plotj," has duplicate samples in at least one year of window ",paste(yearwdk,collapse = ","))
                                         }else{
                                           comijk=comm[which(rownames(comm) %in% sampijk),,drop=FALSE]
                                           outijk=comstab(comijk)
                                         }
                                         outijk
                                       })
                              }))
               if(nrow(stabi)!=length(plot.lev) & nrow(stabi)==1){stabi=t(stabi)}
               rownames(stabi)=plot.lev
               colnames(stabi)=sapply(year.windows[[i]],function(v){paste0("Zeta",zeta.lev[i],paste0(v,collapse = ""))})
               stabi
             }) # Y09 W does not have matches, will be warnings.
stabm=Reduce(cbind,stabl)

treat.rm09=treat[which(treat$Year!="Y09"),,drop=FALSE] # just to annotate the treatment information of each plot, remove the special year
output=data.frame(treat.rm09[match(rownames(stabm),treat.rm09$Plot_ID),1:5,drop=FALSE],stabm,stringsAsFactors = FALSE)
output=output[order(output$Warming,output$Clipping,output$Precip,output$Plot_ID),,drop=FALSE]
head(output)
#save.file(output,filename = "MultiOrderStability.csv")

# plot Figure S7abcde - multiorder compositional stability. Adapted from the concept of Zeta diveristy in Hui and McGeoch American Naturalist 2014.
cs.mo = output

cs.long <- gather(cs.mo, year, cs, Zeta2Y09Y10:Zeta6Y09Y10Y11Y12Y13Y14, factor_key=TRUE)
cs.long <- cs.long[-which(is.na(cs.long$cs)),]
cs.long$EndYear = as.numeric(gsub(".*Y", "", cs.long$year))

get_text <- function(y, x){
  lm = lm(y~x)
  lm_p = round(summary(lm)$coefficients[2,4],3)
  lm_r = round(summary(lm)$adj.r.squared,3)
  lm_s = round(summary(lm)$coefficients[2,1],3)
  txt = paste(lm_s, " (", lm_r, ", ", lm_p, ")", sep="")
  return(txt)
}

plot_cs_full <- function(xc, yc, xw, yw){
  #  yl=c(min(yc, yw), max(yc, yw))

  plot(yw~xw, ylim=c(0.2,1), xlim=c(10,14), ylab="Compositional stability", xlab="Year", pch=19, cex.axis=0.7, cex.lab=0.7, tck=-0.03, col="#e7211f")
  abline(lm(yw~xw), col="#e7211f")

  points(yc~xc, col="#214da0")
  abline(lm(yc~xc), col="#214da0")

  txt1 = get_text(yw,xw)
  txt2 = get_text(yc,xc)
  mtext(txt1, side=3, line=0.8, cex=0.5, col="#e7211f")
  mtext(txt2, side=3, line=0.1, cex=0.5, col="#214da0")
}

plot_cs_zeta6 <- function(xc, yc, xw, yw){
  #  yl=c(min(yc, yw), max(yc, yw))

  plot(yw~xw, ylim=c(0.2,1), xlim=c(10,14), ylab="Compositional stability", xlab="Year", pch=19, cex.axis=0.7, cex.lab=0.7, tck=-0.03, col="#e7211f")
  points(yc~xc, col="#214da0")

  wt = wilcox.test(yw, yc)
  mtext(paste(wt$statistic, wt$p.value, sep=","), cex=0.5)
}

par(mfrow=c(2,3))
plot_cs_full(xc=cs.long$EndYear[intersect(which(cs.long$Warming=="N"), grep("Zeta2", cs.long$year))],
             yc=cs.long$cs[intersect(which(cs.long$Warming=="N"), grep("Zeta2", cs.long$year))],
             xw=cs.long$EndYear[intersect(which(cs.long$Warming=="W"), grep("Zeta2", cs.long$year))],
             yw=cs.long$cs[intersect(which(cs.long$Warming=="W"), grep("Zeta2", cs.long$year))])

plot_cs_full(xc=cs.long$EndYear[intersect(which(cs.long$Warming=="N"), grep("Zeta3", cs.long$year))],
             yc=cs.long$cs[intersect(which(cs.long$Warming=="N"), grep("Zeta3", cs.long$year))],
             xw=cs.long$EndYear[intersect(which(cs.long$Warming=="W"), grep("Zeta3", cs.long$year))],
             yw=cs.long$cs[intersect(which(cs.long$Warming=="W"), grep("Zeta3", cs.long$year))])

plot_cs_full(xc=cs.long$EndYear[intersect(which(cs.long$Warming=="N"), grep("Zeta4", cs.long$year))],
             yc=cs.long$cs[intersect(which(cs.long$Warming=="N"), grep("Zeta4", cs.long$year))],
             xw=cs.long$EndYear[intersect(which(cs.long$Warming=="W"), grep("Zeta4", cs.long$year))],
             yw=cs.long$cs[intersect(which(cs.long$Warming=="W"), grep("Zeta4", cs.long$year))])

plot_cs_full(xc=cs.long$EndYear[intersect(which(cs.long$Warming=="N"), grep("Zeta5", cs.long$year))],
             yc=cs.long$cs[intersect(which(cs.long$Warming=="N"), grep("Zeta5", cs.long$year))],
             xw=cs.long$EndYear[intersect(which(cs.long$Warming=="W"), grep("Zeta5", cs.long$year))],
             yw=cs.long$cs[intersect(which(cs.long$Warming=="W"), grep("Zeta5", cs.long$year))])

plot_cs_zeta6(xc=cs.long$EndYear[intersect(which(cs.long$Warming=="N"), grep("Zeta6", cs.long$year))],
              yc=cs.long$cs[intersect(which(cs.long$Warming=="N"), grep("Zeta6", cs.long$year))],
              xw=cs.long$EndYear[intersect(which(cs.long$Warming=="W"), grep("Zeta6", cs.long$year))],
              yw=cs.long$cs[intersect(which(cs.long$Warming=="W"), grep("Zeta6", cs.long$year))])


# plot Figure 3d. Compositional stability of adjacent two years

colors = c(rep("#214da0",5), rep("#e7211f", 5))
syms = c(1,1,1,1,1, 19,19,19,19,19)

plot_cs <- function(yw, yc, ebw, ebc){
  x=c(1,2,3,4,5)
  yl=c(min(c(yw-ebw,yc-ebc)), max(c(yw+ebw,yc+ebc)))

  plot(yw~x, ylim=yl, ylab="Compositional stability", xlab="Year of warming", pch=19, cex.axis=0.7, cex.lab=0.7, tck=-0.03, col="red")
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

zeta2 = cs.long[grep("Zeta2", cs.long$year),]
zeta2_means = aggregate(zeta2$cs, list(zeta2$Warming, zeta2$EndYear), FUN=mean)
zeta2_sds = aggregate(zeta2$cs, list(zeta2$Warming, zeta2$EndYear), FUN=sd)

plot_cs(yw=zeta2_means$x[which(zeta2_means$Group.1=="W")],
        yc=zeta2_means$x[which(zeta2_means$Group.1=="N")],
        ebw=zeta2_sds$x[which(zeta2_sds$Group.1=="W")],
        ebc=zeta2_sds$x[which(zeta2_sds$Group.1=="N")])
