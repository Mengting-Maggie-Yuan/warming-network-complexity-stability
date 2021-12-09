# The node overlap calculation follows the concept of Zeta diversity in Hui and McGeoch American Naturalist 2014.
# contributor: Maggie Yuan

library(ggplot2)

setwd("/Users/maggieyuan/Documents/!annual_network/GitHub/ExtendedDataFig5.constancy/")

tb = read.table("input_file/NetworkNodeTable_AllNetworks.txt", header=T, sep="\t")
head(tb)

nwlist = unique(tb$year_trt)
clist = nwlist[c(1,2,4,6,8,10)]
wlist = nwlist[c(1,3,5,7,9,11)]

c2 = combn(c(1:6), 2)
c3 = combn(c(1:6), 3)
c4 = combn(c(1:6), 4)
c5 = combn(c(1:6), 5)
c6 = combn(c(1:6), 6)

overlap <- function(id, tb, list){
  otu_overlap = tb$Name[which(tb$year_trt==list[id[1]])]
  for (i in 2:length(id)){
    otu = tb$Name[which(tb$year_trt==list[id[i]])]
    otu_overlap = intersect(otu_overlap, otu)
  }
  return(length(otu_overlap))
}

# control
overlap2 = apply(c2, 2, FUN=overlap, tb=tb, list=clist)
overlap3 = apply(c3, 2, FUN=overlap, tb=tb, list=clist)
overlap4 = apply(c4, 2, FUN=overlap, tb=tb, list=clist)
overlap5 = apply(c5, 2, FUN=overlap, tb=tb, list=clist)
overlap6 = apply(c6, 2, FUN=overlap, tb=tb, list=clist)

dist2 = apply(c2, 2, FUN=dist)
dist3 = colMeans(apply(c3, 2, FUN=dist))
dist4 = colMeans(apply(c4, 2, FUN=dist))
dist5 = colMeans(apply(c5, 2, FUN=dist))
dist6 = colMeans(apply(c6, 2, FUN=dist))

ol_c_full = data.frame(
  overlap = c(overlap2, overlap3, overlap4, overlap5, overlap6),
  gap = c(dist2, dist3, dist4, dist5, dist6),
  zeta = rep(c(2,3,4,5,6), times=c(length(overlap2), length(overlap3), length(overlap4), length(overlap5), length(overlap6)))
)
ol_c_full$warming = rep("control", nrow(ol_c_full))

mean2 = aggregate(overlap2, by=list(dist2), FUN=mean)
mean3 = aggregate(overlap3, by=list(dist3), FUN=mean)
mean4 = aggregate(overlap4, by=list(dist4), FUN=mean)
mean5 = aggregate(overlap5, by=list(dist5), FUN=mean)
mean6 = aggregate(overlap6, by=list(dist6), FUN=mean)

ol_c = rbind(mean2, mean3, mean4, mean5, mean6)
ol_c$zeta = rep(c(2,3,4,5,6), times=c(nrow(mean2), nrow(mean3), nrow(mean4), nrow(mean5), nrow(mean6)))
ol_c$warming = rep("control", nrow(ol_c))

# warming
overlap2 = apply(c2, 2, FUN=overlap, tb=tb, list=wlist)
overlap3 = apply(c3, 2, FUN=overlap, tb=tb, list=wlist)
overlap4 = apply(c4, 2, FUN=overlap, tb=tb, list=wlist)
overlap5 = apply(c5, 2, FUN=overlap, tb=tb, list=wlist)
overlap6 = apply(c6, 2, FUN=overlap, tb=tb, list=wlist)

dist2 = apply(c2, 2, FUN=dist)
dist3 = colMeans(apply(c3, 2, FUN=dist))
dist4 = colMeans(apply(c4, 2, FUN=dist))
dist5 = colMeans(apply(c5, 2, FUN=dist))
dist6 = colMeans(apply(c6, 2, FUN=dist))

ol_w_full = data.frame(
  overlap = c(overlap2, overlap3, overlap4, overlap5, overlap6),
  gap = c(dist2, dist3, dist4, dist5, dist6),
  zeta = rep(c(2,3,4,5,6), times=c(length(overlap2), length(overlap3), length(overlap4), length(overlap5), length(overlap6)))
)
ol_w_full$warming = rep("warming", nrow(ol_c_full))


mean2 = aggregate(overlap2, by=list(dist2), FUN=mean)
mean3 = aggregate(overlap3, by=list(dist3), FUN=mean)
mean4 = aggregate(overlap4, by=list(dist4), FUN=mean)
mean5 = aggregate(overlap5, by=list(dist5), FUN=mean)
mean6 = aggregate(overlap6, by=list(dist6), FUN=mean)

ol_w = rbind(mean2, mean3, mean4, mean5, mean6)
ol_w$zeta = rep(c(2,3,4,5,6), times=c(nrow(mean2), nrow(mean3), nrow(mean4), nrow(mean5), nrow(mean6)))
ol_w$warming = rep("warming", nrow(ol_w))

ol = rbind(ol_c, ol_w)
ol_full = rbind(ol_c_full, ol_w_full)

write.csv(ol_full, "node_overlap.csv")

# plot Extended Data Fig.5b - node overlap
ggplot(ol_full, aes(x = factor(zeta), y=overlap, fill=warming)) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  geom_jitter(shape=16, size=0.5,  position=position_jitterdodge(jitter.width = 0.01, dodge.width = 0.8)) +
  scale_fill_manual(values=c("#214da0", "#e7211f")) +
  labs(x="Order", y = "Number of overlapping nodes") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        legend.position="none")

# plot Extended Data Fig.5c - node overlap
get_text <- function(y, x){
  lm = lm(y~x)
  lm_p = round(summary(lm)$coefficients[2,4],3)
  lm_r = round(summary(lm)$adj.r.squared,3)
  lm_s = round(summary(lm)$coefficients[2,1],3)
  txt = paste(lm_s, " (", lm_r, ", ", lm_p, ")", sep="")
  return(txt)
}

plot_op_full <- function(xc, yc, xw, yw){
  yl=c(min(yc, yw), max(yc, yw))

  plot(yw~xw, ylim=yl, ylab="Number of overlapping nodes", xlab="Average time gap (year)", pch=19, cex.axis=0.7, cex.lab=0.7, tck=-0.03, col="#e7211f")
  abline(lm(yw~xw), col="#e7211f")

  points(yc~xc, col="#214da0")
  abline(lm(yc~xc), col="#214da0")

  txt1 = get_text(yw,xw)
  txt2 = get_text(yc,xc)
  mtext(txt1, side=3, line=0.8, cex=0.5, col="#e7211f")
  mtext(txt2, side=3, line=0.1, cex=0.5, col="#214da0")
}

plot_op_full(xc=ol_full$gap[which(ol_full$warming=="control" & ol_full$zeta==2)],
             yc=ol_full$overlap[which(ol_full$warming=="control" & ol_full$zeta==2)],
             xw=ol_full$gap[which(ol_full$warming=="warming" & ol_full$zeta==2)],
             yw=ol_full$overlap[which(ol_full$warming=="warming" & ol_full$zeta==2)])

plot_op_full(xc=ol_full$gap[which(ol_full$warming=="control" & ol_full$zeta==3)],
             yc=ol_full$overlap[which(ol_full$warming=="control" & ol_full$zeta==3)],
             xw=ol_full$gap[which(ol_full$warming=="warming" & ol_full$zeta==3)],
             yw=ol_full$overlap[which(ol_full$warming=="warming" & ol_full$zeta==3)])

plot_op_full(xc=ol_full$gap[which(ol_full$warming=="control" & ol_full$zeta==4)],
             yc=ol_full$overlap[which(ol_full$warming=="control" & ol_full$zeta==4)],
             xw=ol_full$gap[which(ol_full$warming=="warming" & ol_full$zeta==4)],
             yw=ol_full$overlap[which(ol_full$warming=="warming" & ol_full$zeta==4)])

plot_op_full(xc=ol_full$gap[which(ol_full$warming=="control" & ol_full$zeta==5)],
             yc=ol_full$overlap[which(ol_full$warming=="control" & ol_full$zeta==5)],
             xw=ol_full$gap[which(ol_full$warming=="warming" & ol_full$zeta==5)],
             yw=ol_full$overlap[which(ol_full$warming=="warming" & ol_full$zeta==5)])

plot_op_full(xc=ol_full$gap[which(ol_full$warming=="control")],
             yc=ol_full$overlap[which(ol_full$warming=="control")],
             xw=ol_full$gap[which(ol_full$warming=="warming")],
             yw=ol_full$overlap[which(ol_full$warming=="warming")])

wilcox.test(ol_full$overlap[which(ol_full$zeta==2 & ol_full$warming=="control")], ol_full$overlap[which(ol_full$zeta==2 & ol_full$warming=="warming")])
wilcox.test(ol_full$overlap[which(ol_full$zeta==3 & ol_full$warming=="control")], ol_full$overlap[which(ol_full$zeta==3 & ol_full$warming=="warming")])
wilcox.test(ol_full$overlap[which(ol_full$zeta==4 & ol_full$warming=="control")], ol_full$overlap[which(ol_full$zeta==4 & ol_full$warming=="warming")])
wilcox.test(ol_full$overlap[which(ol_full$zeta==5 & ol_full$warming=="control")], ol_full$overlap[which(ol_full$zeta==5 & ol_full$warming=="warming")])
wilcox.test(ol_full$overlap[which(ol_full$warming=="control")], ol_full$overlap[which(ol_full$warming=="warming")])
