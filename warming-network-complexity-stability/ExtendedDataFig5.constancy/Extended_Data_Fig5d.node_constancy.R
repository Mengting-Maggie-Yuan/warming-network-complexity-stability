# contributor: Maggie Yuan

setwd("/Users/maggieyuan/Documents/!annual_network/GitHub/ExtendedDataFig5.constancy/")

edge_prop = read.table("input_file/NetworkEdgeTable_AllNetworks.txt", sep="\t", header=T)
edge_prop$edge_sign = paste(edge_prop$edge,edge_prop$V5,sep="_")

total_e = length(unique(edge_prop$edge))
total_e_sign = length(unique(edge_prop$edge_sign)) ##note: signs are consistant

edge_table = data.frame(matrix(NA, nrow=total_e, ncol=11))
row.names(edge_table) = unique(edge_prop$edge)
names(edge_table) = unique(edge_prop$year_trt)

for (i in 1:ncol(edge_table))
{
	edge_prop_year = edge_prop[which(edge_prop$year_trt == names(edge_table[i])),]
	edge_table[which(rownames(edge_table) %in% edge_prop_year$edge),i] = 1
}
edge_table[is.na(edge_table)]<-0

head(edge_table)

edge_table$avg_c = rowMeans(edge_table[,c(1,2,4,6,8,10)])
edge_table$sd_c = apply(edge_table[,c(1,2,4,6,8,10)],1,FUN=sd)
edge_table$avg_w = rowMeans(edge_table[,c(1,3,5,7,9,11)])
edge_table$sd_w = apply(edge_table[,c(1,3,5,7,9,11)],1,FUN=sd)

edge_table$constancy_c = edge_table$avg_c/edge_table$sd_c
edge_table$constancy_w = edge_table$avg_w/edge_table$sd_w

edge_table$constancy_c[is.infinite(edge_table$constancy_c)] = NA
edge_table$constancy_w[is.infinite(edge_table$constancy_w)] = NA

# edge_table[1:100,]

v_constancy_c = edge_table$constancy_c[which(!is.na(edge_table$constancy_c))]
n_constancy_c = sum(!is.na(edge_table$constancy_c))
avg_constancy_c = mean(edge_table$constancy_c, na.rm=T)
sd_constancy_c = sd(edge_table$constancy_c, na.rm=T)

v_constancy_w = edge_table$constancy_w[which(!is.na(edge_table$constancy_w))]
n_constancy_w = sum(!is.na(edge_table$constancy_w))
avg_constancy_w = mean(edge_table$constancy_w, na.rm=T)
sd_constancy_w = sd(edge_table$constancy_w, na.rm=T)

t.test(v_constancy_w, v_constancy_c)

library(ggplot2)

# plot Extended Data Fig. 5d - link constancy
df_lc_uw = data.frame(Treatment = c("control", "warming"), Link_constancy = c(avg_constancy_c, avg_constancy_w), se = c(sd_constancy_c/sqrt(n_constancy_c), sd_constancy_w/sqrt(n_constancy_w)))

ggplot(df_lc_uw, aes(x=Treatment, y= Link_constancy)) +
  geom_bar(stat="identity", fill= c("#214da0", "#e7211f"), width = 0.5) +
  geom_errorbar(aes(ymin = Link_constancy - se, ymax = Link_constancy + se), width=0.2) +
  labs(x="Treatment", y = "Unweighted link constancy") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        legend.position="none")
