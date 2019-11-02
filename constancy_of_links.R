rm(list=ls())

setwd("/Users/maggieyuan/Documents/annual_network/datatables/")

edge_prop = read.table("node_edge_files/16S_all_edge_table_09to14.txt", sep="\t", header=T)
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

trt_group_for_plot = c(rep("AmbientT", length(v_constancy_c)), rep("Warming", length(v_constancy_w)))
boxplot(c(v_constancy_c, v_constancy_w) ~ trt_group_for_plot)

hist(v_constancy_c)
hist(v_constancy_w)

for_plot = data.frame(constancy = c(v_constancy_c, v_constancy_w), trt = trt_group_for_plot)

library(ggplot2)
ggplot(for_plot, aes(x= constancy), colour=factor(trt)) + geom_histogram()+ facet_grid(~trt)+ theme_bw() + scale_y_continuous(name="Count (log-scale)", trans = "log10")

