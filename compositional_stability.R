rm(list=ls())
library(vegan)

setwd("/Users/maggieyuan/Documents/annual_network/datatables/OTUtable")
list.files()

otu_in_network = read.table("X16S_OTU_in_network.txt", sep="\t", header=T)
map = read.table("trt_09to14.txt", header=T, sep="\t")

trt_annual_list <- expand.grid(list(year=unique(map$Year), warming = unique(map$Warming)))
trt_annual_list_w = trt_annual_list[c(1,8:12),]
trt_annual_list_c = trt_annual_list[c(1,14:18),]

subset_ids = list(c(1,2), c(2,3), c(3,4), c(4,5), c(5,6))

trt = trt_annual_list_c
cs_mean_c = c()
cs_sd_c = c()
cs_se_c = c()

for (i in 1:length(subset_ids))
{
	id1 = which(paste(map$Year,map$Warming)==paste(trt$year,trt$warming)[subset_ids[[i]][1]])	
	id2 = which(paste(map$Year,map$Warming)==paste(trt$year,trt$warming)[subset_ids[[i]][2]])
	
	mutual_plot = intersect(map$Plot_ID[id1], map$Plot_ID[id2])
	mutual_plot_id1 = id1[which(map$Plot_ID[id1] %in% mutual_plot)]
	mutual_plot_id2 = id2[which(map$Plot_ID[id2] %in% mutual_plot)]
	
	mutual_plot_id1_order = mutual_plot_id1[order(map$Plot_ID[mutual_plot_id1])]
	mutual_plot_id2_order = mutual_plot_id2[order(map$Plot_ID[mutual_plot_id2])]
	
	if (all.equal(map$Plot_ID[mutual_plot_id1_order], map$Plot_ID[mutual_plot_id2_order]))
	{
		paired_cs = c()
		for (j in 1:length(mutual_plot_id1_order))
		{
			abd1 = otu_in_network[,mutual_plot_id1_order[j]]
			abd2 = otu_in_network[,mutual_plot_id2_order[j]]
			toadd = sqrt(1-vegdist(rbind(abd1, abd2)))			
			paired_cs = c(paired_cs, toadd)
		}
		add_mean = mean(paired_cs)
		add_sd = sd(paired_cs)
		add_se = add_sd/sqrt(length(paired_cs))
		
		cs_mean_c = c(cs_mean_c, add_mean)
		cs_sd_c = c(cs_sd_c, add_sd)
		cs_se_c = c(cs_se_c, add_se)
			
	}else{
		print("ids do not match")
	}
}


trt = trt_annual_list_w
cs_mean_w = c()
cs_sd_w = c()
cs_se_w = c()

for (i in 1:length(subset_ids))
{
	id1 = which(paste(map$Year,map$Warming)==paste(trt$year,trt$warming)[subset_ids[[i]][1]])	
	id2 = which(paste(map$Year,map$Warming)==paste(trt$year,trt$warming)[subset_ids[[i]][2]])
	
	mutual_plot = intersect(map$Plot_ID[id1], map$Plot_ID[id2])
	mutual_plot_id1 = id1[which(map$Plot_ID[id1] %in% mutual_plot)]
	mutual_plot_id2 = id2[which(map$Plot_ID[id2] %in% mutual_plot)]
	
	mutual_plot_id1_order = mutual_plot_id1[order(map$Plot_ID[mutual_plot_id1])]
	mutual_plot_id2_order = mutual_plot_id2[order(map$Plot_ID[mutual_plot_id2])]
	
	if (all.equal(map$Plot_ID[mutual_plot_id1_order], map$Plot_ID[mutual_plot_id2_order]))
	{
		paired_cs = c()
		for (j in 1:length(mutual_plot_id1_order))
		{
			abd1 = otu_in_network[,mutual_plot_id1_order[j]]
			abd2 = otu_in_network[,mutual_plot_id2_order[j]]
			toadd = sqrt(1-vegdist(rbind(abd1, abd2)))			
			paired_cs = c(paired_cs, toadd)
		}
		add_mean = mean(paired_cs)
		add_sd = sd(paired_cs)
		add_se = add_sd/sqrt(length(paired_cs))
		
		cs_mean_w = c(cs_mean_w, add_mean)
		cs_sd_w = c(cs_sd_w, add_sd)
		cs_se_w = c(cs_se_w, add_se)
			
	}else{
		print("ids do not match")
	}
}

result = data.frame(cs_mean_c, cs_sd_c, cs_se_c, cs_mean_w, cs_sd_w, cs_se_w)
row.names(result) = c("Y0-Y1", "Y1-Y2", "Y2-Y3", "Y3-Y4", "Y4-Y5")

write.table(result, "/Users/maggieyuan/Documents/annual_network/stability/compositional_stability.txt", sep="\t")


xaxis = c(1,2,3,4,5)

plot(result$cs_mean_w~ xaxis, col="red", main="16S-ITS")
abline(lm(result$cs_mean_w~ xaxis), col="red")
points(result$cs_mean_c~ xaxis, col="black")
abline(lm(result$cs_mean_c~ xaxis), col="black")

summary(lm(result$cs_mean_w~ xaxis))
summary(lm(result$cs_mean_c~ xaxis))



