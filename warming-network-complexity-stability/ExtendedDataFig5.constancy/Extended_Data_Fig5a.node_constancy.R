# contributor: Maggie Yuan

library(ggplot2)
library(gridExtra)

setwd("/Users/maggieyuan/Documents/!annual_network/GitHub/ExtendedDataFig5.constancy/")

otu = read.table("input_file/OTUtable_NetworkedOTUs_AllSamples.txt", sep="\t", header=T, row.names=1)
map = read.table("input_file/SampleMap_AllSamples.txt", sep="\t", header=T)

id_09 = rep(which(map$Year == "Y09"), each=2)
id_10 = which(map$Year == "Y10")
id_11 = which(map$Year == "Y11")
id_12 = which(map$Year == "Y12")
id_13 = which(map$Year == "Y13")
id_14 = which(map$Year == "Y14")

# check plot order
data.frame(map$Plot_full_name[id_09],
           map$Plot_full_name[id_10],
           map$Plot_full_name[id_11],
           map$Plot_full_name[id_12],
           map$Plot_full_name[id_13],
           map$Plot_full_name[id_14])
warming_trt = map$Warming[id_14]

# separate OTU table
otu_09 = as.data.frame(otu[, id_09])
otu_10 = as.data.frame(otu[, id_10])
otu_11 = as.data.frame(otu[, id_11])
otu_12 = as.data.frame(otu[, id_12])
otu_13 = as.data.frame(otu[, id_13])
otu_14 = as.data.frame(otu[, id_14])

# check OTU table order
sum(names(otu_09) != map$Sample[id_09]) # 24 doesn't match because auto changed column names
sum(names(otu_10) != map$Sample[id_10])
sum(names(otu_11) != map$Sample[id_11])
sum(names(otu_12) != map$Sample[id_12])
sum(names(otu_13) != map$Sample[id_13])
sum(names(otu_14) != map$Sample[id_14])

# calculate constancy
otu_mean = (otu_09+otu_10+otu_11+otu_12+otu_13+otu_14)/6
otu_sd = sqrt(((otu_09-otu_mean)^2 + (otu_10-otu_mean)^2 + (otu_11-otu_mean)^2 + (otu_12-otu_mean)^2 + (otu_13-otu_mean)^2 + (otu_14-otu_mean)^2)/5)
otu_constancy = otu_mean/otu_sd

write.table(otu_constancy, "observed_constancy_of_each_node.csv", sep=",")

otu_constancy_w = otu_constancy[,which(warming_trt == "W")]
otu_constancy_c = otu_constancy[,which(warming_trt == "N")]

otu_constancy_w_avg = rowMeans(otu_constancy_w)
otu_constancy_c_avg = rowMeans(otu_constancy_c)

con_w = otu_constancy_w_avg[is.finite(otu_constancy_w_avg)]
con_c = otu_constancy_c_avg[is.finite(otu_constancy_c_avg)]


# plot Extended Data Fig5a - node constancy
nc_df = data.frame(warming = c(rep("control", length(con_c)), rep("warming", length(con_w))),
                   nc = c(con_c, con_w))

ggplot(nc_df, aes(x=warming, y=nc, fill=warming)) +
  geom_boxplot(alpha=1, width=0.4, outlier.shape = NA) +
  geom_jitter(shape=16, size=0.5, position=position_jitterdodge(jitter.width = 0.01, dodge.width = 0.8)) +
  scale_fill_manual(values=c("#214da0", "#e7211f")) +
  labs(y = "Node constancy") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        legend.position="none")
