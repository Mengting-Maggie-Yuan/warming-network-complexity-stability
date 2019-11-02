rm(list=ls())
library(vegan)

setwd("/Users/maggieyuan/Documents/annual_network/env_filtering/find_link")

test_link = function(link_table_row, OTUabd, geodist){
	OTU1 = as.character(as.matrix(link_table_row[1]))
	OTU2	 = as.character(as.matrix(link_table_row[3]))
	OTU1abd = as.numeric(OTUabd[which(row.names(OTUabd)== OTU1),])
	OTU2abd = as.numeric(OTUabd[which(row.names(OTUabd)== OTU2),])
	distOTU1 = vegdist(cbind(OTU1abd,rep(0,length(OTU1abd))), method="bray")
	distOTU2 = vegdist(cbind(OTU2abd,rep(0,length(OTU2abd))), method="bray")
	
	cor1 = cor.test(distOTU1, geodist)
	cor2 = cor.test(distOTU2, geodist)
	
	return(c(OTU1, OTU2, cor1$estimate, cor1$p.value, cor2$estimate, cor2$p.value))
}

dispersal_limitation_link = function(edge_table, OTU_table, distance, p, cutoff){
	edges = read.table(edge_table, sep=" ")
	OTUs = read.table(OTU_table, sep="\t", header=T)
	distances = read.table(distance, head=T, row.names=1)
	
	OTUs_order = OTUs[,order(names(OTUs))]
	distances_order = distances[order(row.names(distances)),]

	if(sum(names(OTUs_order)!=row.names(distances_order))!=0){print("check sample names!")}else{
		geodist = vegdist(distances_order, method="euclid")
	}
	
	test_result = t(apply(edges, 1, FUN=test_link, OTUabd=OTUs_order, geodist=geodist)) # OK - warnings: In vegdist(cbind(OTU2abd, rep(0, length(OTU2abd))), method = "bray") : you have empty rows: their dissimilarities may be meaningless in method “bray”.
	test_result_pick = test_result[which((test_result[,3]>=cutoff) & (test_result[,5]>=cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]
	
	total_link=nrow(edges)
	result= test_result_pick
	
	return(list(total_link, result))
}


dispersal_limitation_link(edge_table="edges_reformat/16S_2009_80.txt", OTU_table="OTUtable/Y09_NW.txt", distance="distance/2009.tsv", p=0.05, cutoff=0)
dispersal_limitation_link(edge_table="edges_reformat/16S_2010_control_80.txt", OTU_table="OTUtable/Y10_N.txt", distance="distance/2010_control.tsv", p=0.05, cutoff=0)
dispersal_limitation_link(edge_table="edges_reformat/16S_2010_warming_80.txt", OTU_table="OTUtable/Y10_W.txt", distance="distance/2010_warming.tsv", p=0.05, cutoff=0)
dispersal_limitation_link(edge_table="edges_reformat/16S_2011_control_80.txt", OTU_table="OTUtable/Y11_N.txt", distance="distance/2011_control.tsv", p=0.05, cutoff=0)
dispersal_limitation_link(edge_table="edges_reformat/16S_2011_warming_80.txt", OTU_table="OTUtable/Y11_W.txt", distance="distance/2011_warming.tsv", p=0.05, cutoff=0)
dispersal_limitation_link(edge_table="edges_reformat/16S_2012_control_81.txt", OTU_table="OTUtable/Y12_N.txt", distance="distance/2012_control.tsv", p=0.05, cutoff=0)
dispersal_limitation_link(edge_table="edges_reformat/16S_2012_warming_80.txt", OTU_table="OTUtable/Y12_W.txt", distance="distance/2012_warming.tsv", p=0.05, cutoff=0)
dispersal_limitation_link(edge_table="edges_reformat/16S_2013_control_80.txt", OTU_table="OTUtable/Y13_N.txt", distance="distance/2013_control.tsv", p=0.05, cutoff=0)
dispersal_limitation_link(edge_table="edges_reformat/16S_2013_warming_80.txt", OTU_table="OTUtable/Y13_W.txt", distance="distance/2013_warming.tsv", p=0.05, cutoff=0)
dispersal_limitation_link(edge_table="edges_reformat/16S_2014_control_80.txt", OTU_table="OTUtable/Y14_N.txt", distance="distance/2014_control.tsv", p=0.05, cutoff=0)
dispersal_limitation_link(edge_table="edges_reformat/16S_2014_warming_80.txt", OTU_table="OTUtable/Y14_W.txt", distance="distance/2014_warming.tsv", p=0.05, cutoff=0)