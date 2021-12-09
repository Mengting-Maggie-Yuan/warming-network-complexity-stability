# contributor: Maggie Yuan

find_overlap = function(mods, bigtable){
	vec1 = bigtable[,which(names(bigtable)==mods[1])]
	vec2 = bigtable[,which(names(bigtable)==mods[2])]
	return(sum(vec1*vec2==1))
}

find_only_in_1 = function(mods, bigtable){
	vec1 = bigtable[,which(names(bigtable)==mods[1])]
	vec2 = bigtable[,which(names(bigtable)==mods[2])]
	return(sum(vec1==1 & vec2==0))
}

find_only_in_2 = function(mods, bigtable){
	vec1 = bigtable[,which(names(bigtable)==mods[1])]
	vec2 = bigtable[,which(names(bigtable)==mods[2])]
	return(sum(vec2==1 & vec1==0))
}

find_N = function(mods, mapping, bigtable){
	nwk1 = bigtable[, which(mapping$year_warming== mapping $year_warming[which(mapping $module==mods[1])])]
	nwk2 = bigtable[, which(mapping$year_warming== mapping $year_warming[which(mapping $module==mods[2])])]
	match_nwk1_nwk2 = sum((rowSums(nwk1) + rowSums(nwk2))>0)
	return(match_nwk1_nwk2)
}

fisher_test = function(x){
	contingency_table <- matrix(unlist(matrix(data.frame(x[3:6]), nrow=2)), nrow=2)
	test_p = fisher.test(contingency_table, alternative = "greater")$p.value
	return(test_p)
}

setwd('/Users/maggieyuan/Documents/!annual_network/GitHub/ExtendedDataFig3.module_preservation/')

node = read.table("input_file/NetworkNodeTable_AllNetworks.txt", header=T, sep="\t")
node$module = paste(node$year_trt, "_M", node$No..module, sep="")
node_table = node[,c(1,17,18)]

module_list = unique(node_table$module)

mytable = data.frame(matrix(NA, nrow=0, ncol=length(module_list)))
colnames(mytable) = module_list
to_add = data.frame(matrix(NA, nrow=1, ncol=length(module_list)))
colnames(to_add) = module_list

map = data.frame(module=module_list, year=gsub("_.*", "", module_list), year_warming=gsub("_M.*", "", module_list))
map$warming=gsub(".*_", "", map$year_warming)

for (i in 1:nrow(node_table))
{
	if (node_table$Name[i] %in% row.names(mytable)){
		mytable[which(row.names(mytable)== node_table$Name[i]), which(names(mytable)== node_table$module[i])] <-1
	}else{
		mytable = rbind(mytable, to_add)
		rownames(mytable)[nrow(mytable)] <-as.character(node_table$Name[i])
		mytable[nrow(mytable), which(names(mytable)== node_table$module[i])] <-1
	}
}
mytable[is.na(mytable)]<-0
nrow(node_table)
sum(mytable)
dim(mytable)

mytable_kp = mytable[,-which(colSums(mytable)<5)]
map_kp = map[-which(colSums(mytable)<5),]
dim(mytable_kp)

###

network_pair_c = t(combn(unique(grep("control", map_kp $year_warming, value=T)), m=2))
network_pair_w = t(combn(unique(grep("warm", map_kp $year_warming, value=T)), m=2))

 network_pair_year = matrix(c("2010_control", "2011_control", "2012_control", "2013_control", "2014_control", "2010_warming", "2011_warming", "2012_warming", "2013_warming", "2014_warming"), ncol=2)

network_pair = rbind(network_pair_c, network_pair_w, network_pair_year)

total_mod_pairs = matrix(NA, nrow=nrow(network_pair), ncol=3)
for (i in 1:nrow(network_pair))
{
	module_pair = as.matrix(expand.grid(map_kp $module[which(map_kp $year_warming==network_pair[i,1])], map_kp $module[which(map_kp $year_warming==network_pair[i,2])]))
	total_mod_pairs[i,] = c(network_pair[i,], nrow(module_pair))
}
total_mod_pairs

sig_mod_pairs = matrix(NA, nrow=0, ncol=4)
sig_detailed_table = c("module1", "module2", "both", "P1A2", "P2A1", "A1A2", "p_raw", "p_adj")
for (i in 1:nrow(network_pair))
{
	module_pair = as.matrix(expand.grid(map_kp $module[which(map_kp $year_warming==network_pair[i,1])], map_kp $module[which(map_kp $year_warming==network_pair[i,2])]))
	overlap = apply(module_pair, 1, FUN= find_overlap, bigtable= mytable_kp)
	only1 = apply(module_pair, 1, FUN= find_only_in_1, bigtable= mytable_kp)
	only2 = apply(module_pair, 1, FUN= find_only_in_2, bigtable= mytable_kp)
	denominator = apply(module_pair, 1, FUN= find_N, mapping=map, bigtable= mytable)
	none = denominator-(overlap + only1 + only2)
	
	count_table = data.frame(module1 = module_pair[,1], module2 = module_pair[,2], Both=overlap, P1A2=only1, P2A1=only2, A1A2=none)
	
	p_raw=c()
	for (tt in 1:nrow(count_table))
	{
		x=count_table[tt,]
		p = fisher_test(x)
		p_raw = c(p_raw, p)
	}
	
	count_table$p_raw = p_raw
	count_table$p_adj = p.adjust(count_table$p_raw, method = "bonferroni")

	network1 = network_pair[i,1]
	network2 = network_pair[i,2]
	sig_count = sum(count_table$p_adj<=0.05)
	
	if(sig_count>0){
		sig_pairs_table = count_table[which(count_table$p_adj<=0.05),c(1:2)]
		sig_pairs_linked = paste(sig_pairs_table[,1], "-", sig_pairs_table[,2], sep="")
		sig_pairs = paste(sig_pairs_linked, collapse=",")
		
		sig_pairs_count_table = count_table[which(count_table$p_adj<=0.05),]
		row.names(sig_pairs_count_table) = sig_pairs_linked
	}else{
		sig_pairs = "None"
	}
	
	add_one_row = c(network1, network2, sig_count, sig_pairs)
	sig_mod_pairs = rbind(sig_mod_pairs, add_one_row)
	
	sig_detailed_table = rbind(sig_detailed_table, sig_pairs_count_table)
	
	print(i)
}
sig_mod_pairs
dim(sig_detailed_table)
head(sig_detailed_table)
write.table(sig_detailed_table[-1,], "preserved_module_pairs.txt", sep="\t")
# module1: module in the first network
# module2: module in the second network
# Both: node number present in both module1 and module2
# P1A2: node number present in module1 but absent from module2
# P2A1: node number present in module2 but absent from module1
# A1A2: node number absent from module1 and module2, but present in either networks
# p_raw and p_adj: raw and Bonferroni adjusted p values for Fisher's exact test


