
module_names = mod_names(cemres)
hubs_top20 = get_hubs(cem=cemres, n=20) # Returns n genes in each module with high connectivity, Default: "adjacency"
hubs = get_hubs(cem=cemres, n="all", method = "kME")
moduleinfo = matrix(nrow=length(module_names), ncol=39, byrow=T)
for(i in 1:length(module_names)) {
	mgenes = module_genes(cemres, module=module_names[i])
	mgenesvec = as_tibble(mgenes$genes)
	mgenesvec

	eigen = eig %>% dplyr::select(c("subjectID",module_names[i]))
	eigen = eigen %>% inner_join(phenodata, by = c("subjectID" = "SampleID"))
	eigen
	res = wilcox.test(eval(parse(text=module_names[i])) ~ dx, data = eigen, exact = FALSE)
	u = res$p.value

	kME = hubs %>% pluck(module_names[i]) %>% as_tibble(rownames="Genes")

	a = lncRNA %>% semi_join(mgenesvec, by = c("idents" = "value"))
	a

	e = pc %>% semi_join(mgenesvec, by = c("idents" = "value"))
	e

	####################################################################
	#
	# pvalue < 0.05
	#
	####################################################################

	b = lncRNA_p0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	b
	b_out = dim(lncRNA_p0.05)[1] - dim(b)[1]

	k = lncRNA_pgreater0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	k
	k_out = dim(lncRNA_pgreater0.05)[1] - dim(k)[1]

	c = lncRNA_p0.05_up %>% semi_join(mgenesvec, by = c("idents" = "value"))
	c

	d = lncRNA_p0.05_down %>% semi_join(mgenesvec, by = c("idents" = "value"))
	d

	f = pc_p0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	f
	f_out = dim(pc_p0.05)[1] - dim(f)[1]

	w = pc_pgreater0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	w
	w_out = dim(pc_pgreater0.05)[1] - dim(w)[1]

	g = pc_p0.05_up %>% semi_join(mgenesvec, by = c("idents" = "value"))
	g

	h = pc_p0.05_down %>% semi_join(mgenesvec, by = c("idents" = "value"))
	h

	####################################################################
	#
	# qvalue < 0.05
	#
	####################################################################

	l = lncRNA_q0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	l
	l_out = dim(lncRNA_q0.05)[1] - dim(l)[1]

	m = lncRNA_qgreater0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	m
	m_out = dim(lncRNA_pgreater0.05)[1] - dim(m)[1]

	n = lncRNA_q0.05_up %>% semi_join(mgenesvec, by = c("idents" = "value"))
	n

	o = lncRNA_q0.05_down %>% semi_join(mgenesvec, by = c("idents" = "value"))
	o

	q = pc_q0.05 %>% semi_join(mgenesvec, by = c("idents" = "value"))
	q

	r = pc_q0.05_up %>% semi_join(mgenesvec, by = c("idents" = "value"))
	r

	s = pc_q0.05_down %>% semi_join(mgenesvec, by = c("idents" = "value"))
	s

	t = fisher.test(matrix(c(dim(b)[1], b_out, dim(k)[1], k_out),2,2), alternative='greater')$p.value

	v = fisher.test(matrix(c(dim(f)[1], f_out, dim(w)[1], w_out),2,2), alternative='greater')$p.value

	x = fisher.test(matrix(c(dim(f)[1]+dim(b)[1], f_out+b_out, dim(w)[1]+dim(k)[1], w_out+k_out),2,2), alternative='greater')$p.value

	####################################################################
	#
	# Mean module memberships (KME) of gene sets/list of genes
	#
	####################################################################

	kME_a = a %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_a = mean(kME_a$value)
	mu_kME_a

	kME_b = b %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_b = mean(kME_b$value)
	mu_kME_b

	kME_c = c %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_c = mean(kME_c$value)
	mu_kME_c

	kME_d = d %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_d = mean(kME_d$value)
	mu_kME_d

	kME_e = e %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_e = mean(kME_e$value)
	mu_kME_e

	kME_f = f %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_f = mean(kME_f$value)
	mu_kME_f

	kME_g = g %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_g = mean(kME_g$value)
	mu_kME_g

	kME_h = h %>% inner_join(kME, by = c("idents" = "Genes"))
	mu_kME_h = mean(kME_h$value)
	mu_kME_h

	moduleinfo[i,1] = bregion # brain region name
	moduleinfo[i,2] = module_names[i] # module names
	moduleinfo[i,3] = dim(mgenesvec)[1] # number of genes in the module
	moduleinfo[i,4] = dim(a)[1] # lncRNA in the module
	moduleinfo[i,5] = dim(b)[1] # lncRNA_p0.05 in the module
	moduleinfo[i,6] = dim(c)[1] # lncRNA_p0.05_up in the module
	moduleinfo[i,7] = dim(d)[1] # lncRNA_p0.05_down in the module
	moduleinfo[i,8] = b_out # lncRNA_p0.05 outside the module
	moduleinfo[i,9] = dim(k)[1] # lncRNA_pgreater0.05 in the module
	moduleinfo[i,10] = k_out # lncRNA_pgreater0.05 outside the module
	moduleinfo[i,11] = t # fisher test for overrepresentation of lncRNA_p0.05 in the module
	moduleinfo[i,12] = v # fisher test for overrepresentation of pc_p0.05 in the module
	moduleinfo[i,13] = x # fisher test for overrepresentation of pc_p0.05+lncRNA_p0.05 in the module
	moduleinfo[i,14] = u # wilcoxon rank sum test for testing the difference in median of eigen gene expression between cases and controls for each module 
	moduleinfo[i,15] = dim(e)[1] # pc in the module
	moduleinfo[i,16] = dim(f)[1] # pc_p0.05 in the module
	moduleinfo[i,17] = dim(g)[1] # pc_p0.05_up in the module
	moduleinfo[i,18] = dim(h)[1] # pc_p0.05_down in the module
	moduleinfo[i,19] = f_out # pc_p0.05 outside the module
	moduleinfo[i,20] = dim(w)[1] # pc_pgreater0.05 in the module
	moduleinfo[i,21] = w_out # pc_pgreater0.05 outside the module
	moduleinfo[i,22] = dim(l)[1] # lncRNA_q0.05 in the module
	moduleinfo[i,23] = dim(n)[1] # lncRNA_q0.05_up in the module
	moduleinfo[i,24] = dim(o)[1] # lncRNA_q0.05_down in the module
	moduleinfo[i,25] = l_out # lncRNA_q0.05 outside the module
	moduleinfo[i,26] = dim(m)[1] # lncRNA_qgreater0.05 in the module
	moduleinfo[i,27] = m_out # lncRNA_qgreater0.05 outside the module
	moduleinfo[i,28] = dim(q)[1] # pc_q0.05 in the module
	moduleinfo[i,29] = dim(r)[1] # pc_q0.05_up in the module
	moduleinfo[i,30] = dim(s)[1] # pc_q0.05_down in the module
	moduleinfo[i,31] = mu_kME_a # mean module membership kME lncRNA in the module
	moduleinfo[i,32] = mu_kME_b # mean module membership kME lncRNA_p0.05 in the module
	moduleinfo[i,33] = mu_kME_c # mean module membership kME lncRNA_p0.05_up in the module
	moduleinfo[i,34] = mu_kME_d # mean module membership kME lncRNA_p0.05_down in the module
	moduleinfo[i,35] = mu_kME_e # mean module membership kME pc in the module
	moduleinfo[i,36] = mu_kME_f # mean module membership kME pc_p0.05 in the module
	moduleinfo[i,37] = mu_kME_g # mean module membership kME pc_p0.05_up in the module
	moduleinfo[i,38] = mu_kME_h # mean module membership kME pc_p0.05_down in the module
	moduleinfo[i,39] = paste(names(hubs_top20[[module_names[i]]]), collapse="/") #top 20 most connected genes in each module by adjacency method
}

moduleinfo = as.data.frame(moduleinfo)
colnames(moduleinfo) = c("Brain_Region",
		"ModuleName",
		"number_genes_in_module",
		"lncRNA",
		"lncRNA_p0.05",
		"lncRNA_p0.05_up",
		"lncRNA_p0.05_down",
		"lncRNA_p0.05_out",
		"lncRNA_pgreater0.05",
		"lncRNA_pgreater0.05_out",
		"lncRNA_p0.05_fisher_test_pvalue",
		"pc_p0.05_fisher_test_pvalue",
		"lncRNA_pc_p0.05_fisher_test_pvalue",
		"wilcoxon_test_pvalue",
		"pc",
		"pc_p0.05",
		"pc_p0.05_up",
		"pc_p0.05_down",
		"pc_p0.05_out",
		"pc_pgreater0.05",
		"pc_pgreater0.05_out",
		"lncRNA_q0.05",
		"lncRNA_q0.05_up",
		"lncRNA_q0.05_down",
		"lncRNA_q0.05_out",
		"lncRNA_qgreater0.05",
		"lncRNA_qgreater0.05_out",
		"pc_q0.05",
		"pc_q0.05_up",
		"pc_q0.05_down",
		"mu_kME_lncRNA",
		"mu_kME_lncRNA_p0.05",
		"mu_kME_lncRNA_p0.05_up",
		"mu_kME_lncRNA_p0.05_down",
		"mu_kME_pc",
		"mu_kME_pc_p0.05",
		"mu_kME_pc_p0.05_up",
		"mu_kME_pc_p0.05_down",
		"top20_hub_genes")
write.csv(moduleinfo, file="Modules_Identified.csv", row.names=F)

orares = ora_data(cem=cemres) # Retrieve over representation analysis (ORA) results
moduleenrichmentinfo = data.frame(head(orares[orares$Module %in% module_names[1],], n=3))
for(i in 2:length(module_names)) {
	module_name = module_names[i]
	moduleenrichmentinfo = rbind(moduleenrichmentinfo, head(orares[orares$Module %in% module_name,], n=3))
}
colnames(moduleenrichmentinfo) = c("Module","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count")
moduleenrichmentinfo = as.data.frame(moduleenrichmentinfo, colClasses=c("character","character","character","numeric","numeric","numeric","numeric","numeric","character","numeric"), stringsAsFactors=FALSE, make.names=T)
moduleenrichmentinfo_moduleinfo_merged = merge(x=moduleenrichmentinfo, y=moduleinfo, by.x="Module", by.y="ModuleName", all.x=T)
write.csv(moduleenrichmentinfo, file="Modules_Identified_Enrichment_top3_pathways.csv", row.names=F)
write.csv(moduleenrichmentinfo_moduleinfo_merged, file="Modules_Identified_Enrichment_top3_pathways_with_module_info.csv", row.names=F)
write.csv(orares, file="Modules_Identified_Enrichment_all_pathways.csv")

orares_moduleinfo_merged = merge(x=orares, y=moduleinfo, by.x="Module", by.y="ModuleName", all.x=T)
write.csv(orares_moduleinfo_merged, file="Modules_Identified_Enrichment_all_pathways_with_module_info.csv")
