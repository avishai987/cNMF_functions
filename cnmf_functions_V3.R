require(reticulate)

sum2one <- function(df) {
  for (col_num in 1:ncol(df)){
    program_sum = sum(df[,col_num])
    norm_program = df[,col_num]/program_sum
    df[,col_num] = norm_program
  }
  return(df)
  
}



#assign program to each cell after add_prgorams_score
#' @param larger_by how much a program should be greater then the second one to be assign as this program of this cell.
program_assignment <- function(dataset,larger_by = 1,program_names) {
  
  scores_metadata = FetchData(object = dataset,vars = c(program_names))
  assignment_df =  data.frame(row.names = rownames(scores_metadata))
  
  find_max_program <- function(x,larger_by) { #find max program index, else: "NA"
    n <- length(x)
    max = sort(x)[n] #get max value
    second_max = sort(x)[n-1] #get 2nd max value
    if (max > larger_by*second_max){program_names[which.max(x)]} else {"NA"} #if max is larger than 2nd max*larger_by, write max index. else, write "NA"
  }
  assignment_df["assigment"] = apply(X = scores_metadata,MARGIN = 1,FUN = find_max_program,larger_by = larger_by) 
  
  dataset = AddMetaData(object = dataset,metadata = assignment_df,col.name = "program.assignment")
  return(dataset)
}


programs_dotplot <- function(seurat_obj,treatment_var,program_names) {
  # return dotplot of z score expression and fraction of assigned cells
  #run program_assignment first
  
  #get program_assignment_data:
    program_assignment_data = FetchData(object = seurat_obj,vars = c(treatment_var,"program.assignment"))%>%
    group_by_at(treatment_var) %>% 
  summarize(
      program.1  = sum(program.assignment == "Program.1"),
      program.2  = sum(program.assignment == "Program.2"),
      program.3  = sum(program.assignment == "Program.3"),
      program.4  = sum(program.assignment == "Program.4"),
      program.5  = sum(program.assignment == "Program.5"),
    ) %>% 
    ungroup() %>% 
    t() %>% as.data.frame() %>% janitor::row_to_names(1) %>% rownames_to_column("program") %>% 
    pivot_longer(!program, names_to = treatment_var, values_to = "assigned_cells") %>% 
    mutate(assigned_cells = as.numeric(assigned_cells))
   
   #get program expression data:
  program_exprs_data = FetchData(object = seurat_obj,vars = c(treatment_var, program_names))%>%
  group_by_at(treatment_var) %>% 
    summarize(across(everything(),mean)) %>% 
    ungroup() %>% 
    t() %>% as.data.frame() %>% janitor::row_to_names(1) %>% rownames_to_column("program") %>% 
    pivot_longer(!program, names_to = treatment_var, values_to = "z_score_level") %>% 
    mutate(program = gsub(program, pattern = "_scaled", replacement = "")) %>% 
    mutate(program = gsub(program, pattern = "Program", replacement = "program"))%>% 
    mutate(z_score_level = as.numeric(z_score_level))
  
  #join and add cells fraction data
  n_cells_treatment = table(seurat_obj[[treatment_var]])
  all_data = full_join(program_assignment_data,program_exprs_data,by = c("program" ,  treatment_var))
  all_data = all_data %>% mutate(cells_in_treatment = as.vector(n_cells_treatment[all_data[[treatment_var]]])) %>% mutate(cells_fraction = assigned_cells/cells_in_treatment)
  
  # plot
  p = ggplot(data = all_data, mapping = aes_string(x = "program", 
      y = treatment_var)) + geom_point(mapping = aes_string(size = "cells_fraction", 
      color = 'z_score_level')) + scale_size(range = c(1, 8)) + theme(axis.title.x = element_blank(), 
      axis.title.y = element_blank()) + guides(size = guide_legend(title = "assigned cells fraction")) + 
      labs(x = "Features", y = treatment_var) + cowplot::theme_cowplot()+
    scale_color_gradient(low = "lightblue", high = "darkblue")
  # change x axis names to program names
   p+scale_x_discrete(labels= paste0(program_names,"\n(",c("IFNa","TNFa-NFKb","HIF","Cell_Cycle"),")"))
  return(p)
}
cell_percentage = function(dataset,time.point_var, by_program = F, by_tp = F,x_order = NULL) {
  if(by_program){
    data =FetchData(object = dataset,vars = c("program.assignment",time.point_var))
    data = data %>% dplyr::count(program.assignment, .[time.point_var]) %>%  dplyr::add_count(program.assignment, wt = n, name = "overall")%>% 
      mutate(proportion = n / overall)   
    
    plt_list = list()
    time.point_var = ensym(time.point_var)
    for (program_name in unique(data$program.assignment)) {
      program_data = data[data$program.assignment == program_name,]
      p = ggplot(data=program_data, aes(x=!!time.point_var, y=proportion)) +geom_bar(stat="identity")+ylab("precentage") +
        ggtitle("program" %>% paste(program_data$program.assignment %>% unique() %>% as.character()))+
        scale_y_continuous(limits = c(0,1))  + scale_x_discrete(limits = x_order)
      plt_list[[program_name]] = p
    }
    p = ggarrange(plotlist = plt_list )
    return(p)
  }
  else if( by_tp){
    data =FetchData(object = dataset,vars = c("program.assignment",time.point_var))
    data = data %>% dplyr::count(program.assignment, .[time.point_var]) %>%  dplyr::add_count(.[time.point_var], wt = n, name = "overall")%>% 
      mutate(proportion = n / overall)   
    
    plt_list = list()
    for (tp in unique(data[,time.point_var])) {
      program_data = data[data[,time.point_var] == tp,]
      p = ggplot(data=program_data, aes(x=program.assignment, y=proportion)) +geom_bar(stat="identity")+ylab("precentage") +
        ggtitle(program_data[,time.point_var] %>% unique() %>% as.character())+
        scale_y_continuous(limits = c(0,1))   + scale_x_discrete(limits = x_order)
      plt_list[[tp]] = p
    }
    p = ggarrange(plotlist = plt_list)
    return(p)
  }
  
}
  

#union programs (i.e all cell cycle programs). input example: groups_list = c(5,3,4)
# will combine metagenes 5+3+4 in the df @all_metagenes
union_programs <- function(groups_list,all_metagenes) {
  all_programs = as.list(1:ncol(all_metagenes)) # create list of all programs
  all_programs = all_programs[!all_programs %in% groups_list] # remove groups that wil be combined
  all_programs[["cobined"]] <- groups_list # add combined groups
  
  unioned_metagenes = data.frame(row.names = rownames(all_metagenes)) #create final df 
  for (group in all_programs) {
    name =group %>% as.character() %>% paste(collapse = ".") #name of col, i.e "1.3"
    name = paste0("gep",name) #change to gep1.3
    col = all_metagenes[,group,drop=F] #take cols
    if (ncol(col)!=1){col = rowMeans(col)} # mean if more than 1 metagenes
    new_metagene  = data.frame(new =  col,row.names = rownames(all_metagenes))
    names(new_metagene) = name #set name
    unioned_metagenes = cbind(unioned_metagenes,new_metagene) # add to final df
  }
  return(unioned_metagenes)
}


#' statistically compare mean of metagenes
#'
#' @param dataset 
#' @param time.point_var time point metadata name in dataset
#' @param prefix adding to patient names i.e "model" result "model 1069"
#' @param patient.ident_var patient id metadata name in dataset
#' @param pre_on names of pre and on cells i.e c("OSI","NT")
#'
#' @return none
#' @export
#'
#' @examples
metagenes_mean_compare <- function(dataset,time.point_var,prefix = "",patient.ident_var,pre_on = c("OSI","NT"),axis.text.x = 11,test = "t.test", programs = c("Hypoxia","TNFa","Cell_cycle"), with_split = T, without_split = T){
  
  for (metegene in programs) {
    #create data:
    genes_by_tp = FetchData(object = dataset,vars = metegene) %>% rowSums() %>% as.data.frame() #mean expression
    names(genes_by_tp)[1] = "Metagene_mean"
    genes_by_tp = cbind(genes_by_tp,FetchData(object = dataset,vars = c(patient.ident_var,time.point_var))) # add id and time points
    
    
    genes_by_tp_forPlot =  genes_by_tp %>% mutate(!!ensym(patient.ident_var) := paste(prefix,genes_by_tp[,patient.ident_var])) #add "model" before  each model/patient
    fm <- as.formula(paste("Metagene_mean", "~", time.point_var)) #make formula to plot
    
    #plot and split by patient:   
    stat.test = compare_means(formula = fm ,data = genes_by_tp_forPlot,method = test,group.by = patient.ident_var,p.adjust.method = "fdr")%>% # Add pairwise comparisons p-value
      dplyr::filter(group1 == pre_on[1] & group2 == pre_on[2])  #filter for pre vs on treatment only
    
    plt = ggboxplot(genes_by_tp_forPlot, x = time.point_var, y = "Metagene_mean", color = time.point_var) + #plot
      stat_pvalue_manual(stat.test, label = "p = {p.adj}",  #add p value
                         y.position = max(genes_by_tp_forPlot$Metagene_mean))+ # set position at the top value
      grids()+  
      ylab(paste(metegene,"mean"))+
      theme(axis.text.x = element_text(size = axis.text.x))+
      ylim(0, max(genes_by_tp_forPlot$Metagene_mean)*1.2) # extend y axis to show p value
    
    plt = facet(plt, facet.by = patient.ident_var) #split by patients
    print_tab(plt = plt,title = c(metegene,"per patient")) 
    
    
    #plot = without split by patient:
    if(without_split){
      stat.test = compare_means(formula = fm ,data = genes_by_tp_forPlot,comparisons = my_comparisons,method = test)%>% 
        dplyr::filter(group1 == pre_on[1] & group2 == pre_on[2]) # Add pairwise comparisons p-value
      
      plt = ggboxplot(genes_by_tp_forPlot, x = time.point_var, y = "Metagene_mean", color = time.point_var) +
        stat_pvalue_manual(stat.test, label = "p = {p.adj}",  #add p value
                           y.position = max(genes_by_tp_forPlot$Metagene_mean))+ # set position at the top value
        grids()+  
        ylab(paste(metegene,"mean"))+
        ylim(0, max(genes_by_tp_forPlot$Metagene_mean)*1.2) # extend y axis to show p value
      
      
      print_tab(plt = plt,title = metegene)
    }
    
  }
  
  
}
#
require(facefuns)
metagenes_violin_compare = function(dataset,prefix = "",pre_on = c("OSI","NT"),axis.text.x = 11,test = "t.test", programs = c("Hypoxia","TNFa","Cell_cycle"),
                                    patient.ident_var = "orig.ident", return_list = F){
  plt.lst = list()
  for (metegene in programs) {
    #create data:
    genes_by_tp = FetchData(object = dataset,vars =  c(patient.ident_var,"treatment",metegene)) %>% filter(treatment %in% pre_on)  %>% as.data.frame() #mean expression
    names(genes_by_tp)[3] = "Metagene_mean"
    
    fm <- as.formula(paste("Metagene_mean", "~", "treatment")) #make formula to plot
    
    #plot and split by patient:   
    stat.test = compare_means(formula = fm ,data = genes_by_tp,method = test,group.by = patient.ident_var,p.adjust.method = "fdr")%>% # Add pairwise comparisons p-value
      dplyr::filter(group1 == pre_on[1] & group2 == pre_on[2])  #filter for pre vs on treatment only
    
    stat.test$p.format =stat.test$p.adj #modift 0 pvalue to be lowest possible float
    stat.test$p.format[!stat.test$p.format == 0 ] <- paste("=",stat.test$p.format[!stat.test$p.format == 0 ])
    stat.test$p.format[stat.test$p.format == 0 ] <- paste("<",.Machine$double.xmin %>% signif(digits = 3))
    
    plt = ggplot(genes_by_tp, aes(x = !!data_sym(patient.ident_var), y = Metagene_mean, fill = treatment)) + geom_split_violin(scale = 'width')+ylab(metegene)+ 
      geom_boxplot(width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0)+
      ylim(min(genes_by_tp$Metagene_mean),max(genes_by_tp$Metagene_mean)*1.25)
    plt = plt +stat_pvalue_manual(stat.test, label = "p {p.format}",  #add p value
                                  y.position = max(genes_by_tp$Metagene_mean)*1.08,x = patient.ident_var,inherit.aes = F,size = 3.3) # set position at the top value
    
    plt.lst[[metegene]] = plt
    if (!return_list) {
      print(plt)
    }
  }
  if (return_list) {
    return(plt.lst)
  }
}

metagenes_violin_combine_patients = function(dataset,prefix = "",pre_on = c("OSI","NT"),axis.text.x = 11,test = "t.test", programs = c("Hypoxia","TNFa","Cell_cycle"),return_list = F){
  require(facefuns)
  plt.lst = list()
    genes_by_tp = FetchData(object = dataset,vars =  c("treatment",programs)) %>% filter(treatment %in% pre_on)  %>% as.data.frame() #mean expression
    formula <- as.formula( paste("c(", paste(programs, collapse = ","), ")~ treatment ") )
    
    #plot and split by patient:   
    stat.test = compare_means(formula = formula ,data = genes_by_tp,method = test,p.adjust.method = "fdr")%>% # Add pairwise comparisons p-value
      dplyr::filter(group1 == pre_on[1] & group2 == pre_on[2])  #filter for pre vs on treatment only
    
    stat.test$p.format =stat.test$p.adj #modift 0 pvalue to be lowest possible float
    stat.test$p.format[!stat.test$p.format == 0 ] <- paste("=",stat.test$p.format[!stat.test$p.format == 0 ])
    stat.test$p.format[stat.test$p.format == 0 ] <- paste("<",.Machine$double.xmin %>% signif(digits = 3))
    
    
    genes_by_tp = reshape2::melt(genes_by_tp, id.vars = c("treatment"),value.name = "score")
    plt = ggplot(genes_by_tp, aes(x = variable, y = score,fill = treatment)) + geom_split_violin(scale = 'width')+ 
      geom_boxplot(width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0)+
      ylim(min(genes_by_tp$score),max(genes_by_tp$score)*1.25)
    plt = plt +stat_pvalue_manual(stat.test, label = "p {p.format}",  #add p value
                                  y.position = max(genes_by_tp$score)*1.08,inherit.aes = F,size = 3.3,x = ".y.") # set position at the top value
    return(plt)
}

# These functions are aim to calculate to usage, but with other count matrix. usage usually calculated with: usage = NMF(counts,gep_scores), so here we can make it with any 
# count matrix. This is good for the Harmony case, when we want to calculate gep with the corrected values, but infer the usage based on the original count matrix
# so we will be sure the Harmony correction is not creating bias on the usage. 

get_norm_counts = "def get_norm_counts(counts, tpm,high_variance_genes_filter): #from cnmf.py
      import numpy as np
      import scipy.sparse as sp
      norm_counts = counts[:, high_variance_genes_filter].copy()
      
      ## Scale genes to unit variance
      if sp.issparse(tpm.X):
        sc.pp.scale(norm_counts, zero_center=False)
      if np.isnan(norm_counts.X.data).sum() > 0:
        print('Warning NaNs in normalized counts matrix')                       
      else:
        norm_counts.X /= norm_counts.X.std(axis=0, ddof=1)
      if np.isnan(norm_counts.X).sum().sum() > 0:
        print('Warning NaNs in normalized counts matrix')                    
      
      
      ## Check for any cells that have 0 counts of the overdispersed genes
      zerocells = norm_counts.X.sum(axis=1)==0
      if zerocells.sum()>0:
        examples = norm_counts.obs.index[zerocells]
        print('Warning: %d cells have zero counts of overdispersed genes. E.g. %s' % (zerocells.sum(), examples[0]))
        print('Consensus step may not run when this is the case')
      
      return(norm_counts)"

#Calculate usage matrix like cNMF, with any expression matrix
get_usage_from_score = 
  "def get_usage_from_score(counts,tpm, genes,cnmf_obj,k, sumTo1 = True,do_norm_counts = True): #based on 'consensus' method
      import anndata as ad
      import scanpy as sc
      import numpy as np
      from sklearn.decomposition import non_negative_factorization
      import pandas as pd
      counts_adata = ad.AnnData(counts)
      tpm_adata = ad.AnnData(tpm)
      
      #get matrices
      if(do_norm_counts):
        norm_counts = get_norm_counts(counts=counts_adata,tpm=tpm_adata,high_variance_genes_filter=np.array(genes)) #norm counts like cnmf
      else:
        norm_counts = ad.AnnData(counts)
      spectra_original = get_median_spectra(cnmf_obj= cnmf_obj,k=k) #get score 
      
      # filter 
      spectra = spectra_original[spectra_original.columns.intersection(genes)] #remove genes not in @genes
      norm_counts = norm_counts[:, spectra.columns].copy() #remove genes not in spectra
      spectra = spectra.T.reindex(norm_counts.to_df().columns).T #reorder spectra genes like norm_counts
      
      # calculate usage
      norm_counts.X = norm_counts.X.astype(np.float32)
      usage_by_calc,_,_ = non_negative_factorization(X=norm_counts.X.astype(np.float64), H = spectra.values, update_H=False,n_components = k,max_iter=1000,init ='random')
      usage_by_calc = pd.DataFrame(usage_by_calc, index=counts.index, columns=spectra.index) #insert to df+add names
      
      #normalize
      if(sumTo1):
          usage_by_calc = usage_by_calc.div(usage_by_calc.sum(axis=1), axis=0) # sum rows to 1 and assign to main df
      
      # reorder
        # get original order
      original_norm_counts = sc.read(cnmf_obj.paths['normalized_counts'])
      usage_by_calc_original,_,_ = non_negative_factorization(X=original_norm_counts.X, H = spectra_original.values, update_H=False,n_components = k,max_iter=1000,init ='random')
      usage_by_calc_original = pd.DataFrame(usage_by_calc_original, index=original_norm_counts.obs.index, columns=spectra_original.index)  
      norm_original_usages =usage_by_calc_original.div(usage_by_calc_original.sum(axis=1), axis=0)      
      reorder = norm_original_usages.sum(axis=0).sort_values(ascending=False)
        #apply
      usage_by_calc = usage_by_calc.loc[:, reorder.index]
      return(usage_by_calc)"

compute_tpm = 
  "def compute_tpm(input_counts): #from cnmf.py
    tpm = input_counts.copy()
    sc.pp.normalize_per_cell(tpm.to_numpy(), counts_per_cell_after=1e6,copy=True)
    return(tpm)"

get_median_spectra = 
   "
import numpy as np
import pandas as pd
import os, errno
import datetime
import uuid
import itertools
import yaml
import subprocess
import scipy.sparse as sp

from scipy.spatial.distance import squareform
from sklearn.decomposition import non_negative_factorization
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.utils import sparsefuncs

from fastcluster import linkage
from scipy.cluster.hierarchy import leaves_list

import matplotlib.pyplot as plt

import scanpy as sc

def save_df_to_npz(obj, filename):
	np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)

def save_df_to_text(obj, filename):
	obj.to_csv(filename, sep='\t')

def load_df_from_npz(filename):
	with np.load(filename, allow_pickle=True) as f:
		obj = pd.DataFrame(**f)
	return obj


def get_median_spectra(cnmf_obj, k, density_threshold=0.5, local_neighborhood_size = 0.30,show_clustering = True, skip_density_and_return_after_stats = False, close_clustergram_fig=False, refit_usage=True):
	merged_spectra = load_df_from_npz(cnmf_obj.paths['merged_spectra']%k)
	norm_counts = sc.read(cnmf_obj.paths['normalized_counts'])
	density_threshold_str = str(density_threshold)
	if skip_density_and_return_after_stats:
		density_threshold_str = '2'
	density_threshold_repl = density_threshold_str.replace('.', '_')
	n_neighbors = int(local_neighborhood_size * merged_spectra.shape[0]/k)
	# Rescale topics such to length of 1.
	l2_spectra = (merged_spectra.T/np.sqrt((merged_spectra**2).sum(axis=1))).T
	if not skip_density_and_return_after_stats:
		# Compute the local density matrix (if not previously cached)
		topics_dist = None
		if os.path.isfile(cnmf_obj.paths['local_density_cache'] % k):
			local_density = load_df_from_npz(cnmf_obj.paths['local_density_cache'] % k)
		else:
			#   first find the full distance matrix
			topics_dist = euclidean_distances(l2_spectra.values)
			#   partition based on the first n neighbors
			partitioning_order  = np.argpartition(topics_dist, n_neighbors+1)[:, :n_neighbors+1]
			#   find the mean over those n_neighbors (excluding self, which has a distance of 0)
			distance_to_nearest_neighbors = topics_dist[np.arange(topics_dist.shape[0])[:, None], partitioning_order]
			local_density = pd.DataFrame(distance_to_nearest_neighbors.sum(1)/(n_neighbors),
										 columns=['local_density'],
										 index=l2_spectra.index)
			# save_df_to_npz(local_density, self.paths['local_density_cache'] % k)
			del(partitioning_order)
			del(distance_to_nearest_neighbors)
		density_filter = local_density.iloc[:, 0] < density_threshold
		l2_spectra = l2_spectra.loc[density_filter, :]
	kmeans_model = KMeans(n_clusters=k, n_init=10, random_state=1)
	kmeans_model.fit(l2_spectra)
	kmeans_cluster_labels = pd.Series(kmeans_model.labels_+1, index=l2_spectra.index)
	# Find median usage for each gene across cluster
	median_spectra = l2_spectra.groupby(kmeans_cluster_labels).median()
	# Normalize median spectra to probability distributions.
	median_spectra = (median_spectra.T/median_spectra.sum(1)).T
	return (median_spectra)"
add_python_funs <- function() {
	py_run_string(code = get_norm_counts) #load function to python
	py_run_string(code = get_usage_from_score) #load function to python
	py_run_string(code = compute_tpm) #load function to python
	py_run_string(code = get_median_spectra) #load function to python
}	


