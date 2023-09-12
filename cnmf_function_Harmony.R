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

expression_mult = function(gep_scores,dataset, top_genes = F,max_genes = F, z_score = F,min_max = F,sum2one = F, hallmark_genes  = NULL,top_genes_num = 200) {
  if (top_genes){ #for every metagene ,multiple only the top genes
    cell_usage = data.frame(row.names =colnames(dataset)) #create empty df to store results
    for (col_num in 1:ncol(gep_scores)) {
      top_200 = gep_scores %>% select(col_num) %>%  arrange(desc(gep_scores[col_num])) %>% head(top_genes_num)  #take top 200 rows
      if (!is_null(hallmark_genes) ){ #intersect with genes from hallmark, i.e take only hypoxia genes
        genes_in_hallmark = intersect(rownames(top_200),hallmark_genes[[col_num]])
        top_200 = top_200[rownames(top_200) %in% genes_in_hallmark,,drop = F]
      }
      top_200 = top_200 %>% t() %>%  as.matrix()
      expression = dataset@assays$RNA@data %>% as.matrix()
      expression = expression[rownames(expression) %in% colnames(top_200),,drop=F]  #remove rows not in top_genes
      top_200= top_200[,colnames(top_200) %in% rownames(expression),drop=F] #remove rows not in expression
      expression = expression[match(colnames(top_200), rownames(expression)),] #order expression rows like gep
      expression = 2**(expression) #convert from log(tpm+1) to tpm
      expression = expression-1
      
      my_usage = top_200%*%expression
      metagene = my_usage %>% t() %>% as.data.frame()
      cell_usage = cbind(cell_usage,metagene)
    }
    cell_usage = cell_usage %>% setNames(colnames(gep_scores)) 
    
  }else if(max_genes){
    require(NMF,quietly = T)
    top_features = extractFeatures(object = gep_scores %>% data.matrix(),method ="max")
    for (i in 1:length(top_features)) {
      top_features[[i]]= rownames(gep_scores)[top_features[[i]]]
    }
    
    cell_usage = data.frame(row.names =colnames(dataset)) #create empty df to store results
    for (i in 1:ncol(gep_scores)) {
      top = top_features[i] %>% unlist()
      expression = dataset@assays$RNA@data %>% as.matrix()
      top_df = gep_scores[rownames(gep_scores) %in% top,i,drop=F] %>% t() %>%  as.matrix()
      
      expression = expression[rownames(expression) %in% colnames(top_df),,drop=F]  #remove rows not in top_genes
      top_df= top_df[,colnames(top_df) %in% rownames(expression),drop=F] #remove rows not in expression
      
      expression = expression[match(colnames(top_df), rownames(expression)),] #order expression rows like gep
      my_usage = top_df%*%expression
      metagene = my_usage %>% t() %>% as.data.frame()
      cell_usage = cbind(cell_usage,metagene)
    }
    cell_usage = cell_usage %>% setNames(colnames(gep_scores)) 
    
    
  }else{
    gep_scores = gep_scores  %>% t() %>%  as.matrix()
    expression = dataset@assays$RNA@data %>% as.matrix()
    expression = expression[rownames(expression) %in% colnames(gep_scores),] #remove rows not in gep_scores
    gep_scores= gep_scores[,colnames(gep_scores) %in% rownames(expression)] #remove rows not in expression
    expression = expression[match(colnames(gep_scores), rownames(expression)),] #order expression rows like gep
    
    cell_usage = gep_scores%*%expression #multiply 
    cell_usage = cell_usage %>% t() %>% as.data.frame()
  }
  #normalize:
  if (z_score) {
    cell_usage = scale (cell_usage) %>% as.data.frame()
  }
  else if(min_max){
    cell_usage = apply(cell_usage, MARGIN = 2, FUN = min_max_normalize)%>% as.data.frame()
  }
  else if(sum2one){
    cell_usage = sum2one(cell_usage)
  }
  
  return(cell_usage)
}

#' create metagenes by expression and gep scores by:
#' metagenes * gep = expression -> metagenes = (gep)-1* expression
#' so (gep)-1 is the left inversion claculated by ginv
#' @param gep_scores 
#' @param dataset 
#'
#' @return
#' @export
#'
#' @examples
expression_inversion <- function(gep_scores,dataset) {
      expression = dataset@assays$RNA@data %>% as.matrix()
      expression = expression[rownames(expression) %in% rownames(gep_scores),,drop=F]  #remove rows not in top_genes
      gep_scores = gep_scores[rownames(gep_scores) %in% rownames(expression),,drop=F]  #remove rows not in top_genes
      expression = expression[match(rownames(gep_scores), rownames(expression)),] #order expression rows like gep
      left_inversion = MASS::ginv(gep_scores)
      res = left_inversion %*% expression
      return(res)
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



sum2one <- function(df) {
  for (col_num in 1:ncol(df)){
    program_sum = sum(df[,col_num])
    norm_program = df[,col_num]/program_sum
    df[,col_num] = norm_program
  }
      return(df)

}

min_max_normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

#take common genes from each score, and mean them
merge_scores_with_common <- function(score_1,score_2,normalization) {
  common_genes = intersect(rownames(score_1), rownames(score_2))
  
  score_1 = score_1[rownames(score_1) %in% common_genes,]
  score_2 = score_2[rownames(score_2) %in% common_genes,]
  
  if (normalization == "z_score") {
    score_1 = scale(score_1)
    score_2 = scale(score_2)
  }
  else if (normalization == "sum2one") {
    score_1 = sum2one(score_1)
    score_2 = sum2one(score_2)
  }
  
  else if (normalization == "min_max_normalize") {
    score_1 = apply(score_1, MARGIN = 2, FUN = min_max_normalize)%>% as.data.frame()
    score_2 = apply(score_2, MARGIN = 2, FUN = min_max_normalize)%>% as.data.frame()
  }
  else if(normalization == "none"){}

  score_2 = score_2[match(rownames(score_1), rownames(score_2)),] #order score_1 rows like score_2
  cor_res = cor(x = score_1,y = score_2)

  print(
    pheatmap(cor_res,breaks=seq(-1, 1, length.out=100),display_numbers = T,color = colorRampPalette(c("blue", "white", "red"))(100))
  )
  return(list(score_1,score_2))
}
#take common genes from each score, and mean them, and then add the not common genes from both scores
merge_scores_with_adding <- function(score_1,score_2,z_score = F,min_max = F,sum2one = F) {

  #scale before merging
  if (z_score) {
    score_1 = scale (score_1) %>% as.data.frame()
    score_2 = scale(score_2) %>% as.data.frame()
  }
  else if(min_max){
    score_1 = apply(score_1, MARGIN = 2, FUN = min_max_normalize)%>% as.data.frame()
    score_2 = apply(score_2, MARGIN = 2, FUN = min_max_normalize)%>% as.data.frame()
  }
  else if(sum2one){
    score_1 = sum2one(score_1)
    score_2 = sum2one(score_2)
  }

  
  #find common rows and set 2 scores alike 
  common_genes = intersect(rownames(score_1), rownames(score_2))
  message(paste("common genes:", length(common_genes)))
  common_score_1 = score_1[rownames(score_1) %in% common_genes,]
  common_score_2 = score_2[rownames(score_2) %in% common_genes,]
  common_score_2 = common_score_2[match(rownames(common_score_1), rownames(common_score_2)),] #order lung rows like xeno
  
  # mean and bind common rows
  cell_cycle = common_score_1[,3] %>% cbind(common_score_2[,3]) %>% rowMeans()
  on_treatment = common_score_1[,2] %>% cbind(common_score_2[,2]) %>% rowMeans()
  hypoxia = common_score_1[,1] %>% cbind(common_score_2[,1]) %>% rowMeans()
  unknown = common_score_1[,4] %>% cbind(common_score_2[,4]) %>% rowMeans()
  common_gep_scores = cbind(hypoxia,on_treatment,cell_cycle,unknown) %>% as.data.frame() 
  
  
  #add all other rows
  only_score_1 = score_1[!rownames(score_1) %in% common_genes,]
  only_score_2 = score_2[!rownames(score_2) %in% common_genes,]
  
  merged_gep_scores = rbind(only_score_1,only_score_2) %>%  setNames(c("cell_cycle","on_treatment","hypoxia","unknown")) %>% rbind(common_gep_scores)
  
  #Enrichment analysis 
  plt_list = list()
  for (i in 1:ncol(merged_gep_scores)) {
    top_genes = merged_gep_scores  %>%  arrange(desc(merged_gep_scores[i])) #sort by score a
    top = head(rownames(top_genes),200) #take top top_genes_num
    res = genes_vec_enrichment(genes = top,background = rownames(merged_gep_scores),homer = T,title = 
                                 names(merged_gep_scores)[i],silent = T,return_all = T)
    
    plt_list[[i]] = res$plt
  }
  gridExtra::grid.arrange(grobs = plt_list)
  return(merged_gep_scores)
}

#' @title combine variable genes of 2 datasets
#' @description take genes from dataset1 if there are in most @most_var_genes_num genes in dataset1 && in the @all_var_genes_num of dataset2. same for dataset 2.
#' @param dataset_1 PARAM_DESCRIPTION
#' @param dataset_2 PARAM_DESCRIPTION
#' @param all_var_genes_num PARAM_DESCRIPTION
#' @param most_var_genes_num PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @export 

combine_var_genes <- function(dataset_1, dataset_2, all_var_genes_num, most_var_genes_num_1,most_var_genes_num_2,plot = F) {
  dataset_1 = FindVariableFeatures(object = dataset_1,nfeatures = all_var_genes_num)
  dataset_1_all_vargenes = VariableFeatures(dataset_1)
  dataset_1 = FindVariableFeatures(object = dataset_1,nfeatures = most_var_genes_num_1)
  dataset_1_most_vargenes = VariableFeatures(dataset_1)
  
  dataset_2 = FindVariableFeatures(object = dataset_2,nfeatures = all_var_genes_num)
  dataset_2_all_vargenes = VariableFeatures(dataset_2)
  
  dataset_2 = FindVariableFeatures(object = dataset_2,nfeatures = most_var_genes_num_2)
  dataset_2_most_vargenes = VariableFeatures(dataset_2)
  
  dataset_2_most_vargenes = dataset_2_most_vargenes[dataset_2_most_vargenes %in% dataset_1_all_vargenes]
  dataset_1_most_vargenes = dataset_1_most_vargenes[dataset_1_most_vargenes %in% dataset_2_all_vargenes]
  
  genes_lst = c(dataset_2_most_vargenes,dataset_1_most_vargenes)
  grid.newpage()
  draw.pairwise.venn(area1=length(dataset_1_most_vargenes), area2=length(dataset_2_most_vargenes),cross.area=intersect(dataset_1_most_vargenes,dataset_2_most_vargenes) %>% length(),
                     category=c("xeno_most_vargenes","patients_most_vargenes"),fill=c("Red","Yellow"))
  
  genes_lst = unique(genes_lst) #remove duplicates
  return(genes_lst)
}




#' @title find positive genes
#' @description find genes that have at least num_of_cells positive value after z score scailing
#' @param dataset dataset
#' @param num_of_cells num_of_cells to b positive
#' @return genes list
#' @export 

positive_genes <- function(dataset,num_of_cells) {
  genes_lst = c()
  for (row in 1:nrow(dataset@assays[["RNA"]]@scale.data)) {
    vector = dataset@assays[["RNA"]]@scale.data[row,]
    if (sum(vector>0,na.rm = T)>=num_of_cells){
      genes_lst = c(genes_lst,rownames(dataset@assays[["RNA"]]@scale.data)[row])
    }
  }
  length(genes_lst)
  return(genes_lst)
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
metagenes_mean_compare <- function(dataset,time.point_var,prefix = "",patient.ident_var,pre_on = c("OSI","NT"),axis.text.x = 11,test = "t.test", programs = c("Hypoxia","TNFa","Cell_cycle")) {
  
  for (metegene in programs) {
    #create data:
    genes_by_tp = FetchData(object = dataset,vars = metegene) %>% rowSums() %>% as.data.frame() #mean expression
    names(genes_by_tp)[1] = "Metagene_mean"
    genes_by_tp = cbind(genes_by_tp,FetchData(object = dataset,vars = c(patient.ident_var,time.point_var))) # add id and time points
    
    
    genes_by_tp_forPlot =  genes_by_tp %>% mutate(!!ensym(patient.ident_var) := paste(prefix,genes_by_tp[,patient.ident_var])) #add "model" before  each model/patient
    fm <- as.formula(paste("Metagene_mean", "~", time.point_var)) #make formula to plot
    
    #plot and split by patient:   
    stat.test = compare_means(formula = fm ,data = genes_by_tp_forPlot,method = test,group.by = patient.ident_var)%>% # Add pairwise comparisons p-value
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
"def get_usage_from_score(counts,tpm, genes,cnmf_obj,k, sumTo1 = True): #based on 'consensus' method
      import anndata as ad
      import scanpy as sc
      import numpy as np
      from sklearn.decomposition import non_negative_factorization
      import pandas as pd
      counts_adata = ad.AnnData(counts)
      tpm_adata = ad.AnnData(tpm)
      
      #get matrices
      norm_counts = get_norm_counts(counts=counts_adata,tpm=tpm_adata,high_variance_genes_filter=np.array(genes)) #norm counts like cnmf
      spectra_original = cnmf_obj.get_median_spectra(k=k) #get score 
      
      # filter 
      spectra = spectra_original[spectra_original.columns.intersection(genes)] #remove genes not in @genes
      norm_counts = norm_counts[:, spectra.columns].copy() #remove genes not in spectra
      spectra = spectra.T.reindex(norm_counts.to_df().columns).T #reorder spectra genes like norm_counts
      
      # calculate usage
      usage_by_calc,_,_ = non_negative_factorization(X=norm_counts.X, H = spectra.values, update_H=False,n_components = k,max_iter=1000,init ='random')
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
    sc.pp.normalize_per_cell(tpm, counts_per_cell_after=1e6,copy=True)
    return(tpm)"

py_run_string(code = get_norm_counts) #load function to python
py_run_string(code = get_usage_from_score) #load function to python
py_run_string(code = compute_tpm) #load function to python

#add this to cnmf.py to get median spectra score from cnmf object
# def get_median_spectra(self, k, density_threshold=0.5, local_neighborhood_size = 0.30,show_clustering = True,
#                        skip_density_and_return_after_stats = False, close_clustergram_fig=False,
#                        refit_usage=True):
#   merged_spectra = load_df_from_npz(self.paths['merged_spectra']%k)
# norm_counts = sc.read(self.paths['normalized_counts'])
# density_threshold_str = str(density_threshold)
# if skip_density_and_return_after_stats:
#   density_threshold_str = '2'
# density_threshold_repl = density_threshold_str.replace('.', '_')
# n_neighbors = int(local_neighborhood_size * merged_spectra.shape[0]/k)
# 
# # Rescale topics such to length of 1.
# l2_spectra = (merged_spectra.T/np.sqrt((merged_spectra**2).sum(axis=1))).T
# 
# if not skip_density_and_return_after_stats:
#   # Compute the local density matrix (if not previously cached)
#   topics_dist = None
# if os.path.isfile(self.paths['local_density_cache'] % k):
#   local_density = load_df_from_npz(self.paths['local_density_cache'] % k)
# else:
#   #   first find the full distance matrix
#   topics_dist = euclidean_distances(l2_spectra.values)
# #   partition based on the first n neighbors
# partitioning_order  = np.argpartition(topics_dist, n_neighbors+1)[:, :n_neighbors+1]
# #   find the mean over those n_neighbors (excluding self, which has a distance of 0)
# distance_to_nearest_neighbors = topics_dist[np.arange(topics_dist.shape[0])[:, None], partitioning_order]
# local_density = pd.DataFrame(distance_to_nearest_neighbors.sum(1)/(n_neighbors),
#                              columns=['local_density'],
#                              index=l2_spectra.index)
# save_df_to_npz(local_density, self.paths['local_density_cache'] % k)
# del(partitioning_order)
# del(distance_to_nearest_neighbors)
# 
# density_filter = local_density.iloc[:, 0] < density_threshold
# l2_spectra = l2_spectra.loc[density_filter, :]
# 
# kmeans_model = KMeans(n_clusters=k, n_init=10, random_state=1)
# kmeans_model.fit(l2_spectra)
# kmeans_cluster_labels = pd.Series(kmeans_model.labels_+1, index=l2_spectra.index)
# 
# # Find median usage for each gene across cluster
# median_spectra = l2_spectra.groupby(kmeans_cluster_labels).median()
# 
# # Normalize median spectra to probability distributions.
# median_spectra = (median_spectra.T/median_spectra.sum(1)).T
# return (median_spectra)
