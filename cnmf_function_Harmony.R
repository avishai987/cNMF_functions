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
expression_mult<-function(gep_scores,dataset, top_genes = F,max_genes = F, z_score = F,min_max = F,sum2one = F) {
  if (top_genes){ #for every metagene ,multiple only the top genes
    cell_usage = data.frame(row.names =colnames(dataset)) #create empty df to store results
    for (col_num in 1:ncol(gep_scores)) {
      top_200 = gep_scores %>% select(col_num) %>%  arrange(desc(gep_scores[col_num])) %>% head(200)  #take top 200 rows
      top_200 = top_200 %>% t() %>%  as.matrix()
      expression = dataset@assays$RNA@data %>% as.matrix()
      expression = expression[rownames(expression) %in% colnames(top_200),,drop=F]  #remove rows not in top_genes
      top_200= top_200[,colnames(top_200) %in% rownames(expression),drop=F] #remove rows not in expression
      expression = expression[match(colnames(top_200), rownames(expression)),] #order expression rows like gep
      expression = 2**expression #convert from log(tpm+1) to tpm
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
cell_percentage <- function(dataset,time.point_var) {
    data =FetchData(object = dataset,vars = c("program.assignment",time.point_var))
  data = data %>% dplyr::count(program.assignment, .[time.point_var]) %>%  dplyr::add_count(program.assignment, wt = n, name = "overall")%>% 
  mutate(proportion = n / overall)   
  
  plt_list = list()
    time.point_var = ensym(time.point_var)
  for (program_name in unique(data$program.assignment)) {
    program_data = data[data$program.assignment == program_name,]
    p = ggplot(data=program_data, aes(x=!!time.point_var, y=proportion)) +geom_bar(stat="identity")+ylab("precentage") +
      ggtitle("program" %>% paste(program_data$program.assignment %>% unique() %>% as.character()))+
      scale_y_continuous(limits = c(0,1))   
      plt_list[[program_name]] = p
  }
    p = ggarrange(plotlist = plt_list )
    return(p)
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

#union programs (i.e all cell cycle programs). input example: groups_list = list(c(2,6),c(5,3,4),c(1))
# will combine metagenes 2+6, 5+3+4,1 in the df @all_metagenes
union_programs <- function(groups_list,all_metagenes) {
  unioned_metagenes = data.frame(row.names = rownames(all_metagenes)) #create final df 
  for (group in groups_list) {
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