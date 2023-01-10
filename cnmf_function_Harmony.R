#assign program to each cell after add_prgorams_score
#' @param larger_by how much a program should be greater then the second one to be assign as this program of this cell.
program_assignment <- function(dataset,larger_by = 1,program_names) {
  
  scores_metadata = FetchData(object = dataset,vars = c(program_names))
  assignment_df =  data.frame(row.names = rownames(scores_metadata))
  
  find_max_program <- function(x,larger_by) { #find max program index, else: "NA"
    n <- length(x)
    max = sort(x)[n] #get max value
    second_max = sort(x)[n-1] #get 2nd max value
    if (max > larger_by*second_max){which.max(x)} else {"NA"} #if max is larger than 2nd max*larger_by, write max index. else, write "NA"
  }
  assignment_df["assigment"] = apply(X = scores_metadata,MARGIN = 1,FUN = find_max_program,larger_by = larger_by) 
  
  dataset = AddMetaData(object = dataset,metadata = assignment_df,col.name = "program.assignment")
  return(dataset)
}
expression_mult <- function(gep_scores,dataset, top_genes = F,z_score = F,min_max = F,sum2one = F) {
  if (top_genes){ #for every metagene ,multiple only the top genes
    cell_usage = data.frame(row.names =colnames(dataset)) #create empty df to store results
    for (col_num in 1:ncol(gep_scores)) {
       top_200 = gep_scores %>% select(col_num) %>%  arrange(desc(gep_scores[col_num])) %>% head(200)  #take top 200 rows
       top_200 = top_200 %>% t() %>%  as.matrix()
       expression = lung@assays$RNA@data %>% as.matrix()
       expression = expression[rownames(expression) %in% colnames(top_200),]  #remove rows not in top_genes
       top_200= top_200[,colnames(top_200) %in% rownames(expression)] #remove rows not in expression
      
        my_usage = top_200%*%expression
        metagene = my_usage %>% t() %>% as.data.frame()
        cell_usage = cbind(cell_usage,metagene)
    }
  
  }else{
    gep_scores = gep_scores  %>% t() %>%  as.matrix()
    expression = dataset@assays$RNA@data %>% as.matrix()
    expression = expression[rownames(expression) %in% colnames(gep_scores),] #remove rows not in gep_scores
    gep_scores= gep_scores[,colnames(gep_scores) %in% rownames(expression)] #remove rows not in expression
    
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
  data = data %>% dplyr::count(program.assignment, .[time.point_var]) %>%  dplyr::add_count(.[time.point_var], wt = n, name = "overall")%>% 
  mutate(proportion = n / overall)   
  
  plt_list = list()
    time.point_var = ensym(time.point_var)
  for (program_name in unique(data$program.assignment)) {
    program_data = data[data$program.assignment == program_name,]
    p = ggplot(data=program_data, aes(x=!!time.point_var, y=proportion)) +geom_bar(stat="identity")+ylab("precentage") +ggtitle("program" %>% paste(program_data$program.assignment %>% unique() %>% as.character())) 
      plt_list[[program_name]] = p
  }
  gridExtra::grid.arrange(grobs = plt_list)
}

sum2one <- function(df) {
  for (col_num in 1:ncol(df)){
    program_sum = sum(df[,col_num])
    norm_program = df[,col_num]/program_sum
    df[,col_num] = norm_program
    return(df)
  }
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