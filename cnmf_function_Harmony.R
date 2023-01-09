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
expression_mult <- function(gep_scores,dataset) {
gep_scores = gep_scores  %>% t() %>%  as.matrix()
expression = dataset@assays$RNA@data %>% as.matrix()
expression = expression[rownames(expression) %in% colnames(gep_scores),]
usage = gep_scores%*%expression
all_metagenes = usage %>% t() %>% as.data.frame()
return(all_metagenes)
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