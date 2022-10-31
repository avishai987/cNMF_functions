tpm <- function(counts, lengths) {
  rpk <- counts / lengths
  coef <- sum(rpk) / 1e6
  rpk/coef
}



create_data <- function() {
  #From Yotam code
  setwd("C:/Users/avishaiw.EKMD/Desktop/Data/ACC sc")
  acc.data <- read.delim("all.txt",skip=0,header=T,sep="\t",stringsAsFactors=F,row.names=1)
  rownames(acc.data)=make.unique(acc.data$gene_name)
  lengths = acc.data[,5]
  omitgenes=startsWith(rownames(acc.data),"MT-")|startsWith(rownames(acc.data),"ERCC-")
  acc.tpm=apply(acc.data[!omitgenes,7:dim(acc.data)[2]], 2, function(x) tpm(x, lengths[!omitgenes]) )
  acc.tpm = as.data.frame(acc.tpm)
  
  acc.data=acc.data[,7:ncol(acc.data)]
  cell.labels <- gsub(".gene_counts.tsv","",colnames(acc.data))
  cell.labels <- gsub(".sort.bam","",cell.labels)
  cell.labels <- gsub("_",".",cell.labels)
  cell.labels <- gsub(".S[0-9]*","",cell.labels)
  well=regmatches(cell.labels, regexpr("[A-H][0-9]*$", cell.labels))
  plate=substr(cell.labels,1,nchar(cell.labels)-nchar(well)-1)
  
  colnames(acc.data)=paste(plate,well,sep="_")
  colnames(acc.tpm)=colnames(acc.data)
  lst = list (data = acc.data,tpm = acc.tpm)
  return (lst)
}

get_common_genes <- function(dataset, genes, metadata_name = "patient.ident") {
  #search for subset of genes, that are at list in every cell
  patient_vector = as.vector(unique(dataset@meta.data[[metadata_name]]))
  
  #find genes that are most var and at least one cell of every patient
  for (patient in patient_vector){
    #filter according to cancer cells:
    # ACC_patient = subset(x = dataset, subset = metadata_name == patient) #not working
    expr <- FetchData(dataset, vars = metadata_name)
    ACC_patient <- dataset[, which(expr == patient)]
    acc_patient_raw_counts = as.data.frame(ACC_patient@assays[["RNA"]]@counts)
    acc_patient_raw_counts = acc_patient_raw_counts[rownames(acc_patient_raw_counts) %in% genes,] #take only rows that is in genes
    seurat_raw_count <- CreateSeuratObject(counts = acc_patient_raw_counts, project = "ACC", min.cells = 1, min.features = 1) #remove genes that not exits in at list one cell
    genes = rownames (seurat_raw_count) #update  genes
  }
  return(genes)
}

write_patient_counts <- function(dataset, acc_counts,genes,name, all_genes = F) {
  #acc_counts: input from fread
  #dataset: seurat object of all cells
  
  data = acc_counts
  patient_vector = as.vector(unique(dataset@meta.data[["patient.ident"]]))
  
  #make 1 line as rownames, and convert to vector
  data = as.data.frame(data)
  rownames(data) <- data[,1] #gene names are in the 1 col
  data[1] <-NULL #delete 1 col
  data = as.data.frame(t(data))
  
  
  if (Sys.info()["nodename"] == "AVISHAI"){
    setwd("G:/האחסון שלי/תואר שני/ACC/original ACC/Cell Reports Submission/cnmf")
    
  } else{
    setwd("C:/Users/avishaiw.EKMD/My Drive/תואר שני/ACC/original ACC/Cell Reports Submission/cnmf")
    
  }
  
  
  for (patient in patient_vector){
    #filter according to cancer cells:
    ACC_patient = subset(x = acc_cancercells, subset = patient.ident == patient)
    acc_patient_data = filter_matrix(data = data,rows = genes, cols = colnames(ACC_patient))
    acc_patient_data = as.data.frame(t(acc_patient_data)) #cNMF require cells X genes 
    
    #create one of these:
    if (name == "tpm" & all_genes == T){
      write.table(x = acc_patient_data, file = paste0("./acc patients tpm all genes/",patient,"_tpm.txt"),
                  sep = '\t',row.names=T,col.names=NA) 
    }
    
    if (name == "counts"  & all_genes == T){
      write.table(x = acc_patient_data, file = paste0("./acc patients raw counts all genes/",patient,"_raw_counts.txt"),
                  sep = '\t',row.names=T,col.names=NA)
    }
    
    
    if (name == "tpm"  & all_genes == F){
      write.table(x = acc_patient_data, file = paste0("./acc patients tpm/",patient,"_tpm.txt"),
                  sep = '\t',row.names=T,col.names=NA) 
    }
    
    if (name == "counts"  & all_genes == F){
      write.table(x = acc_patient_data, file = paste0("./acc patients raw counts/",patient,"_raw_counts.txt"),
                  sep = '\t',row.names=T,col.names=NA)
    }
    
  }
  
}


filter_matrix <- function(data,rows,cols, remove_zero = T) {
  #insert into seurat to modify gene names 
  data <- CreateSeuratObject(counts = data, project = "ACC", min.cells = 0, min.features = 0)
  data = as.data.frame(data@assays[["RNA"]]@counts)
  
  #filter according cols
  data = data[,colnames(data) %in% cols]
  
  #filter according to gene filtration that perform on original data:
  data = data[rownames(data) %in% rows,]
  
  if (remove_zero == T){
    data <- CreateSeuratObject(counts = data, project = "ACC", min.cells = 1, min.features = 1) #remove zero features or cells
    data = as.data.frame(data@assays[["RNA"]]@counts)
    
  }
  # ''/````
  return (data)
}

#' @title read cnmf
#' @description read the cNMF result
#' @param cNMF_k k that the cNMF was run with
#' @param patients which patients to read
#' @param sum_to_1 ,divide each col by the sum, so that the sum will be 1 (sort of normalization): Default: T
#' @param sum_rows_to_1 divide each row by the sum (performed for every patient), so that the sum will be 1, Default: F
#' @param filter_uncommon_genes use correlation only with genes that in the top 200 genes of at at least one patient, Default: F
#' @return result table- genes X programs_of_all_patients.i.e, if  k ==3, and length(patinets)==4, the result will be genes X 12(3 progrmas per patients)
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export 

read_cnmf <- function(cNMF_k, patients, sum_to_1 = T, sum_rows_to_1 = F, filter_uncommon_genes = F, gene_list = NULL,
                      units = "TPM", directory = "./",quiet = F) {
  #should setwd to the folder with the results and genes file
  genes = scan("./common_genes.txt", character(), quiet =  quiet)
  genes = genes[order(genes)]
  
  #choose k that ran on cNMF
  result = data.frame(gene_names = genes)
  
  for (name in patients){
    #load patient programs scores
    folder = paste0("./",name,"_cNMF")
    if (units == "TPM"){
      file  = paste0(name,"_cNMF.gene_spectra_tpm.k_",cNMF_k,".dt_0_1.txt")
    }else if(units == "Z_score") {
      file  = paste0(name,"_cNMF.gene_spectra_score.k_",cNMF_k,".dt_0_1.txt")
    }
    table = fread(paste0(folder,"/",file)) %>%  as.data.frame()
    table[1] <- NULL
    table = as.data.frame(t(table))
    table = table[rownames(table) %in% genes,]
    #paste patient and program to colnames
    for (program_num in 1:cNMF_k) {
      colnames(table)[program_num] = paste0(name," program",program_num)
    }
    #order by rownames
    table = table[order(rownames(table)),]
    
    
    
    for (program_num in 1:cNMF_k) {
      # if(sum_to_1 == T){  old 11.9
      # program_sum = sum(table[,program_num])
      # norm_program = table[,program_num]/program_sum 
      # }else { norm_program = table[,program_num]}
      norm_program = table[,program_num]
      
      if (table(rownames(table) == result$gene_names)["TRUE"] != length(genes)){
        print ("wrong calculation, check genes order")
      }
      result[,paste0(name," program",program_num)] = norm_program #paste in result table
    }
    
    
  }
  
  
  #make gene names as rownames:
  rownames(result) = result$gene_names
  result= result[ , -which(names(result) %in% c("gene_names")) ]
  
  
  if (filter_uncommon_genes == T) {
    all_common_genes = c()
    for (col_name in names(result)){
      col = result[,col_name,drop = F]
      col = col %>% ungroup() %>% arrange(desc(col[1]))
      top_genes = head(col,200) %>% rownames() # "num_of_top_genes" genes from top_genes
      all_common_genes = c(all_common_genes, top_genes)
    }
    all_common_genes = unique(all_common_genes)
    result = result[(row.names(result) %in% all_common_genes), ]
  }
  
  
  if ( is.null(gene_list) == F ) { #filter by gene list
    result = result %>% filter((row.names(result) %in% gene_list))
  }
  
  
  if (sum_rows_to_1 == T){ #sum rows to 1,instead of cols. Yotam said that this could be better and improve clustering of the programs
    for (row_num in 1:nrow(table)){
      row_sum = sum(result[row_num,])
      norm_row = result[row_num,]/row_sum 
      result[row_num,] = norm_row
      
    }
  }
  
  if(sum_to_1 == T){ #new 11.9
    for (col_num in 1:ncol(table)){
      program_sum = sum(result[,col_num])
      norm_program = result[,col_num]/program_sum 
      result[,col_num] = norm_program
    }
  }
  
  
  return (result)
}


cluster_programs <- function(num_of_clusters, pht, cluster_col_or_row = "tree_row") {
  myannotation = as.data.frame(cutree(pht[[cluster_col_or_row]], k = num_of_clusters))
  names(myannotation)[1] = "cluster"
  myannotation$cluster = as.factor(myannotation$cluster)
  
  palette1 <- palette(brewer.pal(num_of_clusters, "Paired"))
  palette1 <- palette(brewer.pal(num_of_clusters, "Paired")) #need to do twice to work
  
  names(palette1) = unique(myannotation$cluster)
  ann_colors = list (cluster = palette1)
  return(list(ann_colors = ann_colors, myannotation = myannotation))
}

#analyse cluster quality for each cluster- num of progrmas/ unique patients. 1 indicates that no patient has 2 programs in this cluster
analayze_clusters <- function(annotation) {
  all_clusters = annotation[["myannotation"]]
  all_cluster_indexes = 1:(all_clusters$cluster %>% as.vector() %>% max())
  # cluster_analysis = data.frame(cluster = paste("cluster", all_cluster_indexes))
  
  for (cluster_index in all_cluster_indexes){ #for 1 ,2 3...k clusters 
    programs = all_clusters %>% filter(cluster == cluster_index) %>% rownames()
    first_space = unlist(gregexpr(' ', programs))
    patients_names = programs
    str_sub(patients_names, first_space, stri_length(patients_names)) <- "" #remove from space to the end to get patient name
    patients_names = unique(patients_names)
    final_quality = length(programs)/length(patients_names)
    print("cluster " %s+% cluster_index  %s+% " quality: "%s+%final_quality)
  }
}


cnmf_enrichment <- function(num_of_clusters,annotation, result, num_of_top_genes, db, write_top_genes = F, normalization = "average",
                            top_genes_num = 200 ,common_gene_num = 2, write_plots = T, silent = F,
                            print_significant = T) {
  if (normalization != "average" & normalization != "common_genes"){
    stop("normalization != average & normalization != common_genes")
  }
  
  analayze_clusters(annotation = annotation)
  genes = scan("C:/Users/avishaiw.EKMD/My Drive/תואר שני/EGFR/Single cell project (Dana)/cNMF/common_genes.txt", character(), quote = "")
  myannotation = annotation[["myannotation"]]
  hypyxia_programs = c()
  result_lst = list()
  for (cluster in 1:num_of_clusters) {
    col_name = paste("program",cluster)
    
    #If cluster has only 1 program, skip enrichment analysis
    lst = myannotation$cluster == cluster
    sum = sum(lst, na.rm = TRUE)
    if (sum <= 2){next}
    
    program_rows = rownames(myannotation)[myannotation$cluster == cluster] #get rows of the cluster
    
    if (normalization == "average"){
      program_concensus =data.frame(score = rowMeans(result[,program_rows]))  #create df with mean of cluster
      top_genes_indexes = sort(program_concensus$score,index.return=TRUE, decreasing = T)$ix #sort by score and get top genes indexes
      top_genes = rownames(program_concensus)[head(top_genes_indexes, num_of_top_genes)] # "num_of_top_genes" genes from top_genes
      if (write_top_genes == T){
        fwrite(x = as.list(top_genes), file = paste("program",cluster,"top genes.txt" ),sep = "\n")
      }
    }
    
    if (normalization == "common_genes"){ #take genes to enrichment if that in at list 2 columns in the program 
      fun <- function(x, character = FALSE) {
        top_genes = sort(x, decreasing = T) #sort by score a
        top = head(names(top_genes),top_genes_num) #take top top_genes_num
        top #export
      }
      
      program = result[,program_rows] #take cols of the program
      all_genes = apply(program, 2, fun)  #get top genes for each column
      all_genes = all_genes %>% t() %>% as.vector() %>% table() #convert into 1 vector and make table
      top_genes  = names(all_genes)[all_genes > common_gene_num] #take genes that have more than common_gene_num columns for enrichment analysis
      if(length(top_genes) == 0){next} #if there are no genes that in common
    }
    
    #perform enrichment analysis
    enrichment_result = genes_vec_enrichment(genes = top_genes, background = genes, gene_sets = db, title = paste("cluster",cluster),
                                             add_bg = F, silent = silent)
    if (print_significant == TRUE){
      print ("cluster" %s+% cluster %s+% ":") #return only significant values
      print (enrichment_result %>% filter(enrichment_result$p.adjust <0.05)) #return only significant values
      
    }
    name = "cluster" %s+% cluster
    current_result = enrichment_result %>% filter(enrichment_result$p.adjust <0.05)
    result_lst [[name]] = current_result
    # write.csv(x = enrichment_res,file = paste("program",cluster,"enrichment analysis.csv" ))
    #check which programs are classified as hypoxia:
    # if (enrichment_res[enrichment_res$pathway_name == "HALLMARK_HYPOXIA", "p.adjust"] > 0.05){
    #   hypyxia_programs = c(hypyxia_programs, program_rows)
    # }
  }
  return (result_lst)
  
  # return(hypyxia_programs)
}

find_clusters <- function(result,num_of_clusters, clustering_distance, write = T) {
  cor_martix = cor (result)
  patients = colnames(cor_martix)
  patients = gsub(pattern = " program.*", replacement = "",x = patients)
  patient_annontation = data.frame(patient = patients, row.names = colnames(cor_martix))
  
  pht = pheatmap(mat = cor_martix, breaks=seq(-1, 1, length.out=101),
                 annotation_col = patient_annontation,
                 clustering_distance_rows = clustering_distance,
                 clustering_distance_cols = clustering_distance,silent = T)
  
  
  
  annotation = cluster_programs(num_of_clusters = num_of_clusters, pht = pht)
  
  p = pheatmap(mat = cor_martix, breaks=seq(-1, 1, length.out=101),annotation_col =  annotation[["myannotation"]],
               annotation_colors = annotation[["ann_colors"]], clustering_distance_rows = clustering_distance,
               clustering_distance_cols = clustering_distance)
  if(write == T){
    pdf(file = "correlation plot.pdf",width = 10,height = 8)
    print(p)
    dev.off()
  }
  print(p)
  return(annotation)
  
}

#Instead of average programs in each cluster, analyze every single one of them and then cluster
analyze_by_single_program  <- function(cNMF_k,patients_vector, db,sum_rows_to_1 , sum_to_1 , units,clustering_distance ) {
  setwd("C:/Users/avishaiw.EKMD/Desktop/Data/EGFR/cNMF/input/EGFR_cNMF/")
  all_table = data.frame(names = unique(db$gs_name), row.names = unique(db$gs_name))
  
  for (patient in patients_vector){
    result = read_cnmf(cNMF_k = cNMF_k,patients = patient, sum_rows_to_1 = sum_rows_to_1, sum_to_1 = sum_to_1,units = units)
    for (col_num in 1:ncol(result)){
      all_table[names(result)[col_num]] = 0
      col = result[,col_num,drop = F] 
      col = col %>% ungroup() %>% arrange(desc(col[1]))
      top_genes = head(col,200) %>% rownames() # "num_of_top_genes" genes from top_genes
      genes = scan("C:/Users/avishaiw.EKMD/My Drive/תואר שני/EGFR/Single cell project (Dana)/cNMF/common_genes.txt", character(), quote = "")
      enrichment_res = genes_vec_enrichment(genes = top_genes, background = genes, gene_sets = db, title = names(result)[col_num],
                                            add_bg = F, silent = T)
      enrichment_res = enrichment_res %>% arrange(factor(enrichment_res$pathway_name, levels = all_table$names)) #sort enrichment res by result table
      all_table[names(result)[col_num]] = enrichment_res$p.adjust
    }
  }
  all_table$names <- NULL
  p = sig_heatmap(all_patients_result = all_table,title = "all patients separately k =3",clustering_distance = clustering_distance)
  return (p)
}

add_prgorams_score <- function(result,num_of_clusters,annotation,cNMF_k,normalization,num_of_top_genes,
                               dataset = NULL) {
  myannotation = annotation[["myannotation"]]
  
  for (cluster in 1:num_of_clusters) {
    col_name = paste("metaProgram",cluster)
    
    program_rows = rownames(myannotation)[myannotation$cluster == cluster] #get rows of the cluster
    
    if (normalization == "average"){
      program_consensus =data.frame(score = rowMeans(result[,program_rows]))  #create df with mean of cluster
      program_consensus = program_consensus %>% arrange(desc(score)) %>% top_n(200) #sort by score and take top 200 rows
      program_consensus = program_consensus %>% rename(!!col_name := score) #rename "score" to col_name
    }else{ stop ("only average is available at the moment")}
    gene_expression = dataset@assays[["RNA"]]@data [ rownames(dataset@assays[["RNA"]]@data) %in% rownames(program_consensus),] %>% 
      as.data.frame()
    # Sort d1 and d2 with columns A and B
    gene_expression <- gene_expression[order(match(rownames(gene_expression),rownames(program_consensus))),]
    final_score = program_consensus[,1,drop = T] * gene_expression
    final_score_average = colMeans(final_score)
    final_score_average = data.frame(score = final_score_average, row.names = colnames(final_score))
    final_score_average = final_score_average %>% rename(!!col_name := score) #rename "score" to col_name
    
    dataset = AddMetaData(object = dataset,metadata = final_score_average)
  }
  return(dataset)
}
