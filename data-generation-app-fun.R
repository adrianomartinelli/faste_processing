#This script is used to create artificial FASTA sequences to test the extract counts C++ program
write_code_lists = function(code_list, file){
  for(i in code_list){
    write(i, file, append = T)
  }
}

#Create code
create_code = function(len, set){
  return(paste0(sample(set, len, replace = T), collapse = ''))
}

#Create code list
create_code_list = function(len_list = 10, len_code = 6, set = c("A", "T", "G", "C")){
  stopifnot(length(set)^len_code > len_list)
  
  res = vector("character", length = len_list)
  for(i in 1:len_list){
    repeat{
      res[i] = create_code(len_code, set) 
      
      if(!any(res[i] == res[-i])){
        break
      }
    }
  }
  return(res)
}

#Create fasta line
create_fasta_line = function(S_list, B_list, C_list, S_idx, B_idx){
  seq = S_list[[1]][S_idx[1]]
  
  #Append all building block codes separated by a const region
  for(i in 1:length(B_list)){
    seq = paste0(seq, C_list[i], B_list[[i]][B_idx[i]])
  }
  
  #Append remaining selection codes
  l = length(B_list)
  if(length(S_list) > 1){
    for(i in 2:length(S_list)){
      seq = paste0(seq, C_list[l+i-1], S_list[[i]][S_idx[i]])
    }
  }
  return(seq)
  #seq = paste0(seq, C_list[length(C_list)])
}

#Create fasta file
create_fasta_file = function(file, file_size,
                             set_code, 
                             set_const, 
                             len_S, len_S_list, n_S_list, 
                             len_B, len_B_list, n_B_list, 
                             len_C, len_C_list, n_C_list){
  
  info_line = ">This is random information"
  S_combinations = vector("list", file_size)
  B_combinations = vector("list", file_size)
  
  #start = seq(1, by = floor((seq_length-len_code*n_codes_lists)/(n_codes_lists-1)), length.out = n_codes_lists)
  S_list = lapply(vector("list", length = n_S_list), function(x) create_code_list(len_S_list, len_S, set_code))
  B_list = lapply(vector("list", length = n_B_list), function(x) create_code_list(len_B_list, len_B, set_code))
  C_list = create_code_list(n_C_list, len_C, set = set_const)
  
  S_names = paste0("S_Code", c(1:n_S_list), ".txt")
  B_names = paste0("B_Code", c(1:n_B_list), ".txt")
  C_names = paste0("Const", c(1:n_C_list), ".txt")
  
  for(i in 1:length(S_names)){
    write_code_lists(S_list[[i]], S_names[i])
  }
  
  for(i in 1:length(B_names)){
    write_code_lists(B_list[[i]], B_names[i])
  }
  
  for(i in 1:length(C_names)){
    write(C_list[i], C_names[i])
  }
  
  tbl = matrix(nrow = file_size, ncol = n_S_list + n_B_list + 1)
  for(i in 1:file_size){
    if(!(i %% 10^4)){
      print(i)
    }
    B_idx = sample(len_B_list, n_B_list, replace = T)
    S_idx = sample(len_S_list, n_S_list, replace = T)
    tbl[i,] = c(1,S_idx, B_idx)
    
    line = create_fasta_line(S_list, B_list, C_list, S_idx, B_idx)
    write(info_line, file, append=TRUE)
    write(line, file, append=TRUE)
  }
  return(tbl)
}

create_struct_file = function(struct_file, faste_file,
                              len_S, len_B, len_C,
                              n_S_list, n_B_list, n_C_list){
  
  S_start = vector("numeric", n_S_list)
  S_start[1] = 1
  if(n_S_list > 1){
    for(i in 2:n_S_list){
      S_start[i] = len_S + len_B*n_B_list + len_C*n_B_list + (i-1)*len_C + (i-2)*len_S + 1
    }
  }
  S_end = S_start + len_S - 1
  
  B_start = vector("numeric", n_B_list)
  
  for(i in 1:n_B_list){
    B_start[i] = len_S + (i)*len_C + (i-1)*len_B + 1
  }
  B_end = B_start + len_B - 1
  
  C_start = vector("numeric", n_C_list)
  for(i in 1:n_C_list){
    #first const region after first S_identifier, then follow at max n_B_list identifier, finaly all remaining seleciton identifier
    C_start[i] = len_S + min(n_B_list,(i-1))*len_B + (i-1)*len_C + max(i-n_B_list-1,0)*len_S +1
  }
  C_end = C_start + len_C - 1
  
  write(faste_file, struct_file, append = F)
  write("Evaluation/", struct_file, append = T)
  
  for(i in 1:n_C_list){
    line = paste(C_start[i], C_end[i], 'C', paste0("Const", i, ".txt"), sep = '\t')
    write(line, struct_file, append = T)
  }
  
  for(i in 1:n_S_list){
    line = paste(S_start[i], S_end[i], 'S', paste0("S_Code", i, ".txt"), sep = '\t')
    write(line, struct_file, append = T)
  }
  
  for(i in 1:n_B_list){
    line = paste(B_start[i], B_end[i], 'B', paste0("B_Code", i, ".txt"), sep = '\t')
    write(line, struct_file, append = T)
  }
}

create_summary = function(tbl){
  tbl = as_tibble(tbl)
  
  #Sum up all lines with the same building blocks and selection identifiers
  tbl_sum = group_by_at(tbl, 2:ncol(tbl)) %>% summarise(count = sum(V1))
  
  #Group data by selections
  tbl_sel = group_by_at(tbl_sum, 1:(n_S_list)) %>% group_data()
  
  #Reorder tbl_sum and drop selection identifiers
  #col_idx_B = (n_S_list+1):(ncol(tbl_sum)-1)
  tbl_sum2 = ungroup(tbl_sum) %>% select(-c(1:n_S_list)) %>% select(last_col(), everything())
  colnames(tbl_sum2) = c("Counts", paste0("Code", 1:(ncol(tbl_sum2)-1)))
  
  #Print selection-wise
  for(i in 1:nrow(tbl_sel)){
    rows = unlist(tbl_sel[i,]$.rows)
    sel = select(tbl_sel[i,], -last_col())
    sel = paste0(sel, collapse = "_")
    sel = paste0("./Data/selection_", sel, "_.txt")
    write.table(tbl_sum2[rows,], sel, row.names = FALSE, sep = "\t", quote = FALSE)
  }
}