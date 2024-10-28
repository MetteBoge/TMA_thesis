library(tidyverse)

data_path <- "/home/people/metped/data/TCGA_pancan_expected_count.tsv"

exp_count <- read_tsv(data_path, 
         show_col_types = FALSE) %>% 
  mutate(across(-c(sample), function(x) 2^x)) # Reverse log2

lib_size <- exp_count %>% 
  select(-sample) %>% 
  colSums()

lib_size <- lib_size %>% 
  as.data.frame() %>% 
  rownames_to_column() 

colnames(lib_size) <- c("sample", "lib.size")

write_delim(lib_size,
          "/home/people/metped/data/lib_size.txt",
          delim = "\t")


lib_size %>% 
  select(lib.size) %>% 
  min()
