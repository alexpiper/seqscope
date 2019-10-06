
#Process seqtab to simplify
seqtab <- readRDS("seqscope/data/dros_data/dros_seqtab.rds")

library(tidyverse)

seqtab_filt <- seqtab %>%
  as.data.frame() %>%
  mutate(SampleID = rownames(.)) %>%
  filter(str_detect(SampleID, "_Rep1"))%>%
  mutate(SampleID = str_replace(SampleID, pattern="_Rep1", replacement = "")) %>%
  filter(!str_detect(SampleID, "SynMock")) %>%
  magrittr::set_rownames(.$SampleID) %>%
  select(-SampleID) %>%
  select(which(colSums(.) > 0)) %>%
  as.matrix()

write_rds(seqtab_filt, path="seqscope/data/demo_seqtab.rds")

#Process samdf
samdf <- read.csv(file="seqscope/data/dros_data/dros_samdf.csv")

samdf_filt <- samdf %>%
  filter(str_detect(SampleID, "_Rep1"))%>%
  mutate(SampleID = str_replace(SampleID, pattern="_Rep1", replacement = "")) %>%
  mutate(sample_id = str_replace(sample_id, pattern="_Rep1", replacement = "")) %>%
  filter(!str_detect(SampleID, "SynMock"))

write_csv(samdf_filt, path="seqscope/data/demo_samdf.csv")

#Process taxtab to simplify
taxtab <- readRDS("seqscope/data/dros_data/dros_taxtab.rds")

library(tidyverse)

taxtab_filt <- taxtab %>%
  as.data.frame() %>%
  mutate(seq = rownames(.)) %>%
  filter(seq %in% colnames(seqtab_filt))%>%
  magrittr::set_rownames(.$seq)  %>%
  select(-seq) %>%
  as.matrix()

write_rds(taxtab_filt, path="seqscope/data/demo_taxtab.rds")
