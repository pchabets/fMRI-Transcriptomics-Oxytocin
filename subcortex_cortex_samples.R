########################################################################################################################################
## Author: P.C. Habets ## Department: Internal Medicine, Division of Endocrinology  ####
####################################################################################################################################
## Institution: Leiden Medical Center Utrecht ## Script: select only cortical and subcortical AHBA samples and write to file #####
####################################################################################################################################

# ontology file is the same for each donor, showing structural paths annotated with distinct numbers: 

# 4008 - cortical samples
# 4275 - subcortical samples
# 4392 - thalamus samples
# 4249 - hippocampal samples


library(tidyverse)

setwd("/path/to/local/folder/data") #set working directory to local folder containing the AHBA downloaded files. See readme file for setup of local folder with files.

#read in ontology file (same for all donors)
ontology <- read.csv("normalized_microarray_donor01/Ontology.csv")

##AHBA donor brain sample info
donorNames <- c("donor01", "donor02", "donor03", "donor04", "donor05", "donor06")
dnr <- lapply(donorNames, function(x){
  path <- paste0("normalized_microarray_", x, "/SampleAnnot.csv")
  dn <- read.csv(path, header = T)
  dn$donor <- x
  d <- paste0(x, "_sample_%d")
  dn$sample <- sprintf(d, seq(1:nrow(dn)))
  dn
})
dnr <- bind_rows(dnr)
dnr$donor <- factor(dnr$donor)
dnr$sample <- factor(dnr$sample)

#select only ontology IDs that have 4008 in structural path, but not 4275, 4392 or 4249
exclude <- paste(c("4275", "4392", "4249"),collapse = '|')
corticalID <- ontology %>% 
  filter(str_detect(structure_id_path, "4008")) %>% 
  filter(str_detect(structure_id_path, exclude, negate = TRUE)) %>% 
  dplyr::rename(structure_id = id) %>% 
  dplyr::rename(structure_acronym = acronym)

#remove cerebellum (CB) and brainstem (BS) samples from concatenated sample dataframe
dnrCX <- dnr %>% 
  filter(slab_type == "CX") 

#split remaining samples into cerebral cortical samples versus subcortical samples
dnrCortical <- inner_join(dnrCX, corticalID)
dnrSubCortical <- anti_join(dnrCX, corticalID)

#save splitted samples to files
write_csv(dnrCortical, "corticalSamplesOnly.csv")
write_csv(dnrSubCortical, "SubCorticalSamplesOnly.csv")





