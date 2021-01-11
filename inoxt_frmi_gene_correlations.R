########################################################################################################################################
## Author: P.C. Habets, C Mclain ## Department: Internal Medicine, Division of Endocrinology  ####
####################################################################################################################################
## Institution: Leiden Medical Center Utrecht ## Script: differential gene expression analysis using four different masks #####
####################################################################################################################################

library(dplyr)
library(oro.nifti)
library(parallel)
library(pracma)
library(neurobase)
library(plotly)
library(AnalyzeFMRI)
library(oce)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(limma)
library(Hmisc)
library(psych)

setwd("/path/to/local/folder/data") #set working directory to local folder containing the AHBA downloaded files. See readme file for setup of local folder with separate 'data' directory that holds all needed files.
cores <- detectCores()

######################################################################################################
## Step 1: remove batch effects across all 6 brains separately for cortical and subcortical samples
######################################################################################################

# select highest expression probes (selected based on the output of line 68 of the "inoxt_fmri_gene_expression.R" script: View(exprmeans)
probes <- c("1053419", "1053417", "1058445", "1013120", "1058917", "1058916")

cort <- read.csv("corticalSamplesonly.csv", header = TRUE) # all cortical samples
cortical <- as.vector(unique(cort[,5])) # cortical acronyms

donors <- c(1, 2, 3, 4, 5, 6) 
# get indexes of cortical samples
cort_ind <- lapply(donors, function(d) {
  # load data file
  FileSample <- paste0("normalized_microarray_donor0", d, "/SampleAnnot.csv")
  sample_info <- as.data.frame(read.csv(FileSample, header = TRUE))
  strucac <- sample_info$structure_acronym
  indac <- vector() # to store indexes for each acronym
  ind <- vector() # to store all indexes
  for (i in 1:length(cortical)) {
    indac <- which(strucac==cortical[i]) 
    ind <- append(ind, indac)
  }
  ind
})

# cortical  expression for each gene of interest 
cort_expr <- mclapply(donors, mc.cores=cores, function(d) {
  FileMicroarray <- paste0("normalized_microarray_donor0", d, "/MicroarrayExpression.csv")
  e <- read.csv(FileMicroarray, header = FALSE)
  rownames(e) <- e[,1]
  e <- e[,-1] # get rid of probe IDs
  ind <- cort_ind[[d]] 
  e <- e[probes,ind] # only columns corresponding to cortical samples of genes of interest
  e[7,] <- d 
  e
})
cort_expression <- Reduce(cbind, cort_expr)
limmacort <- data.frame(t(removeBatchEffect(cort_expression, batch = cort_expression[7,])))

# find well_id for each cortical sample to enable inner_join later
cort_annot <- lapply(donors, function(d) { sample_annot <- read.csv(paste0("normalized_microarray_donor0", d, "/SampleAnnot.csv"), header = TRUE)
sample_annot <- sample_annot[cort_ind[[d]],] 
sample_annot$well_id })

well_id <- c(unlist(cort_annot)) 
limmacort$well_id <- well_id
colnames(limmacort) <- c("OXT", "OXTR", "CD38", "AVPR2", "AVPR1A", "AVPR1B", "donor", "well_id")

# subcortical expression for each gene of interest 
subcort_expr <- mclapply(donors, mc.cores=cores, function(d) {
  FileMicroarray <- paste0("normalized_microarray_donor0", d, "/MicroarrayExpression.csv")
  e <- read.csv(FileMicroarray, header = FALSE)
  rownames(e) <- e[,1]
  e <- e[,-1] # get rid of probe IDs
  # remove cortical, brain stem and cerebellum samples
  indcort <- cort_ind[[d]] 
  sample_annot <- read.csv(paste0("normalized_microarray_donor0", d, "/SampleAnnot.csv"), header = TRUE)
  indBSCB <- c(which(sample_annot$slab_type == "BS"), which(sample_annot$slab_type =="CB"))
  indthrow <- c(indcort, indBSCB)
  e <- e[probes, -indthrow]
  e[7,] <- d
  e
})

subcort_expression <- Reduce(cbind, subcort_expr)
limmasubcort <- data.frame(t(removeBatchEffect(subcort_expression, batch = subcort_expression[7,])))

# find well_id for each sample to enable inner_join later
subcort_annot <- lapply(donors, function(d) { sample_annot <- read.csv(paste0("normalized_microarray_donor0", d, "/SampleAnnot.csv"), header = TRUE)
# remove cortical, brain stem and cerebellum samples
indcort <- cort_ind[[d]] 
sample_annot <- read.csv(paste0("normalized_microarray_donor0", d, "/SampleAnnot.csv"), header = TRUE)
indBSCB <- c(which(sample_annot$slab_type == "BS"), which(sample_annot$slab_type =="CB"))
indthrow <- c(indcort, indBSCB)
sample_annot <- sample_annot[-indthrow,]
sample_annot$well_id
}) 

well_id <- c(unlist(subcort_annot))
limmasubcort$well_id <- well_id
colnames(limmasubcort) <- c("OXT", "OXTR", "CD38", "AVPR2", "AVPR1A", "AVPR1B", "donor", "well_id")

#####################################################################
## Step 2: assign samples to affected/unaffected areas on the fMRI map
#####################################################################

##AHBA donor brains
donorNames <- c("donor01", "donor02", "donor03", "donor04","donor05", "donor06")
donor <- lapply(donorNames, function(x){
  path <- paste0("normalized_microarray_", x, "/SampleAnnot.csv")
  dn <- read.csv(path, header = T)
  dn$mask <- "no"
  dn$donor <- x
  dn
})
dnr <- bind_rows(donor)
dnr$donor <- factor(dnr$donor)
dnr$mask <- factor(dnr$mask)

# read in overall fMRI file
allTasks_mask <- "OXT-PBO_ALL_P.nii" # We used the all-task mask for correlation analysis
z <- readnii(allTasks_mask) 

dz <- img_color_df(z, zlim = NULL, breaks = NULL, col = gray(0:64/64))
dz <- dz[dz$value < 0.05,]
tB <- apply(dz[,1:3], 1, function(x){
  translateCoordinate(c(x[1], x[2], x[3]),z , verbose = F)
})
tB <- as.data.frame(t(tB))
tB <- tB %>% dplyr::distinct(V1, V2, V3, .keep_all = T)
tB$mask <- "yes"
tB$donor <- "mask"
colnames(tB) <- colnames(dnr[,11:15])

##plotting
plot <- bind_rows(dnr[,11:15], tB)
plot$mask <- factor(plot$mask)
plot$donor <- factor(plot$donor)
plot_ly(plot, x=~mni_x, y=~mni_y, z=~mni_z) %>% add_markers(color = ~donor, colors = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#696969", "#010101"), alpha = 0.8, size = 50)

z <- NULL; dz <- NULL; plot <- NULL; tB <- NULL

## select samples within the mask 
f <- readnii(allTasks_mask)
#if value <0.05, turn value into 1 (create binary mask)
mask <- f<0.05 #turns into logical
img <- img_data(mask) #turns logical .nii into logical array
img <- img*1 # turns logical into numerical array
output = array(0, dim=c(77, 96, 79)) ##dim(f), same for all the fMRI files
for (x in c(1:dim(img)[1])) {
  for (y in c(1:dim(img)[2])) {
    for (z in c(1:dim(img)[3])) {
      if(img[x,y,z] == 1) {
        output[x,y,z] = 1
      }
    }
  }
}
img <- output; output <- NULL

#analyse img for 3Dinterpolation (img = effectmask as 3D array)
L <- f.read.nifti.header(allTasks_mask)
M <- t(as.matrix(dnr[,11:13])) ##convert to matrix and transpose x,y,z for xyz2ijk function 
VS <- as.data.frame(xyz2ijk(xyz = M, method = 2, L)) ##transform donordata to voxelspace of mask 
tVS <- as.data.frame(t(VS))## transpose for selecting column as vector
intp <- approx3d(x = c(1:77), y = c(1:96), z = c(1:79), f = img, xout = tVS$V1, yout = tVS$V2, zout = tVS$V3)##produces vector with all interpolated values per row of tVS in same order
plot(x = c(1:length(intp)), y = intp, type = 'l') ## quick plot results
names(intp) <- c(1:length(intp))
i <- as.data.frame(intp[intp >= 0.2]) ##threshold for including samples: at least 0.2 as interpolated mask-value
i <- as.numeric(rownames(i))

##all samples included in mask
iC <- as.numeric(); iC <- append(iC, i); iC <- sort(iC); iC <- unique(iC)

############################################################################################################################
## Step 4: compare relative gene expression levels in affected vs. unaffected cortical & subcortical areas
###########################################################################################################################

############# results across 6 donors ################

# samples included in mask 
dnrinc <- dnr[iC,]
dnrinc$dnrSample <- rownames(dnrinc); dnrinc <- dnrinc[colnames(dnrinc)[c(16, 1:15)]]
dnrinc$mask <- 'yes'
# samples excluded from mask 
dnrexc <- dnr[-iC,]
dnrexc$dnrSample <- rownames(dnrexc); dnrexc <- dnrexc[colnames(dnrexc)[c(16, 1:15)]] 

# create 4 dfs with info from all 6 donors 
maskcort <- inner_join(limmacort, dnrinc, by = c("well_id")) # joins by well_id, which is unique to each sample
masksubcort <- inner_join(limmasubcort, dnrinc, by = c("well_id"))
exccort <- inner_join(limmacort, dnrexc, by = c("well_id"))
excsubcort <- inner_join(limmasubcort, dnrexc, by = c("well_id"))

############ within-mask correlations ###############

## CORTICAL
cortcorrmask1 <- cor(maskcort[,1:6], method = "pearson") #r values in matrix
cortcorrmask <- signif(cortcorrmask, digits = 2) #2 decimals for plotting
cortcorrmask_p <- rcorr(as.matrix(maskcort[,1:6]), type = "pearson")$P #p values of calculated r values in matrix

## SUBCORTICAL 
subcortcorrmask1 <- cor(masksubcort[,1:6], method = "pearson")
subcortcorrmask <- signif(subcortcorrmask, digits = 2)
subcortcorrmask_p <- rcorr(as.matrix(masksubcort[,1:6]), type = "pearson")$P

############ outside mask correlations ###############

## CORTICAL
cortcorrexc1 <- cor(exccort[,1:6], method = "pearson")
cortcorrexc <- signif(cortcorrexc, digits = 2)
cortcorrexc_p <- rcorr(as.matrix(exccort[,1:6]), type = "pearson")$P

## SUBCORTICAL 
subcortcorrexc1 <- cor(excsubcort[,1:6], method = "pearson")
subcortcorrexc <- signif(subcortcorrexc, digits = 2)
subcortcorrexc_p <- rcorr(as.matrix(excsubcort[,1:6]), type = "pearson")$P


##Correct p-values and add "*" to r values if the correlation is significant after correction. 

#function to calculate the possible combinations of pairwise correlations
comb <- function(n, r) {
  #n = number of total variables
  #r = number of variables included in combination
  factorial(n) / factorial(n-r) / factorial(r)
}
combinations <- comb(6,2) #6 genes with pairwise (r=2) correlations

#function to add * if significant (<0.05)
significance <- function(x, y) {
  for (col in 1:ncol(x)) {
    for (row in 1:nrow(x)) {
      if(x[row,col]<0.05 & !is.na(x[row,col])){
        y[row,col] <- paste0(y[row,col], '*')
      }
    }
  }
  y
}

########### create matrix with significance ########
## within-mask correlations 

#cortical
cortcorrmask_pAdjusted <- apply(cortcorrmask_p, 2, p.adjust, method = "bonferroni", n = combinations) #Bonferroni correct p-values of Pearson's r
cortcorrmaskp <- significance(cortcorrmask_pAdjusted, cortcorrmask) #add * for each corrected p <0.05

#subcortical
subcortcorrmask_pAdjusted <- apply(subcortcorrmask_p, 2, p.adjust, method = "bonferroni", n = combinations) #Bonferroni correct p-values of Pearson's r
subcortcorrmaskp <- significance(subcortcorrmask_pAdjusted, subcortcorrmask) #add * for each corrected p <0.05

## outside-mask correlations 

#cortical
cortcorrexc_pAdjusted <- apply(cortcorrexc_p, 2, p.adjust, method = "bonferroni", n = combinations) #Bonferroni correct p-values of Pearson's r
cortcorrexcp <- significance(cortcorrexc_pAdjusted, cortcorrexc) #add * for each corrected p <0.05

#subcortical
subcortcorrexc_pAdjusted <- apply(subcortcorrexc_p, 2, p.adjust, method = "bonferroni", n = combinations) #Bonferroni correct p-values of Pearson's r
subcortcorrexcp <- significance(subcortcorrexc_pAdjusted, subcortcorrexc) #add * for each corrected p <0.05


#############################################################################################################################################################################
## STEP 5: test whether correlation matrices are equal for affected vs unaffected, and test for specific significant differences in correlation inside and outside fMRI mask.
#############################################################################################################################################################################

#cortical correlations (based on 1762 cortical samples from 6 donor brains: 138 samples in all-experiments mask; 1624 outside mask)
cortCortest <- cortest.normal(R1 = as.matrix(cortcorrmask1), R2 = as.matrix(cortcorrexc1), n1 = 138, n2 = 1624, fisher = TRUE)
cortCortest #no significant difference

#subcortical correlations (based on 986 subcortical samples from 6 donor brains: 208 samples in all-experiments mask; 778 outside mask)
subcortCortest <- cortest.normal(R1 = as.matrix(subcortcorrmask1), R2 = as.matrix(subcortcorrexc1), n1 = 208, n2 = 778, fisher = TRUE)
subcortCortest #significantly different

#####################################################################################################################
## subcortical OXTR correlations inside mask vs outside mask: no significant correlation in both situations for OXTR/AVPR2 and OXTR/AVPR1B
# subcortical correlation differences (inside vs outside mask) in OXTR and OXT/CD38/AVPR1a
subcortcorrmask_pAdjusted
subcortcorrexc_pAdjusted

#subcortical correlation difference of OXTR and OXT inside vs outside mask
subcorticalOXTR_OXT <- r.test(n = 208, r12 = 0.36750026, n2 = 778, r34 = 0.2020628646, twotailed = TRUE)
subcorticalOXTR_OXT$p
p.adjust(subcorticalOXTR_OXT$p, method = "bonferroni", n = 4)
0.36750026-0.2020628646

#subcortical correlation difference of OXTR and CD38 inside vs outside mask
subcorticalOXTR_CD38 <- r.test(n = 208, r12 = 0.42704158, n2 = 778, r34 = 0.4586527219, twotailed = TRUE)
subcorticalOXTR_CD38$p
p.adjust(subcorticalOXTR_CD38$p, method = "bonferroni", n = 4)
0.42704158-0.4586527219

#subcortical correlation difference of OXTR and AVPR1a inside vs outside mask
subcorticalOXTR_AVPR1a <- r.test(n = 208, r12 = 0.25346791, n2 = 778, r34 = -0.00832344, twotailed = TRUE)
subcorticalOXTR_AVPR1a$p
p.adjust(subcorticalOXTR_AVPR1a$p, method = "bonferroni", n = 4)
0.25346791 - -0.00832344

#subcortical correlation difference of OXTR and MC4R inside vs outside mask
subcortical_MC4R <- r.test(n = 208, r12 = 0.18545363, n2 = 778, r34 = 0.1541427282, twotailed = TRUE)
subcortical_MC4R$p
p.adjust(subcortical_MC4R$p, method = "bonferroni", n = 4)

# OXTR and AVPR1B, and OXTR and AVPR2 don't have significant correlations in either unaffected or affected subcortical samples,
# so Fisher's z-test for correlation differences is superfluous. These differences are not significant. 
0.033261924-0.04374125 #difference for OXTR and AVPR1B
-0.12495506--0.0003843006 ##difference for OXTR and AVPR2


########### heatmaps ##############

#subcortical affected samples
hmsubcortcorrmask <- Heatmap(subcortcorrmask, name = 'Subcortical\nCorrelations',
                             col = colorRamp2(c(-1,0,1), c("blue", "white", "red")),
                             cluster_rows = FALSE,
                             cluster_columns = FALSE,
                             row_names_gp = gpar(fontsize = 10),
                             row_names_side = c("left"),
                             row_title_rot = 0,
                             column_names_side = c("top"),
                             column_names_gp = gpar(fontsize = 10),
                             column_names_rot = 30,
                             width = unit(ncol(subcortcorrmask)*3, "lines"), 
                             height = unit(nrow(subcortcorrmask)*3, "lines"),
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.text(sprintf("%s", subcortcorrmaskp[i, j]), x, y, gp = gpar(fontsize = 8)) }
                             
)

draw(hmsubcortcorrmask)

#subcortical unaffected samples
hmsubcortcorrexc <- Heatmap(subcortcorrexc, name = 'Subcortical\nCorrelations',
                         cluster_rows = FALSE,
                         cluster_columns = FALSE,
                         row_names_gp = gpar(fontsize = 10),
                         row_names_side = c("left"),
                         row_title_rot = 0,
                         column_names_side = c("top"),
                         column_names_gp = gpar(fontsize = 10),
                         column_names_rot = 30,
                         width = unit(ncol(subcortcorrexc)*3, "lines"), 
                         height = unit(nrow(subcortcorrexc)*3, "lines"),
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           grid.text(sprintf("%s", subcortcorrexcp[i, j]), x, y, gp = gpar(fontsize = 8)) 
                          })
draw(hmsubcortcorrexc)

#cortical affected samples
hmcortcorrmask <- Heatmap(cortcorrmask, name = 'Cortical\nCorrelations',
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          row_names_gp = gpar(fontsize = 10),
                          row_names_side = c("left"),
                          row_title_rot = 0,
                          column_names_side = c("top"),
                          column_names_gp = gpar(fontsize = 10),
                          column_names_rot = 30,
                          width = unit(ncol(cortcorrmask)*3, "lines"), 
                          height = unit(nrow(cortcorrmask)*3, "lines"),
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%s", cortcorrmaskp[i, j]), x, y, gp = gpar(fontsize = 8)) }
)

draw(hmcortcorrmask)

#cortical unaffected samples
hmcortcorrexc <- Heatmap(cortcorrexc, name = 'Cortical\nCorrelations',
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      row_names_gp = gpar(fontsize = 10),
                      row_names_side = c("left"),
                      row_title_rot = 0,
                      column_names_side = c("top"),
                      column_names_gp = gpar(fontsize = 10),
                      column_names_rot = 30,
                      width = unit(ncol(cortcorrexc)*3, "lines"), 
                      height = unit(nrow(cortcorrexc)*3, "lines"),
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%s", cortcorrexcp[i, j]), x, y, gp = gpar(fontsize = 8)) }
)

draw(hmcortcorrexc)


                         