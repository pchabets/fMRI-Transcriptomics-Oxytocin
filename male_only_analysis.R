#########################################################################################################################################################
## Author: C Mclain, P.C. Habets ## Department: Internal Medicine, Division of Endocrinology  ####
############################################################################################################################################################
## Institution: Leiden Medical Center Utrecht ## Script: differential gene expression analysis using only male donor brains and the male-only fMRI file#####
############################################################################################################################################################

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
library(limma)
library(GO.db)

setwd("/path/to/local/folder/data") #set working directory to local folder containing the AHBA downloaded files. See readme file for setup of local folder with separate 'data' directory that holds all needed files.
cores <- detectCores()

######################################################################
## Step 1: find the highest-intensity probe for each gene of interest
######################################################################

# find all (reannotated) probes for each gene of interest

# oxytocin pathway genes 
OXT <- c("1053421", "1053420","1053419")
OXTR <- c("1053418", "1053417")
CD38 <- c("1058445","1058446", "1061956")

# vasopressin receptor genes
AVPR2 <- c("1012380", "1011597", "1018065", "1013120")
AVPR1A <- c("1058917", "1058918", "1058919", "1066666")
AVPR1B <- c("1058914", "1058915", "1058916")

probes <- c(OXT, OXTR, CD38, AVPR2, AVPR1A, AVPR1B)

# find the expression of each probe across all samples 
exprList <- mclapply(1:6, mc.cores = cores, function(d) {
  FileMicroarray <- paste0("normalized_microarray_donor0", d, "/MicroarrayExpression.csv")
  e <- read.csv(FileMicroarray, header = FALSE)
  rownames(e) <- e[,1]
  e <- e[,-1] # get rid of ProbeIDs
  
  FileSample <- paste0("normalized_microarray_donor0", d, "/SampleAnnot.csv")
  sample_annot <- read.csv(FileSample, header = TRUE)
  
  # remove cerebellum and brainstem samples
  indBS <- which(sample_annot$slab_type == "BS")
  indCB <- which(sample_annot$slab_type =="CB")
  BSCB <- c(indBS, indCB)
  
  # take the expression just for probes of interest
  e <- e[probes,-BSCB]
  e
})
expression <- Reduce(cbind, exprList)

# take the mean of each 
exprmeans <- data.frame(rowMeans(expression))
View (exprmeans)

# take the probe for each gene with the highest expression 

probes <- c("1053419", "1053417", "1058445", "1013120", "1058917", "1058916")
# same results as with highest RNA seq correlation 

######################################################################################################
## Step 2: remove batch effects across all 6 brains separately for cortical and subcortical samples
######################################################################################################

cort <- read.csv("corticalSamplesonly.csv", header = TRUE) # all cortical samples
cortical <- as.vector(unique(cort[,5])) # cortical acronyms

donors <- c(1, 2, 3, 4, 6) # remove the 5 for male-only analysis 

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
  if (d!=6) {
    ind <- cort_ind[[d]] }
  else {
    ind <- cort_ind[[5]] }
  e <- e[probes,ind] # only columns corresponding to cortical samples of genes of interest
  e[7,] <- d
  rownames(e) <- c("OXT", "OXTR", "CD38", "AVPR2", "AVPR1A", "AVPR1B", "donor")
  e 
})
cort_expression <- Reduce(cbind, cort_expr)
limmacort <- removeBatchEffect(cort_expression, batch = cort_expression[7,])

# subcortical expression for each gene of interest 
subcort_expr <- mclapply(donors, mc.cores=cores, function(d) {
  FileMicroarray <- paste0("normalized_microarray_donor0", d, "/MicroarrayExpression.csv")
  e <- read.csv(FileMicroarray, header = FALSE)
  rownames(e) <- e[,1]
  e <- e[,-1] # get rid of probe IDs
  # remove cortical, brain stem and cerebellum samples
  if (d!=6) {
    indcort <- cort_ind[[d]] }
  else {
    indcort <- cort_ind[[5]] }
  sample_annot <- read.csv(paste0("normalized_microarray_donor0", d, "/SampleAnnot.csv"), header = TRUE)
  indBSCB <- c(which(sample_annot$slab_type == "BS"), which(sample_annot$slab_type =="CB"))
  indthrow <- c(indcort, indBSCB)
  e <- e[probes, -indthrow]
  e[7,] <- d
  rownames(e) <- c("OXT", "OXTR", "CD38", "AVPR2", "AVPR1A", "AVPR1B", "donor")
  e
})

subcort_expression <- Reduce(cbind, subcort_expr)
limmasubcort <- removeBatchEffect(subcort_expression, batch = subcort_expression[7,])

###########################################################
## Step 3:  z-normalize for each gene,
## separately for subcortical and cortical samples
###########################################################

## z-normalize cortical expression 
expr <- limmacort
meanoxt <- mean(expr[1,]); sdoxt <- sd(unlist(expr[1,]))
meanoxtr <- mean(expr[2,]); sdoxtr <- sd(unlist(expr[2,]))
meancd <- mean(expr[3,]); sdcd <- sd(unlist(expr[3,]))
meanavpr2 <- mean(expr[4,]); sdavpr2 <- sd(unlist(expr[4,]))
meanavpr1a <- mean(expr[5,]); sdavpr1a <- sd(unlist(expr[5,]))
meanavpr1b <- mean(expr[6,]); sdavpr1b <- sd(unlist(expr[6,]))

zoxt <- (expr[1,]- meanoxt)/sdoxt 
zoxtr <- (expr[2,]- meanoxtr)/sdoxtr
zcd <- (expr[3,]- meancd)/sdcd 
zavpr2 <- (expr[4,]- meanavpr2)/sdavpr2
zavpr1a <- (expr[5,]- meanavpr1a)/sdavpr1a
zavpr1b <- (expr[6,]- meanavpr1b)/sdavpr1b

cort_annot <- lapply(donors, function(d) { sample_annot <- read.csv(paste0("normalized_microarray_donor0", d, "/SampleAnnot.csv"), header = TRUE)
if (d!=6) {
sample_annot <- sample_annot[cort_ind[[d]],] }
else {
sample_annot <- sample_annot[cort_ind[[5]],]
}
sample_annot$well_id })

well_id <- c(unlist(cort_annot))  # to enable inner_join later

cortzscores <- data.frame(t(rbind(zoxt, zoxtr, zcd,zavpr2, zavpr1a, zavpr1b, well_id)))
colnames(cortzscores) <- c("OXT", "OXTR", "CD38", "AVPR2", "AVPR1A", "AVPR1B", "well_id")

## z-normalize subcortical expression 
expr <- limmasubcort
meanoxt <- mean(expr[1,]); sdoxt <- sd(unlist(expr[1,]))
meanoxtr <- mean(expr[2,]); sdoxtr <- sd(unlist(expr[2,]))
meancd <- mean(expr[3,]); sdcd <- sd(unlist(expr[3,]))
meanavpr2 <- mean(expr[4,]); sdavpr2 <- sd(unlist(expr[4,]))
meanavpr1a <- mean(expr[5,]); sdavpr1a <- sd(unlist(expr[5,]))
meanavpr1b <- mean(expr[6,]); sdavpr1b <- sd(unlist(expr[6,]))

zoxt <- (expr[1,]- meanoxt)/sdoxt 
zoxtr <- (expr[2,]- meanoxtr)/sdoxtr
zcd <- (expr[3,]- meancd)/sdcd 
zavpr2 <- (expr[4,]- meanavpr2)/sdavpr2
zavpr1a <- (expr[5,]- meanavpr1a)/sdavpr1a
zavpr1b <- (expr[6,]- meanavpr1b)/sdavpr1b


subcort_annot <- lapply(donors, function(d) { sample_annot <- read.csv(paste0("normalized_microarray_donor0", d, "/SampleAnnot.csv"), header = TRUE)
# remove cortical, brain stem and cerebellum samples
if (d!=6) {
indcort <- cort_ind[[d]] }
else {
  indcort <- cort_ind[[5]] 
}
sample_annot <- read.csv(paste0("normalized_microarray_donor0", d, "/SampleAnnot.csv"), header = TRUE)
indBSCB <- c(which(sample_annot$slab_type == "BS"), which(sample_annot$slab_type =="CB"))
indthrow <- c(indcort, indBSCB)
sample_annot <- sample_annot[-indthrow,]
sample_annot$well_id
}) 

well_id <- c(unlist(subcort_annot))

subcortzscores <- data.frame(t(rbind(zoxt, zoxtr, zcd,zavpr2, zavpr1a, zavpr1b, well_id))) # wellID row to enable inner_join later
colnames(subcortzscores) <- c("OXT", "OXTR", "CD38", "AVPR2", "AVPR1A", "AVPR1B", "well_id")

#####################################################################
## Step 4: assign samples to affected/unaffected areas on the fMRI map
#####################################################################

##AHBA donor brains
donorNames <- c("donor01", "donor02", "donor03", "donor04", "donor06")
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

# read in file
z <- readnii("MALE_OXT-PBO_P.nii") 

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
f <- readnii("MALE_OXT-PBO_P.nii")
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
L <- f.read.nifti.header("MALE_OXT-PBO_P.nii")
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
## Step 5: compare relative gene expression levels in affected vs. unaffected cortical & subcortical areas
###########################################################################################################################

############################# results across 6 donors #################################

# samples included in mask 
dnrinc <- dnr[iC,]
dnrinc$dnrSample <- rownames(dnrinc); dnrinc <- dnrinc[colnames(dnrinc)[c(16, 1:15)]]
dnrinc$mask <- 'yes'
# samples excluded from mask 
dnrexc <- dnr[-iC,]
dnrexc$dnrSample <- rownames(dnrexc); dnrexc <- dnrexc[colnames(dnrexc)[c(16, 1:15)]] 

# create 4 dfs with info from all 5 donors 
maskcort <- inner_join(cortzscores, dnrinc, by = c("well_id")) # joins by well_id, which is unique to each sample
masksubcort <- inner_join(subcortzscores, dnrinc, by = c("well_id"))
exccort <- inner_join(cortzscores, dnrexc, by = c("well_id"))
excsubcort <- inner_join(subcortzscores, dnrexc, by = c("well_id"))

# mean expression for each gene by cortical/subcortical
cortical_affected <- c(mean(maskcort$OXT), mean(maskcort$OXTR), mean(maskcort$CD38), mean(maskcort$AVPR2), mean(maskcort$AVPR1A), mean(maskcort$AVPR1B))
subcortical_affected <- c(mean(masksubcort$OXT), mean(masksubcort$OXTR), mean(masksubcort$CD38), mean(masksubcort$AVPR2), mean(masksubcort$AVPR1A), mean(masksubcort$AVPR1B))
p_values_cort <- c(wilcox.test(maskcort$OXT,exccort$OXT)$p.value, wilcox.test(maskcort$OXTR,exccort$OXTR)$p.value, wilcox.test(maskcort$CD38,exccort$CD38)$p.value, 
                   wilcox.test(maskcort$AVPR2,exccort$AVPR2)$p.value, wilcox.test(maskcort$AVPR1A, exccort$AVPR1A)$p.value, wilcox.test(maskcort$AVPR1B, exccort$AVPR1B)$p.value) 

p_values_cort_adj <- p.adjust(p_values_cort, method = "bonferroni")

cortical_unaffected <- c(mean(exccort$OXT), mean(exccort$OXTR), mean(exccort$CD38), mean(exccort$AVPR2), mean(exccort$AVPR1A), mean(exccort$AVPR1B))
subcortical_unaffected <- c(mean(excsubcort$OXT), mean(excsubcort$OXTR), mean(excsubcort$CD38), mean(excsubcort$AVPR2), mean(excsubcort$AVPR1A), mean(excsubcort$AVPR1B))
p_values_subcort <- c(wilcox.test(masksubcort$OXT,excsubcort$OXT)$p.value, wilcox.test(masksubcort$OXTR, excsubcort$OXTR)$p.value, wilcox.test(masksubcort$CD38, excsubcort$CD38)$p.value,
                      wilcox.test(masksubcort$AVPR2,excsubcort$AVPR2)$p.value, wilcox.test(masksubcort$AVPR1A, excsubcort$AVPR1A)$p.value, wilcox.test(masksubcort$AVPR1B, excsubcort$AVPR1B)$p.value) 

p_values_subcort_adj <- p.adjust(p_values_subcort, method = "bonferroni")

meanzscores <- data.frame(cortical_affected, cortical_unaffected, p_values_cort, p_values_cort_adj, subcortical_affected, subcortical_unaffected, p_values_subcort, p_values_subcort_adj)
rownames(meanzscores) <- c("OXT", "OXTR", "CD38", "AVPR2","AVPR1A", "AVPR1B")
View(meanzscores)

## checking for number of right-hemisphere samples

right = sum(maskcort$mni_x>0, masksubcort$mni_x>0)
left =  sum(maskcort$mni_x<0, masksubcort$mni_x<0)


## means within affected structures for oxytocin genes
cort_aff_oxy <- aggregate(maskcort[,1:3], list(maskcort$structure_acronym), mean) 
cort_aff_oxy <- cort_aff_oxy[order(cort_aff_oxy$OXTR),] # ordering by OXTR 
cort_aff_oxy$affected <- as.factor("yes")

subcort_aff_oxy <- aggregate(masksubcort[,1:3], list(masksubcort$structure_acronym), mean)
subcort_aff_oxy <- subcort_aff_oxy[order(subcort_aff_oxy$OXTR),] 
subcort_aff_oxy$affected <- as.factor("yes")

## means within unaffected structures for oxytocin genes
cort_unaff_oxy <- aggregate(exccort[,1:3], list(exccort$structure_acronym), mean) 
cort_unaff_oxy <- cort_unaff_oxy[order(cort_unaff_oxy$OXTR),] 
cort_unaff_oxy$affected <- as.factor("no")

subcort_unaff_oxy <- aggregate(excsubcort[,1:3], list(excsubcort$structure_acronym), mean)
subcort_unaff_oxy <- subcort_unaff_oxy[order(subcort_unaff_oxy$OXTR),] 
subcort_unaff_oxy$affected <- as.factor("no")

## means within affected structures for vasopressin genes

cort_aff_vaso <- aggregate(maskcort[,4:6], list(maskcort$structure_acronym), mean) 
cort_aff_vaso <- cort_aff_vaso[order(cort_aff_vaso$AVPR2),] # ordering by AVPR2 because it shows the largest difference
cort_aff_vaso$affected <- as.factor("yes")

subcort_aff_vaso <- aggregate(masksubcort[,4:6], list(masksubcort$structure_acronym), mean)
subcort_aff_vaso <- subcort_aff_vaso[order(subcort_aff_vaso$AVPR2),] 
subcort_aff_vaso$affected <- as.factor("yes")

## means within unaffected structures for oxytocin genes
cort_unaff_vaso <- aggregate(exccort[,4:6], list(exccort$structure_acronym), mean) 
cort_unaff_vaso <- cort_unaff_vaso[order(cort_unaff_vaso$AVPR2),] 
cort_unaff_vaso$affected <- as.factor("no")

subcort_unaff_vaso <- aggregate(excsubcort[,4:6], list(excsubcort$structure_acronym), mean)
subcort_unaff_vaso <- subcort_unaff_vaso[order(subcort_unaff_vaso$AVPR2),] 
subcort_unaff_vaso$affected <- as.factor("no")


####################### heatmaps #########################
## run these independently based on the desired outcome ## 

###### affected samples only ######

## oxytocin genes

cort_aff_oxy <- cort_aff_oxy[,c(3,2,4,1)]
cortmatoxy <- data.matrix(cort_aff_oxy)
rownames(cortmatoxy) <- cort_aff_oxy$Group.1

hmcortoxy <- Heatmap(cortmatoxy[,1:3],name = 'z-score',
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     row_names_gp = gpar(fontsize = 6),
                     row_names_side = c("left"),
                     row_title_rot = 0,
                     column_names_side = c("top"),
                     column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                     column_names_rot = 45,
                     width = unit(ncol(cort_aff_oxy)*1.5, "lines"), 
                     height = unit(nrow(cort_aff_oxy)*0.5, "lines")
)

draw(hmcortoxy)

subcort_aff_oxy <- subcort_aff_oxy[,c(3,2,4,1)]
subcortmatoxy <- data.matrix(subcort_aff_oxy)
rownames(subcortmatoxy) <- subcort_aff_oxy$Group.1

hmsubcortoxy <- Heatmap(subcortmatoxy[,1:3],name = 'zscore',
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        row_names_gp = gpar(fontsize = 6),
                        row_names_side = c("left"),
                        row_title_rot = 0,
                        column_names_side = c("top"),
                        column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                        column_names_rot = 45,
                        width = unit(ncol(cort_aff_oxy)*1.5, "lines"), 
                        height = unit(nrow(cort_aff_oxy)*0.5, "lines")
)

draw(hmsubcortoxy)


####### affected and unaffected samples #######
#in order for this to work, run lines 322 - 357 again

## oxytocin genes
cort_aff_oxy <- cort_aff_oxy[,c(3,2,4,5,1)]
cort_unaff_oxy <- cort_unaff_oxy[,c(3,2,4,5,1)]
cortmatoxy <- data.matrix(rbind(cort_aff_oxy[,1:4], cort_unaff_oxy[,1:4]))
rownames(cortmatoxy) <- c(cort_aff_oxy$Group.1, cort_unaff_oxy$Group.1)

hmcortoxy <- Heatmap(cortmatoxy[,1:3], name = 'Cortical\nZ-Score\nexpression',
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     row_names_gp = gpar(fontsize = 3),
                     row_names_side = c("left"),
                     row_title_rot = 0,
                     column_names_side = c("top"),
                     column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                     column_names_rot = 45,
                     width = unit(ncol(cort_aff_oxy)*1.5, "lines"), 
                     height = unit(nrow(cort_aff_oxy)*0.7, "lines"),
                     row_split = cortmatoxy[,4],
                     row_title = c("1" = "Affected", "2" = "Unaffected")
)

draw(hmcortoxy)


subcort_aff_oxy <- subcort_aff_oxy[,c(3,2,4,5,1)]
subcort_unaff_oxy <- subcort_unaff_oxy[,c(3,2,4,5,1)]
subcortmatoxy <- data.matrix(rbind(subcort_aff_oxy[,1:4], subcort_unaff_oxy[,1:4]))
rownames(subcortmatoxy) <- c(subcort_aff_oxy$Group.1, subcort_unaff_oxy$Group.1)

hmsubcortoxy <- Heatmap(subcortmatoxy[,1:3], name = 'Subcortical\nZ-Score\nexpression',
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        row_names_gp = gpar(fontsize = 2.5),
                        row_names_side = c("left"),
                        row_title_rot = 0,
                        column_names_side = c("top"),
                        column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                        column_names_rot = 45,
                        width = unit(ncol(subcort_aff_oxy)*1.5, "lines"), 
                        height = unit(nrow(subcort_aff_oxy)*0.45, "lines"),
                        row_split = subcortmatoxy[,4],
                        row_title = c("1" = "Affected", "2" = "Unaffected")
)

draw(hmsubcortoxy)

## vasopressin genes 
cortmatvaso <- data.matrix(rbind(cort_aff_vaso[,2:5], cort_unaff_vaso[,2:5]))
rownames(cortmatvaso) <- c(cort_aff_vaso$Group.1, cort_unaff_vaso$Group.1)

hmcortvaso <- Heatmap(cortmatvaso[,1:3], name = 'Cortical\nZ-Score\nexpression',
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      row_names_gp = gpar(fontsize = 3),
                      row_names_side = c("left"),
                      row_title_rot = 0,
                      column_names_side = c("top"),
                      column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                      column_names_rot = 45,
                      width = unit(ncol(cort_aff_vaso)*1.5, "lines"), 
                      height = unit(nrow(cort_aff_vaso)*0.7, "lines"),
                      row_split = cortmatvaso[,4],
                      row_title = c("1" = "Affected", "2" = "Unaffected")
)

draw(hmcortvaso)

subcortmatvaso <- data.matrix(rbind(subcort_aff_vaso[,2:5], subcort_unaff_vaso[,2:5]))
rownames(subcortmatvaso) <- c(subcort_aff_vaso$Group.1, subcort_unaff_vaso$Group.1)

hmsubcortvaso <- Heatmap(subcortmatvaso[,1:3], name = 'Subcortical\nZ-Score\nexpression',
                         cluster_rows = FALSE,
                         cluster_columns = FALSE,
                         row_names_gp = gpar(fontsize = 2.5),
                         row_names_side = c("left"),
                         row_title_rot = 0,
                         column_names_side = c("top"),
                         column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                         column_names_rot = 45,
                         width = unit(ncol(subcort_aff_vaso)*1.5, "lines"), 
                         height = unit(nrow(subcort_aff_vaso)*0.4, "lines"),
                         row_split = subcortmatvaso[,4],
                         row_title = c("1" = "Affected", "2" = "Unaffected")
)

draw(hmsubcortvaso)

######## combined oxytocin and vasopressin heatmaps ############
######## affected and unaffected #########
#in order for this to work, run lines 322 - 357 again

cort_aff_all <- merge(cort_aff_oxy, cort_aff_vaso, by = "Group.1")
cort_aff_all <- cort_aff_all[order(cort_aff_all$OXTR),] # ordering by OXTR 
cort_unaff_all <- merge(cort_unaff_oxy, cort_unaff_vaso, by = "Group.1")
cort_unaff_all <- cort_unaff_all[order(cort_unaff_all$OXTR),] 
cortmatall <- data.matrix(rbind(cort_aff_all[,c(3,2,4,6:9)], cort_unaff_all[,c(3,2,4,6:9)]))
rownames(cortmatall) <- c(cort_aff_all$Group.1, cort_unaff_all$Group.1)

hmcortall <- Heatmap(cortmatall[,1:6], name = 'z-score',
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     row_names_gp = gpar(fontsize = 4),
                     row_names_side = c("left"),
                     row_title_rot = 0,
                     column_names_side = c("top"),
                     column_names_gp = gpar(fontsize = 7, fontface = "italic"),
                     column_names_rot = 45,
                     width = unit(ncol(cort_aff_oxy)*1.5, "lines"), 
                     height = unit(nrow(cort_aff_oxy)*0.7, "lines"),
                     row_split = cortmatall[,7],
                     row_title = c("1" = "Affected", "2" = "Unaffected")
)

draw(hmcortall)


subcort_aff_all <- merge(subcort_aff_oxy, subcort_aff_vaso, by = "Group.1")
subcort_aff_all <- subcort_aff_all[order(subcort_aff_all$OXTR),] # ordering by OXTR 
subcort_unaff_all <- merge(subcort_unaff_oxy, subcort_unaff_vaso, by = "Group.1")
subcort_unaff_all <- subcort_unaff_all[order(subcort_unaff_all$OXTR),] 
subcortmatall <- data.matrix(rbind(subcort_aff_all[,c(3,2,4,6:9)], subcort_unaff_all[,c(3,2,4,6:9)]))
rownames(subcortmatall) <- c(subcort_aff_all$Group.1, subcort_unaff_all$Group.1)

hmsubcortall <- Heatmap(subcortmatall[,1:6], name = 'z-score',
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        row_names_gp = gpar(fontsize = 4),
                        row_names_side = c("left"),
                        row_title_rot = 0,
                        column_names_side = c("top"),
                        column_names_gp = gpar(fontsize = 7, fontface = "italic"),
                        column_names_rot = 45,
                        width = unit(ncol(subcort_aff_all)*1, "lines"), 
                        height = unit(nrow(subcort_aff_all)*0.58, "lines"),
                        row_split = subcortmatall[,7],
                        row_title = c("1" = "Affected", "2" = "Unaffected")
)

draw(hmsubcortall)

