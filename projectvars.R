# libraries
library(ggplot2);
library(RColorBrewer);
library(paletteer);
library(scales); # function show_col
library(tidyverse);
library(stringr);
library(rlang); # required for the bang-bang operator "!!" used for variable repalcement in ggplot
library(ComplexHeatmap);
library(circlize);
library(igraph);                                                                                    
library(ggraph);                                                                                    
library(graphlayouts);                                                                              
library(ggforce);                                                                                   
library(ggpubr);                                                                                    
library(extrafont);

sessionInfo();

options(width=120,prompt="HEATMAPS> ");

WD <- Sys.getenv("WD");
ODRed <- Sys.getenv("ODIRed");
ODRaj <- Sys.getenv("ODIRaj");

DS <- c( # vector of dataset IDs
    "rbf_only",
    "rbf_x_allneural",
    "rbf_x_neural-aminergic",
    "rbf_x_neural-colinergic",
    "rbf_x_neural-gabaergic",
    "rbf_x_neural-glutamergic",
    "rbf_x_neural-glycinergic",
    "rbf_x_neural-panneuronal",
    "rbf_x_neural-peptidergic",
    "rbf_x_parenchimal+glial",
    "rbf_x_parenchimal",
    "rbf_x_glial",
    "rbf_x_ligth",
    "rbf_x_allsets"
    );

DG <- read.table(paste0(WD,"/heatmaps_gene_matrix.tbl"), sep="\t", header=FALSE);
colnames(DG) <- c(DS,"geneid","symbol","name","set");
rownames(DG) <- DG$geneid;
dim(DG);
# [1] 41 18

head(DG,2);
#   rbf_only rbf_x_allneural rbf_x_neural-aminergic rbf_x_neural-colinergic
# 1        1               1                      1                       1
# 2        1               1                      1                       1
#   rbf_x_neural-gabaergic rbf_x_neural-glutamergic rbf_x_neural-glycinergic
# 1                      1                        1                        1
# 2                      1                        1                        1
#   rbf_x_neural-panneuronal rbf_x_neural-peptidergic rbf_x_parenchimal+glial
# 1                        1                        1                       1
# 2                        1                        1                       1
#   rbf_x_parenchimal rbf_x_glial rbf_x_ligth rbf_x_allsets               geneid
# 1                 1           1           1             1 dd_Smed_v6_11730_0_1
# 2                 1           1           1             1  dd_Smed_v6_5546_0_1
#   symbol                                       name        set
# 1   Rft1        riboflavin transporter [phagocytes] RF-targets
# 2   Rft2 riboflavin transporter [glial/parenchymal] RF-targets

pfxRed <- paste0(ODRed,"/merged_coexpression_riboflavin");
pfxRaj <- paste0(ODRaj,"/merged_coexpression_riboflavin");
sfxRed <- "_Reddien.tbl";
sfxRaj <- "_Rajewsky.tbl";

DT.Red <- list(
    "rbf_only"                 = read.table(paste0(pfxRed,"_only",                sfxRed), sep="\t", header=TRUE),
    "rbf_x_allneural"          = read.table(paste0(pfxRed,"_x_allneural",         sfxRed), sep="\t", header=TRUE),
    "rbf_x_neural-aminergic"   = read.table(paste0(pfxRed,"_x_neural-aminergic",  sfxRed), sep="\t", header=TRUE),
    "rbf_x_neural-colinergic"  = read.table(paste0(pfxRed,"_x_neural-colinergic", sfxRed), sep="\t", header=TRUE),
    "rbf_x_neural-gabaergic"   = read.table(paste0(pfxRed,"_x_neural-gabaergic",  sfxRed), sep="\t", header=TRUE),
    "rbf_x_neural-glutamergic" = read.table(paste0(pfxRed,"_x_neural-glutamergic",sfxRed), sep="\t", header=TRUE),
    "rbf_x_neural-glycinergic" = read.table(paste0(pfxRed,"_x_neural-glycinergic",sfxRed), sep="\t", header=TRUE),
    "rbf_x_neural-panneuronal" = read.table(paste0(pfxRed,"_x_neural-panneuronal",sfxRed), sep="\t", header=TRUE),
    "rbf_x_neural-peptidergic" = read.table(paste0(pfxRed,"_x_neural-peptidergic",sfxRed), sep="\t", header=TRUE),
    "rbf_x_parenchimal+glial"  = read.table(paste0(pfxRed,"_x_parenchimal+glial", sfxRed), sep="\t", header=TRUE),
    "rbf_x_parenchimal"        = read.table(paste0(pfxRed,"_x_parenchimal",       sfxRed), sep="\t", header=TRUE),
    "rbf_x_glial"              = read.table(paste0(pfxRed,"_x_glial",             sfxRed), sep="\t", header=TRUE),
    "rbf_x_ligth"              = read.table(paste0(pfxRed,"_x_ligth",             sfxRed), sep="\t", header=TRUE),
    "rbf_x_allsets"            = read.table(paste0(pfxRed,"_x_allsets",           sfxRed), sep="\t", header=TRUE)
    );
  
DT.Raj <- list(
    "rbf_only"                 = read.table(paste0(pfxRaj,"_only",                sfxRaj), sep="\t", header=TRUE),
    "rbf_x_allneural"          = read.table(paste0(pfxRaj,"_x_allneural",         sfxRaj), sep="\t", header=TRUE),
    "rbf_x_neural-aminergic"   = read.table(paste0(pfxRaj,"_x_neural-aminergic",  sfxRaj), sep="\t", header=TRUE),
    "rbf_x_neural-colinergic"  = read.table(paste0(pfxRaj,"_x_neural-colinergic", sfxRaj), sep="\t", header=TRUE),
    "rbf_x_neural-gabaergic"   = read.table(paste0(pfxRaj,"_x_neural-gabaergic",  sfxRaj), sep="\t", header=TRUE),
    "rbf_x_neural-glutamergic" = read.table(paste0(pfxRaj,"_x_neural-glutamergic",sfxRaj), sep="\t", header=TRUE),
    "rbf_x_neural-glycinergic" = read.table(paste0(pfxRaj,"_x_neural-glycinergic",sfxRaj), sep="\t", header=TRUE),
    "rbf_x_neural-panneuronal" = read.table(paste0(pfxRaj,"_x_neural-panneuronal",sfxRaj), sep="\t", header=TRUE),
    "rbf_x_neural-peptidergic" = read.table(paste0(pfxRaj,"_x_neural-peptidergic",sfxRaj), sep="\t", header=TRUE),
    "rbf_x_parenchimal+glial"  = read.table(paste0(pfxRaj,"_x_parenchimal+glial", sfxRaj), sep="\t", header=TRUE),
    "rbf_x_parenchimal"        = read.table(paste0(pfxRaj,"_x_parenchimal",       sfxRaj), sep="\t", header=TRUE),
    "rbf_x_glial"              = read.table(paste0(pfxRaj,"_x_glial",             sfxRaj), sep="\t", header=TRUE),
    "rbf_x_ligth"              = read.table(paste0(pfxRaj,"_x_ligth",             sfxRaj), sep="\t", header=TRUE),
    "rbf_x_allsets"            = read.table(paste0(pfxRaj,"_x_allsets",           sfxRaj), sep="\t", header=TRUE)
    );
    
dim(DT.Red$rbf_only);
# [1] 11352  13 -> OK
dim(DT.Raj$rbf_only);
# [1]   996  13 -> OK


HTMP.order <-  c(
     "Head",     #
     "PrePhar",  #
     "Pharynx",  #
     "Trunk",    #
     "Tail",     #
     "Whole"     #
     );
