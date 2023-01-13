## =============================================================================
## Description
#
# This script contains the code for our simulation study.
# There are in total 6 cases considered:
# 4nodes-sparse, 4nodes-dense, 5nodes-sparse,
# 5nodes-dense, 6nodes-sparse, 6nodes-dense.
# We estimated a GGM and PAG (using CCD algorithm) for each case.
# We compute overall density and degree per node for each model for comparison.
# In step 7 and 8, we show the plots comparing density
# and degree centrality across all different conditions.
## =============================================================================


## =============================================================================
## 0) Preparation
#
# First, load all necessary packages and functions.
## =============================================================================
library(qgraph)
library(pcalg)
library(qgraph)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(purrr)
library(magrittr)
library(tidyr)
library(stringr)

source("R/CCD_fnc.R")
source("R/plot_fnc.R")
source("R/dsep_fnc.R")
source("R/searchAM_KP_fnc.R")
source("R/equivset_fnc.R")
source("R/data_generating_fnc.R")
source("R/eval_metric_fnc.R")


## =============================================================================
## 1) Four nodes - sparse model (low density)
## =============================================================================

## Specify B matrix
# set the number of nodes (p)
p = 4
B4 = matrix(c(0, 0, 0, 0,
              1, 0, 0.5, 0,
              0, 0.5, 0, 0.9,
              0, 0, 0, 0), p, p, byrow = T)
colnames(B4) <- c("X1", "X2", "X3", "X4")


## Generate data
# first, equilibrium check (necessary condition for cyclic models to converge)
equilibrium_check(B4)
# generated data with N = 10^6, seed = 1
data4p <- gen_dat(B4, N =1e6, seed = 1)

## Specify layout
layout4 = matrix(c(-1,1,
                   -1,0,
                   1,0,
                   1,1),4,2,byrow = T)

## True cyclic graph
true4p <- qgraph(t(B4), layout=layout4, labels = colnames(B4), theme="colorblind")

## Estimate GGM
ggm4p <- qgraph(cor(data4p), layout=layout4, theme="colorblind")

## Run CCD algorithm
ccd_4p <- ccdKP(df=data4p, dataType = "continuous", alpha = 0.05)
mat4p <- CreateAdjMat(ccd_4p, 4)

## Estimate PAG
pag4p <- plotPAG(ccd_4p, mat4p)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv4p <- semiequiv_dcg(ccd_4p, mat4p)
# save(equiv4p, file="data/equiv4p.RData")
load("data/equiv4p.RData")

## Density comparison
# true model density
truemoddensity(B4)
# GGM density
GGMdensity(ggm4p)
# average density for all equivalent DCGs
DCGdensity(equiv4p)

## Degree comparison
# true model degrees per node
GGMdegree(true4p)
# GGM  degrees per node
GGMdegree(ggm4p)
# average degree per node for all equivalent DCGs
DCGdegree(equiv4p)

## =============================================================================
## 2) Four nodes - dense model (high density)
## =============================================================================

## Specify B matrix
# set the number of nodes (p)
p = 4
B4_high = matrix(c(0, 0, 0, 0,
                   0.9, 0, 0.4, 0,
                   0, 0.5, 0, .5,
                   -0.8, 0, 0, 0), p, p, byrow = T)
colnames(B4_high) <- c("X1", "X2", "X3", "X4")


## Generate data
# first, equilibrium check (necessary condition for cyclic models to converge)
equilibrium_check(B4_high)
# generated data with N = 10^6, seed = 1
data4p_high <- gen_dat(B4_high, N =1e6, seed = 1)

## Specify layout
layout4 = matrix(c(-1,1,
                   -1,0,
                   1,0,
                   1,1),4,2,byrow = T)

## True cyclic graph
true4p_high <- qgraph(t(B4_high), layout=layout4, labels = colnames(B4_high), theme="colorblind")

## Estimate GGM
ggm4p_high <- qgraph(t(cor(data4p_high)), layout=layout4, theme="colorblind")

## Run CCD algorithm
ccd_4p_high <- ccdKP(df=data4p_high, dataType = "continuous", alpha = 0.05)
mat4p_high <- CreateAdjMat(ccd_4p_high, 4)

# Estimate PAG
pag4p <- plotPAG(ccd_4p_high, mat4p_high)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv4p_high <- semiequiv_dcg(ccd_4p_high, mat4p_high)
# save(equiv4p_high, file="data/equiv4p_high.RData")
load("data/equiv4p_high.RData")

## Density comparison
# true model density
truemoddensity(B4_high)
# GGM density
GGMdensity(ggm4p_high)
# average density for all equivalent DCGs
DCGdensity(equiv4p_high)

## Degree comparison
# true model degrees per node
GGMdegree(true4p_high)
# GGM  degrees per node
GGMdegree(ggm4p_high)
# average degree per node for all equivalent DCGs
DCGdegree(equiv4p_high)

## =============================================================================
## 3) Five nodes - sparse model (low density)
## =============================================================================

## Specify B matrix
# set the number of nodes (p)
p = 5
B5 = matrix(c(0, 1, 0, 0, 0,
              0, 0, 0, 0.7, 0,
              0, 0.4, 0, 0, 0,
              0, 0, .5, 0, 0,
              0, 0, 0, -1.5, 0), p, p, byrow = T)
colnames(B5) <- c("X1", "X2", "X3", "X4", "X5")

## Generate data
# first, equilibrium check (necessary condition for cyclic models to converge)
equilibrium_check(B5)
# generated data with N = 10^6, seed = 123
data5p <- gen_dat(B5, N =1e6, seed = 123)

## Specify layout
layout5 = matrix(c(0,1,
                   0,0,
                   1,-1,
                   2,0,
                   2,1),5,2,byrow = T)

## True cyclic graph
true5p <- qgraph(t(B5), layout=layout5, labels = colnames(B5), theme="colorblind")

## Estimate GGM
ggm5p <- qgraph(cor(data5p), layout = layout5, theme="colorblind")

## Run CCD algorithm
ccd_5p <- ccdKP(df=data5p, dataType = "continuous", alpha = 0.05)
mat5p <- CreateAdjMat(ccd_5p, 5)

# Estimate PAG
pag5p <- plotPAG(ccd_5p, mat5p)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv5p <- semiequiv_dcg(ccd_5p, mat5p)
# save(equiv5p, file="data/equiv5p.RData")
load("data/equiv5p.RData")

## Density comparison
# true model density
truemoddensity(B5)
# GGM density
GGMdensity(ggm5p)
# average density for all equivalent DCGs
DCGdensity(equiv5p)

## Degree comparison
# true model degrees per node
GGMdegree(true5p)
# GGM  degrees per node
GGMdegree(ggm5p)
# average degree per node for all equivalent DCGs
DCGdegree(equiv5p)

## =============================================================================
## 4) Five nodes - dense model (high density)
## =============================================================================

## Specify B matrix
# set the number of nodes (p)
p = 5

B5_high = matrix(c(0, 0.9, 0, 0, 0.6,
              0, 0, 0, 0.7, 0,
              0, 0.9, 0, 0, 0,
              0, 0, 0.5, 0, 0,
              0, 0, 0, 1, 0), p, p, byrow = T)
colnames(B5_high) <- c("X1", "X2", "X3", "X4", "X5")

## Generate data
# first, equilibrium check (necessary condition for cyclic models to converge)
equilibrium_check(B5_high)
# generated data with N = 10^6, seed = 1
data5p_high <- gen_dat(B5_high, N =1e6, seed = 1)

## Specify layout
layout5 = matrix(c(0,1,
                   0,0,
                   1,-1,
                   2,0,
                   2,1),5,2,byrow = T)

## True cyclic graph
true5p_high <- qgraph(t(B5_high), layout=layout5, labels = colnames(B5_high), theme="colorblind")

## Estimate GGM
ggm5p_high <- qgraph(cor(data5p_high), layout = layout5, theme="colorblind")

## Run CCD algorithm
ccd_5p_high <- ccdKP(df=data5p_high, dataType = "continuous", alpha = 0.05)
mat5p_high <- CreateAdjMat(ccd_5p_high, 5)

## Estimate PAG
pag5p_high <- plotPAG(ccd_5p_high, mat5p_high)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv5p_high <- semiequiv_dcg(ccd_5p_high, mat5p_high)
# save(equiv5p_high, file="data/equiv5p_high.RData")
load("data/equiv5p_high.RData")

## Density comparison
# true model density
truemoddensity(B5_high)
# GGM density
GGMdensity(ggm5p_high)
# average density for all equivalent DCGs
DCGdensity(equiv5p_high)

## Degree comparison
# true model degrees per node
GGMdegree(true5p_high)
# GGM  degrees per node
GGMdegree(ggm5p_high)
# average degree per node for all equivalent DCGs
DCGdegree(equiv5p_high)

## =============================================================================
## 5) Six nodes - sparse model (low density)
## =============================================================================

## Specify B matrix
# set the number of nodes (p)
p = 6
B6 = matrix(c(0, 0, 0, 0, 0, 0,
              0.3, 0, 0.4, 0, 0, 0,
              0, 0, 0, 0.9, 0, 0,
              0, 0, 0, 0, 0.4, 0,
              0, 0, 1, 0, 0, 0,
              1, 0, 0, 0, 0.5, 0), p, p, byrow = T)
colnames(B6) <- c("X1", "X2", "X3", "X4", "X5", "X6")

## Generate data
# first, equilibrium check (necessary condition for cyclic models to converge)
equilibrium_check(B6)
# generated data with N = 10^6, seed = 123
data6p <- gen_dat(B6, N =1e6, seed = 123)

## Specify layout
layout6 = matrix(c(1, 2,
                   0,1,
                   0,0,
                   1,-1,
                   2,0,
                   2,1),6,2,byrow = T)

## True cyclic graph
true6p <- qgraph(t(B6), layout=layout6, labels = colnames(B6), theme="colorblind")

## Estimate GGM
ggm6p <- qgraph(cor(data6p), layout = layout6, theme="colorblind")

## Run CCD algorithm
ccd_6p <- ccdKP(df=data6p, dataType = "continuous", alpha = 0.05)
mat6p <- CreateAdjMat(ccd_6p, 6)

## Estimate PAG
pag6p <- plotPAG(ccd_6p, mat6p)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv6p <- semiequiv_dcg(ccd_6p, mat6p)
# save(equiv6p, file="data/equiv6p.RData")
load("data/equiv6p.RData")

## Density comparison
# true model density
truemoddensity(B6)
# GGM density
GGMdensity(ggm6p)
# average density for all equivalent DCGs
DCGdensity(equiv6p)

## Degree comparison
# true model degrees per node
GGMdegree(true6p)
# GGM  degrees per node
GGMdegree(ggm6p)
# average degree per node for all equivalent DCGs
DCGdegree(equiv6p)

## =============================================================================
## 6) Six nodes - dense model (high density)
## =============================================================================

## Specify B matrix
# set the number of nodes (p)
p = 6
B6_high = matrix(c(0, 0, 0, 0, 0, 0,
              0.7, 0, 0.4, 0, 0, 0.9,
              0, 0, 0, 0.9, 0, 0,
              0, 0, 0, 0, 0.4, 0,
              0, 0, 1, 0, 0, 0,
              1, 0, 0, 0, 0.5, 0), p, p, byrow = T)
# colnames for B matrix is necessary for running CCD
colnames(B6_high) <- c("X1", "X2", "X3", "X4", "X5", "X6")

## Generate data
# first, equilibrium check (necessary condition for cyclic models to converge)
equilibrium_check(B6_high)
# generated data with N = 10^6, seed = 123
data6p_high<- gen_dat(B6_high, N =1e6, seed = 123)

## True cyclic graph
# we use the layout specified earlier in the 6p-sparse model.
true6p_high <- qgraph(t(B6_high), layout=layout6, labels = colnames(B6), theme="colorblind")

## Estimate GGM
ggm6p_high <- qgraph(cor(data6p_high), layout = layout6, theme="colorblind")

## Run CCD algorithm
ccd_6p_high <- ccdKP(df=data6p_high, dataType = "continuous", alpha = 0.05)
mat6p_high <- CreateAdjMat(ccd_6p_high, 6)

## Estimate PAG
pag6p_high <- plotPAG(ccd_6p_high, mat6p_high)


## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv6p_high <- semiequiv_dcg(ccd_6p_high, mat6p_high)
# save(equiv6p_high, file="data/equiv6p_high.RData")
load("data/equiv6p_high.RData")

## Density comparison
# true model density
truemoddensity(B6_high)
# GGM density
GGMdensity(ggm6p_high)
# average density for all equivalent DCGs
DCGdensity(equiv6p_high)

## Degree comparison
# true model degrees per node
GGMdegree(true6p_high)
# GGM  degrees per node
GGMdegree(ggm6p_high)
# average degree per node for all equivalent DCGs
DCGdegree(equiv6p_high)

## =============================================================================
## 7) Create density plots for all of the models
#
# For a comparison, we create two density plots:
# one for the sparse condition & the other one for the dense condition.
## =============================================================================

## 1. Combine densities for each type of model together for plotting
# density for true models
trueden <- list(B4, B4_high, B5, B5_high, B6, B6_high) %>%
  map( ~ truemoddensity(.)) %>% unlist() %>%  as.data.frame() %>% rename("TRUE"=".")
# density for GGM
ggmden <- list(ggm4p, ggm4p_high, ggm5p, ggm5p_high, ggm6p, ggm6p_high) %>%
  map( ~ GGMdensity(.)) %>% unlist() %>%  as.data.frame() %>% rename("GGM"=".")
# density for true DCG
dcgden <- list(equiv4p, equiv4p_high, equiv5p, equiv5p_high, equiv6p, equiv6p_high) %>%
  map( ~ DCGdensity(.)) %>% unlist() %>%  as.data.frame() %>% rename("DCG"=".")
# bind them together
modeldensities <- bind_cols(trueden, ggmden, dcgden) %>%
  set_rownames(c("4p-sparse", "4p-dense", "5p-sparse", "5p-dense", "6p-sparse", "6p-dense")) %>%
  mutate(model = rownames(.)) %>% pivot_longer(!model, names_to = "type", values_to="density")

## 2. Create plots
# 2-1) Specify my custom theme that I will use for the remainder
MyTheme <-  theme(plot.title = element_blank(),
                  plot.subtitle = element_text(face = "italic", family = "Palatino"),
                  axis.text=element_text(face = "bold"),
                  legend.text = element_text(face = "bold"))

# 2-2) Create density plots for each condition
# select the sparse conditions
low_densities <- modeldensities %>% filter(model %in% c("4p-sparse", "5p-sparse", "6p-sparse"))
# select the dense conditions
high_densities <- modeldensities %>% filter(model %in% c("4p-dense", "5p-dense", "6p-dense"))
# map each objects to the plotting function
densityplots <- list(low_densities, high_densities) %>%
  map( ~ ggplot(., aes(x=model, y=density, group = type, colour = type)) +
      scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
      geom_line(aes(group = type)) +
      geom_point() +
      labs(x="", y="", title = "", subtitle = paste(str_split(.$model, "-")[[1]][2], "condition")) +
      theme_classic() + MyTheme)
# combine plots
ggarrange(plotlist = densityplots, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(top = text_grob("Comparing the overall density of different models",
                                  face = "bold", family = "Palatino"))

## =============================================================================
## 8) Create degree centrality plots for all of the models
#
# For comparisons, we create the degree centrality plot for each model,
# which results in 3 (p = 4,5,6) by 2 (sparse/dense) = 6 plots in total.
## =============================================================================

## 1. Combine degrees for each type of model together for plotting
# list of true model
trueobj <- list(true4p, true4p_high, true5p, true5p_high, true6p, true6p_high)
# list of GGM
ggmobj <- list(ggm4p, ggm4p_high, ggm5p, ggm5p_high, ggm6p, ggm6p_high)
# list of DCG
dcgobj <- list(equiv4p, equiv4p_high, equiv5p, equiv5p_high, equiv6p, equiv6p_high)
modelnames <- c("4node-sparse","4node-dense","5node-sparse",
                "5node-dense","6node-sparse","6node-dense")
# compute the degree for all models and combine them all in a list
deglist <- list()
for(i in seq_along(trueobj)){
  # bind degrees for true, GGM, and DCG per each case
  deglist[[i]] <- bind_cols(GGMdegree(trueobj[[i]]),GGMdegree(ggmobj[[i]]),
                            DCGdegree(dcgobj[[i]])) %>%
    select(contains(c("node...1", "degree"))) %>%
    # rename the columns
    rename(node = node...1, truedegree = degree...2, ggmdegree = degree...4 ,
           dcg_avgdegree = average_degree) %>%
    # convert it to a long-format
    pivot_longer(!node, names_to = "model", values_to = "degree") %>%
    mutate(name = modelnames[i]) %>%
    # suppress messages from renaming columns step
    suppressMessages()
}

## 2. Create degree centrality plots
degplots <- deglist %>%
  # map the list containing all degree info to the plotting function
  map(~ ggplot(data = ., aes(x = degree, y = node, group = model, colour = model)) +
        geom_point() + geom_path(aes(group = model)) +
        labs(x="", y="", subtitle=.$name[1]) +
        scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name = "",
                            labels = c("DCG", "GGM", "TRUE")) + theme_bw() + MyTheme)
# combine plots
ggarrange(plotlist = degplots,
          ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(top = text_grob("Comparing the degree of different models",
                                  face = "bold", family = "Palatino"))
