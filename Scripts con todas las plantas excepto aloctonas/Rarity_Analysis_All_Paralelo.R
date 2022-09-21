library(tidyverse)
library(phylosignal)
library(phylobase)
library(ape)
library(openxlsx)
library(vegan)
library(corrgram)
library(caper)
library(popbio)
library(adiv)
library(picante)
library(phytools)
library(doParallel)
library(parallel)
library(foreach)

# Cargamos los datos
ellenberg<-read.csv('~/Dropbox/Tesis/Rareza/Datos/EIV.csv')
load("/home/hector/Dropbox/DATA__LAB/PHYLO/de_Hector/1000_ArbolesPirineo_31_05_2022.RData")
fp<-read.xlsx('~/Dropbox/DATA__LAB/__FLORA/FLORAPYR_REFERENCIAL/FLORAPYR_REFERENCIAL_10022021.xlsx',
              sheet = 6)
load("/home/hector/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
rm(tab_v2, tre)
load("~/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")

ellenberg<-read.csv('C:/Users/18172844S/Dropbox/Tesis/Rareza/Datos/EIV.csv')
load("C:/Users/18172844S/Dropbox/DATA__LAB/__FLORA/PHYLO/de_Hector/1000_ArbolesPirineo_31_05_2022.RData")
fp<-read.xlsx('C:/Users/18172844S/Dropbox/DATA__LAB/__FLORA/FLORAPYR_REFERENCIAL/FLORAPYR_REFERENCIAL_17032022.xlsx',
              sheet = 6)
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
rm(tab_v2, tre)
load("C:/Users/18172844S/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")

names(ellenberg)[6]<-'TAXON_REF_PYR_MOD'
eiv<-ellenberg[,6:15]
eiv<-eiv[eiv$RANGO_RP == 'ES',]

names(fp)[6]<-'TAXON_REF_PYR_MOD'
eiv<-eiv %>% left_join(dplyr::select(fp, TAXON_REF_PYR_MOD, GENERO_RP, UTM10_RP,
                                     UTM10_SP_FP, ENDEMICA_PIRINEOS, HABITAT_FP)) %>%
  left_join(Rareza) %>% 
  left_join(dplyr::select(fp_fams, -TAXON_RP), 'GENERO_RP') %>%
  filter(!duplicated(TAXON_REF_PYR_MOD)) %>% 
  mutate(ENDEMICA_PIRINEOS = ifelse(is.na(ENDEMICA_PIRINEOS), 'NoEnd', 'End')) %>%
  filter(complete.cases(RA)) %>%
  filter(TAXON_REF_PYR_MOD %in% trees_sp_yule_1000[[1]]$tip.label)

eiv$UTM10_RP<-ifelse(is.na(eiv$UTM10_RP), eiv$RA, eiv$UTM10_RP)
eiv$UTM10_SP_FP<-ifelse(is.na(eiv$UTM10_SP_FP), eiv$RA, eiv$UTM10_SP_FP)

eiv$RA_max<-apply(dplyr::select(eiv, c(UTM10_RP, UTM10_SP_FP, RA)), 1, max) 


############
# Señal filogenetica en todo el arbol 

Datos_signal<-eiv %>% 
  dplyr::select(c(TAXON_REF_PYR_MOD, RA, RA_max, LA, Endemic, n_habs, B, B_h, PS)) %>%
  filter(complete.cases(.)) %>%
  mutate(RA_max = scale(as.numeric(RA_max)),
         B_h = scale(B_h),
         LA = scale(LA), 
         Endemic = factor(Endemic, levels = c('NoEnd', 'End')))
rownames(Datos_signal)<-Datos_signal$TAXON_REF_PYR_MOD


# Montamos el cluster
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)
n_trees<-length(trees_sp_yule_1000)/100

Endemic_PhyloD<-foreach(K = 1:n_trees,
                       .packages=c('phangorn',
                                   'tidyverse',
                                   'castor',
                                   'caper',
                                   'picante')) %dopar% {
                                     comp<-comparative.data(keep.tip(trees_sp_yule_1000[[K]], rownames(Datos_signal)), 
                                                            Datos_signal, 
                                                            TAXON_REF_PYR_MOD)
                                     phylo.d(comp, binvar = Endemic)}
stopCluster(cl)

data<-Datos_signal %>%
  dplyr::select(c(Endemic, RA_max, LA, B_h)) %>%
  mutate(RA_max = as.numeric(RA_max)) %>%
  mutate(Endemic = ifelse(Endemic == 'End', 1, 0))

cl <- makeCluster(no_cores)
registerDoParallel(cl)

filo_all<-foreach(K = 1:n_trees,
                   .packages=c('phangorn',
                               'tidyverse',
                               'castor',
                               'caper',
                               'picante',
                               'phylobase',
                               'phylosignal')) %dopar% {
                                 tree<-phylo4d(keep.tip(trees_sp_yule_1000[[K]], rownames(data)),
                                               tip.data = data)
                                 signal_all<-phyloSignal(tree, 999, method = 'Lambda')
                                 
                                 lipa_all<-lipaMoran(tree)
                                 list(signal_all, lipa_all)
                               }
stopCluster(cl)

signal_all<-lapply(filo_all, function(x) x[[1]])
lipa_all<-lapply(filo_all, function(x) x[[2]])

#### Señal filogenetica angios ####

Datos_signal_angios<-eiv %>% filter(Group == 'angiosperms') %>%
  dplyr::select(c(TAXON_REF_PYR_MOD, RA, RA_max, LA, Endemic, n_habs, B, B_h, PS)) %>%
  filter(complete.cases(.)) %>%
  mutate(RA_max = scale(as.numeric(RA_max)),
         B_h = scale(B_h),
         LA = scale(LA))
rownames(Datos_signal_angios)<-Datos_signal_angios$TAXON_REF_PYR_MOD

cl <- makeCluster(no_cores)
registerDoParallel(cl)

Endemic_PhyloD_ang<-foreach(K = 1:n_trees,
                        .packages=c('phangorn',
                                    'tidyverse',
                                    'castor',
                                    'caper',
                                    'picante')) %dopar% {
                                      comp<-comparative.data(keep.tip(trees_sp_yule_1000[[K]], rownames(Datos_signal_angios)), 
                                                             Datos_signal_angios, 
                                                             TAXON_REF_PYR_MOD)
                                      phylo.d(comp, binvar = Endemic)}
stopCluster(cl)

cl <- makeCluster(no_cores)  
registerDoParallel(cl)

filo_angios<-foreach(K = 1:n_trees,
                  .packages=c('phangorn',
                              'tidyverse',
                              'castor',
                              'caper',
                              'picante',
                              'phylobase',
                              'phylosignal')) %dopar% {
                                tree<-phylo4d(keep.tip(trees_sp_yule_1000[[K]], rownames(data_angios)),
                                              tip.data = data_angios)
                                signal_angios<-phyloSignal(tree, 999, method = 'Lambda')
                                
                                lipa_angios<-lipaMoran(tree)
                                list(signal_angios, lipa_angios)
                              }
stopCluster(cl)

signal_angios<-lapply(filo_angios, function(x) x[[1]])
lipa_angios<-lapply(filo_angios, function(x) x[[2]])


#### Señal filogenetica gimnos + pteridofitos ####
# Señal + LIPA
Datos_signal_gimn<-eiv %>% filter(Group != 'angiosperms') %>%
  dplyr::select(c(TAXON_REF_PYR_MOD, RA_max, LA, Endemic, B_h)) %>%
  filter(complete.cases(.)) %>%
  mutate(RA_max = scale(as.numeric(RA_max)),
         B_h = scale(B_h),
         LA = scale(LA))
rownames(Datos_signal_gimn)<-Datos_signal_gimn$TAXON_REF_PYR_MOD

data_gimn<-Datos_signal_gimn %>%
  dplyr::select(-c(TAXON_REF_PYR_MOD)) %>%
  mutate(RA_max = as.numeric(RA_max)) %>%
  mutate(Endemic = ifelse(Endemic == 'End', 1, 0))

cl <- makeCluster(no_cores)  
registerDoParallel(cl)

signal_gimn<-list()
lipa_gimn<-list()
filo_gimn<-foreach(K = 1:n_trees,
        .packages=c('phangorn',
                    'tidyverse',
                    'castor',
                    'caper',
                    'picante',
                    'phylobase',
                    'phylosignal')) %dopar% {
                      tree<-phylo4d(keep.tip(trees_sp_yule_1000[[K]], rownames(data_gimn)),
                                    tip.data = data_gimn)
                      signal_gimn<-phyloSignal(tree, 999, method = 'Lambda')
                      
                      lipa_gimn<-lipaMoran(tree)
                      list(signal_gimn, lipa_gimn)
                    }
stopCluster(cl)

signal_gimn<-lapply(filo_gimn, function(x) x[[1]])
lipa_gimn<-lapply(filo_gimn, function(x) x[[2]])


#### Filo PCA ####
trees<-lapply(trees_sp_yule_1000, function(x)
  keep.tip(x, Datos_signal$TAXON_REF_PYR_MOD))

dat_pca<-Datos_signal %>%
  dplyr::select(RA_max, B_h, LA)
rownames(dat_pca)<-Datos_signal$TAXON_REF_PYR_MOD

pca_filo_trees<-lapply(trees[1:3], function(x)
  phyl.pca(tree = x, method = 'lambda', dat_pca))

pca_filo_trees_null<-lapply(trees[1], function(x)
  phyl.pca(tree = x, method = 'lambda', dat_pca, 
           opt = 'fixed', 
           lambda = 0))

save(eiv, Endemic_PhyloD, signal_all, lipa_all,
     Endemic_PhyloD_ang, signal_angios, lipa_angios,
     signal_gimn, lipa_gimn, pca_filo_trees,
     file = 'FiloSignal_Rareza_scaled_All.RData')
