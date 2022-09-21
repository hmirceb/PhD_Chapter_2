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

# Cargamos los datos
ellenberg<-read.csv('~/Dropbox/Tesis/Rareza/Datos/EIV.csv')
load("/home/hector/Dropbox/DATA__LAB/__FLORA/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
fp<-read.xlsx('~/Dropbox/DATA__LAB/__FLORA/FLORAPYR_REFERENCIAL/FLORAPYR_REFERENCIAL_17032022.xlsx',
              sheet = 6)
load("/home/hector/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
rm(tab_v2, tre)
load("~/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")

ellenberg<-read.csv('C:/Users/18172844S/Dropbox/Tesis/Rareza/Datos/EIV.csv')
load("C:/Users/18172844S/Dropbox/DATA__LAB/__FLORA/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
fp<-read.xlsx('C:/Users/18172844S/Dropbox/DATA__LAB/__FLORA/FLORAPYR_REFERENCIAL/FLORAPYR_REFERENCIAL_17032022.xlsx',
              sheet = 6)
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
rm(tab_v2, tre)
load("C:/Users/18172844S/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")

# Cambiamos nombre de la variable correspondiente a las especies,
# nos quedamos solo con especies completas y quitamos los duplicados
names(ellenberg)[6]<-'TAXON_REF_PYR_MOD'
eiv<-ellenberg[,6:15]
eiv<-eiv[eiv$RANGO_RP == 'ES',]

# Unimos la tabla de FLORAPYR con la de valores de Ellenberg y los de rareza
names(fp)[6]<-'TAXON_REF_PYR_MOD'
eiv<-eiv %>% left_join(dplyr::select(fp, TAXON_REF_PYR_MOD, GENERO_RP, UTM10_RP,
                              UTM10_SP_FP, ENDEMICA_PIRINEOS, HABITAT_FP)) %>%
  left_join(Rareza) %>% 
  left_join(dplyr::select(fp_fams, -TAXON_RP), 'GENERO_RP') %>%
  filter(!duplicated(TAXON_REF_PYR_MOD)) %>% 
  mutate(ENDEMICA_PIRINEOS = ifelse(is.na(ENDEMICA_PIRINEOS), 'NoEnd', 'End')) %>%
  filter(complete.cases(RA)) %>%
  filter(TAXON_REF_PYR_MOD %in% trees_sp_yule_10[[1]]$tip.label)

# Como el numero de UTM10 varia entre los inventarios y el FP, nos quedamos con el valor
  # mas alto
eiv$UTM10_RP<-ifelse(is.na(eiv$UTM10_RP), eiv$RA, eiv$UTM10_RP)
eiv$UTM10_SP_FP<-ifelse(is.na(eiv$UTM10_SP_FP), eiv$RA, eiv$UTM10_SP_FP)

eiv$RA_max<-as.numeric(apply(dplyr::select(eiv, c(UTM10_RP, UTM10_SP_FP, RA)), 1, max))


# Quitamos las especies que aparezcan en menos de 3 inventarios (singletons y doubletons) y las que no esten en el arbol (Xatardia scabra basicamente)
#eiv<-eiv %>%
 # filter(n_ocurr >= 3) %>%
  #filter(TAXON_REF_PYR_MOD %in% trees_sp_yule_10[[1]]$tip.label)

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

# Señal filogenetica de endemismo con el indice D para variables dicotomicas

Endemic_PhyloD<-list()
for (i in 1:length(trees_sp_yule_10)) {
  comp<-comparative.data(keep.tip(trees_sp_yule_10[[i]], rownames(Datos_signal)), 
                         Datos_signal, 
                         TAXON_REF_PYR_MOD)
  Endemic_PhyloD[[i]]<-phylo.d(comp, binvar = Endemic)
}


data<-Datos_signal %>%
  dplyr::select(c(Endemic, RA_max, LA, B_h)) %>%
  mutate(RA_max = as.numeric(RA_max)) %>%
  mutate(Endemic = ifelse(Endemic == 'End', 1, 0))
  

# LIPA (Local Indicator of Phylogenetic Association) es una extension de la 
  #I de Moran que comprueba la correlacion rama a rama
  # y dice en que zonas las especies estan mas correlacionadas
signal_all<-list()
lipa_all<-list()

for (j in 1:length(trees_sp_yule_10)) {
  tree<-phylo4d(keep.tip(trees_sp_yule_10[[j]], rownames(data)),
                tip.data = data)
  # Senyal con todos los indices
  signal_all[[j]]<-phyloSignal(tree, 999, method = 'Lambda')
  
  # LIPA por cada indice
  lipa_all[[j]]<-lipaMoran(tree)
  
}


#### Señal filogenetica angios ####

Datos_signal_angios<-eiv %>% filter(Group == 'angiosperms') %>%
  dplyr::select(c(TAXON_REF_PYR_MOD, RA, RA_max, LA, Endemic, n_habs, B, B_h, PS)) %>%
  filter(complete.cases(.)) %>%
  mutate(RA_max = scale(as.numeric(RA_max)),
         B_h = scale(B_h),
         LA = scale(LA))
rownames(Datos_signal_angios)<-Datos_signal_angios$TAXON_REF_PYR_MOD

Endemic_PhyloD_ang<-list()
for (i in 1:length(trees_sp_yule_10)) {
  comp<-comparative.data(keep.tip(trees_sp_yule_10[[i]], rownames(Datos_signal_angios)), 
                         Datos_signal_angios, 
                         TAXON_REF_PYR_MOD)
  Endemic_PhyloD_ang[[i]]<-phylo.d(comp, binvar = Endemic)
}

data_angios<-Datos_signal_angios %>%
  dplyr::select(c(Endemic, RA_max, LA, B_h)) %>%
  mutate(RA_max = as.numeric(RA_max)) %>%
  mutate(Endemic = ifelse(Endemic == 'End', 1, 0))

signal_angios<-list()
lipa_angios<-list()

for (j in 1:length(trees_sp_yule_10)) {
  tree<-phylo4d(keep.tip(trees_sp_yule_10[[j]], rownames(data_angios)),
                tip.data = data_angios)
  signal_angios[[j]]<-phyloSignal(tree, 999, method = 'Lambda')
  
  lipa_angios[[j]]<-lipaMoran(tree)
  
}

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

signal_gimn<-list()
lipa_gimn<-list()

for (j in 1:length(trees_sp_yule_10)) {
  tree<-phylo4d(keep.tip(trees_sp_yule_10[[j]], rownames(data_gimn)),
                tip.data = data_gimn)
  signal_gimn[[j]]<-phyloSignal(tree, 999, method = 'Lambda')
  
  lipa_gimn[[j]]<-lipaMoran(tree)
  
}

#### Filo PCA ####
trees<-lapply(trees_sp_yule_10, function(x)
  keep.tip(x, Datos_signal$TAXON_REF_PYR_MOD))

dat_pca<-Datos_signal %>%
  dplyr::select(RA_max, B_h, LA)
rownames(dat_pca)<-Datos_signal$TAXON_REF_PYR_MOD

pca_filo_trees<-lapply(trees, function(x)
  phyl.pca(tree = x, method = 'lambda', dat_pca))

pca_filo_trees_null<-lapply(trees[1], function(x)
  phyl.pca(tree = x, method = 'lambda', dat_pca, 
           opt = 'fixed', 
           lambda = 0))

save(eiv, Endemic_PhyloD, signal_all, lipa_all,
     Endemic_PhyloD_ang, signal_angios, lipa_angios,
     signal_gimn, lipa_gimn, pca_filo_trees,
     file = 'FiloSignal_Rareza_scaled_All.RData')


