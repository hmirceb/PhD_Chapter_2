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

# Cargamos los datos
ellenberg<-read.csv('~/Dropbox/Tesis/Rareza/Datos/EIV.csv')
load("/home/hector/Dropbox/DATA__LAB/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
fp<-read.xlsx('~/Dropbox/DATA__LAB/__FLORA/FLORAPYR_REFERENCIAL/FLORAPYR_REFERENCIAL_10022021.xlsx',
              sheet = 5)
load("/home/hector/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
rm(tab_v2, tre)
load("~/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")

ellenberg<-read.csv('C:/Users/18172844S/Dropbox/Tesis/Rareza/Datos/EIV.csv')
load("C:/Users/18172844S/Dropbox/DATA__LAB/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
fp<-read.xlsx('C:/Users/18172844S/Dropbox/DATA__LAB/__FLORA/FLORAPYR_REFERENCIAL/FLORAPYR_REFERENCIAL_10022021.xlsx',
              sheet = 5)
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
  filter(complete.cases(RA))

# Como el numero de UTM10 varia entre los inventarios y el FP, nos quedamos con el valor
  # mas alto
eiv$UTM10_RP<-ifelse(is.na(eiv$UTM10_RP), eiv$RA, eiv$UTM10_RP)
eiv$UTM10_SP_FP<-ifelse(is.na(eiv$UTM10_SP_FP), eiv$RA, eiv$UTM10_SP_FP)

eiv$RA_max<-apply(dplyr::select(eiv, c(UTM10_RP, UTM10_SP_FP, RA)), 1, max)


# Quitamos las especies que aparezcan en menos de 3 inventarios (singletons y doubletons) y las que no esten en el arbol (Xatardia scabra basicamente)
eiv<-eiv %>%
  filter(n_ocurr >= 3) %>%
  filter(TAXON_REF_PYR_MOD %in% trees_sp_yule_10[[1]]$tip.label)

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


trees<-lapply(trees_sp_yule_10, function(x) keep.tip(phy = x, tip = eiv$TAXON_REF_PYR_MOD))

EDs<-lapply(trees, function(x) evol.distinct(tree = x, type = 'fair.proportion'))

EDs<-data.frame('TAXON_REF_PYR_MOD' = EDs[[1]]$Species, 
                'ED' = apply(do.call('cbind', lapply(EDs, function(x) x$w)), 1, mean))

### Contribucion por nodos todo
#### IMPORTANTE !!!!!!!! #####
# La contribucion de los nodos impica que si un nodo contribuye mucho entonces las especies
# que descienden de el SON MUY DIFERENTES (lo contrario a la senyal filogenetica). Si los 
# nodos que contribuyen mucho estan cerca de la raiz significa que toda la diversidad esta
# en los clados principales y por tanto las especies luego se parecen mas entre si. Por
# contrario, si los nodos que contribuyen mas estan cerca de los tips es que hay mucha variacion
# y esas especies se parecen menos

# Hay que correr el rtestdecdiv. Los resultados son:
#   stat1 (S1): Toda la diversidad esta en un solo nodo
#   stat2 (S2): La contribucion se distribuye de manera uniforme entre nodos
#   stat3 (S3): Los nodos que mas contribuyen se situan mas o menos cerca de la raiz

# En el caso de stat3 hay que elegir como se ordenan los nodos, por complejidad o por distancia a la raiz


trees_all<-lapply(trees_sp_yule_10, function(x) 
  keep.tip(x, rownames(data)))
End_node_all<-lapply(trees_all,
                    function(x)
                      decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                       x$tip.label),
                             dis = dist(dplyr::select(data, Endemic)), 
                             formula = 'EDI'))

RA_node_all<-lapply(trees_all,
                    function(x)
                      decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                       x$tip.label),
                             dis = dist(dplyr::select(data, RA_max)),
                             formula = 'EDI'))

HS_node_all<-lapply(trees_all,
                    function(x)
                      decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                       x$tip.label),
                             dis = dist(dplyr::select(data, B_h)),
                             formula = 'EDI'))

LA_node_all<-lapply(trees_all,
                    function(x)
                      decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                       x$tip.label),
                             dis = dist(dplyr::select(data, LA)),
                             formula = 'EDI'))

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

### Contribucion por nodos angios
trees_ang<-lapply(trees_sp_yule_10, function(x) 
  keep.tip(x, rownames(data_angios)))
End_node_ang<-lapply(trees_ang,
                    function(x)
                      decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                       x$tip.label),
                             dis = dist(dplyr::select(data_angios, Endemic)),
                             formula = 'EDI'))

RA_node_ang<-lapply(trees_ang,
                     function(x)
                       decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                        x$tip.label),
                              dis = dist(dplyr::select(data_angios, RA_max)),
                              formula = 'EDI'))

HS_node_ang<-lapply(trees_ang,
                     function(x)
                       decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                        x$tip.label),
                              dis = dist(dplyr::select(data_angios, B_h)),
                              formula = 'EDI'))

LA_node_ang<-lapply(trees_ang,
                     function(x)
                       decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                        x$tip.label),
                              dis = dist(dplyr::select(data_angios, LA)),
                              formula = 'EDI'))


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


# Contribucion por nodos gimnos
trees_gimn<-lapply(trees_sp_yule_10, function(x) 
  keep.tip(x, rownames(data_gimn)))

End_node_gimn<-lapply(trees_gimn,
                     function(x)
                       decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                        x$tip.label),
                              dis = dist(dplyr::select(data_gimn, Endemic)),
                              formula = 'EDI'))

RA_node_gimn<-lapply(trees_gimn,
                     function(x)
                       decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                        x$tip.label),
                              dis = dist(dplyr::select(data_gimn, RA_max)),
                              formula = 'EDI'))

HS_node_gimn<-lapply(trees_gimn,
                     function(x)
                       decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                        x$tip.label),
                              dis = dist(dplyr::select(data_gimn, B_h)),
                              formula = 'EDI'))

LA_node_gimn<-lapply(trees_gimn,
                     function(x)
                       decdiv(phyl = x, comm = setNames(rep(1, times = length(x$tip.label)),
                                                        x$tip.label),
                              dis = dist(dplyr::select(data_gimn, LA)),
                              formula = 'EDI'))

save(eiv, EDs, Endemic_PhyloD, signal_all, lipa_all, End_node_all, RA_node_all, HS_node_all, LA_node_all,
     Endemic_PhyloD_ang, signal_angios, lipa_angios, End_node_ang, RA_node_ang, HS_node_ang, LA_node_ang,
     signal_gimn, lipa_gimn, End_node_gimn, RA_node_gimn, HS_node_gimn, LA_node_gimn,
     file = 'FiloSignal_Rareza_scaled.RData')
