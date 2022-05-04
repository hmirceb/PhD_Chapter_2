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
fp<-read.xlsx('~/Dropbox/DATA__LAB/__FLORA/FLORAPYR_REFERENCIAL/FLORAPYR_REFERENCIAL_02122020.xlsx',
              sheet = 5)
load("/home/hector/Dropbox/DATA__LAB/Hector_tesis/Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
rm(tab_v2, tre)
angios<-read.csv('~/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.csv')
angios<-angios %>% dplyr::select(-X)

ellenberg<-read.csv('C:/Users/18172844S/Dropbox/Tesis/Rareza/Datos/EIV.csv')
load("C:/Users/18172844S/Dropbox/DATA__LAB/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
fp<-read.xlsx('C:/Users/18172844S/Dropbox/DATA__LAB/__FLORA/FLORAPYR_REFERENCIAL/FLORAPYR_REFERENCIAL_02122020.xlsx',
              sheet = 5)
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
rm(tab_v2, tre)
angios<-read.csv('C:/Users/18172844S/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.csv')
angios<-angios %>% dplyr::select(-X)

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
  left_join(dplyr::select(angios, -TAXON_RP), 'GENERO_RP') %>%
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
         LA = scale(LA))
rownames(Datos_signal)<-Datos_signal$TAXON_REF_PYR_MOD


data<-Datos_signal %>%
  dplyr::select(c(Endemic, RA_max, LA, B_h)) %>%
  mutate(RA_max = as.numeric(RA_max)) %>%
  mutate(Endemic = ifelse(Endemic == 'End', 1, 0))



### Contribucion por nodos todo
#### IMPORTANTE !!!!!!!! #####
# La contribucion de los nodos impica que si un nodo contribuye mucho entonces las especies
# que descienden de el SON MUY DIFERENTES (lo contrario a la senyal filogenetica). Si los 
# nodos que contribuyen mucho estan cerca de la raiz significa que toda la diversidad esta
# en los clados principales y por tanto las especies luego se parecen mas entre si. Por
# contrario, si los nodos que contribuyen mas estan cerca de los tips es que hay mucha variacion
# y esas especies se parecen menos
# Tiene 5 opciones (Donde M = numero de linajes que descienden de un nodo:
# 1. Peso de los nodos/(diversidad entre nodos - diversidad intra nodos)
# 2. (Diversidad entre nodos - diversidad intra nodos)/(1- diversidad intranodos) * M/(M - 1)
# 3. (Diversidad entre nodos - diversidad intra nodos)/(1- diversidad entre nodos)/(M- 1)
# 4. (Diversidad entre nodos - diversidad intra nodos)
# 5. Peso de los nodos

# Hay que correr el rtestdecdiv. Los resultados son:
#   stat1 (S1): Toda la diversidad esta en un solo nodo
#   stat2 (S2): La contribucion se distribuye de manera uniforme entre nodos
#   stat3 (S3): Los nodos que mas contribuyen se situan mas o menos cerca de la raiz

# En el caso de stat3 hay que elegir como se ordenan los nodos, por complejidad o por distancia a la raiz

runs = 99

trees_all<-lapply(trees_sp_yule_10, function(x) 
  keep.tip(x, rownames(data)))


End_node_all_test<-lapply(trees_all,
                          function(x)
                            rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                   x$tip.label),
                                        dis = dist(dplyr::select(data, Endemic)),
                                        formula = 'EDI',
                                        vranking = 'droot',
                                        nrep = runs))

RA_node_all_test<-lapply(trees_all,
                         function(x)
                           rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                  x$tip.label),
                                       dis = dist(dplyr::select(data, RA_max)),
                                       formula = 'EDI',
                                       vranking = 'droot',
                                       nrep = runs))

HS_node_all_test<-lapply(trees_all,
                         function(x)
                           rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                  x$tip.label),
                                       dis = dist(dplyr::select(data, B_h)),
                                       formula = 'EDI',
                                       vranking = 'droot',
                                       nrep = runs))

LA_node_all_test<-lapply(trees_all,
                         function(x)
                           rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                  x$tip.label),
                                       dis = dist(dplyr::select(data, LA)),
                                       formula = 'EDI',
                                       vranking = 'droot',
                                       nrep = runs))

#### Señal filogenetica angios ####

Datos_signal_angios<-eiv %>% filter(Group == 'angiosperms') %>%
  dplyr::select(c(TAXON_REF_PYR_MOD, RA, RA_max, LA, Endemic, n_habs, B, B_h, PS)) %>%
  filter(complete.cases(.)) %>%
  mutate(RA_max = scale(as.numeric(RA_max)),
         B_h = scale(B_h),
         LA = scale(LA))
rownames(Datos_signal_angios)<-Datos_signal_angios$TAXON_REF_PYR_MOD

data_angios<-Datos_signal_angios %>%
  dplyr::select(c(Endemic, RA_max, LA, B_h)) %>%
  mutate(RA_max = as.numeric(RA_max)) %>%
  mutate(Endemic = ifelse(Endemic == 'End', 1, 0))



### Contribucion por nodos angios
trees_ang<-lapply(trees_sp_yule_10, function(x) 
  keep.tip(x, rownames(data_angios)))

End_node_ang_test<-lapply(trees_ang,
                         function(x)
                           rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                  x$tip.label),
                                       dis = dist(dplyr::select(data_angios, Endemic)),
                                       formula = 'EDI',
                                       vranking = 'droot',
                                       nrep = runs))

RA_node_ang_test<-lapply(trees_ang,
                         function(x)
                           rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                  x$tip.label),
                                       dis = dist(dplyr::select(data_angios, RA_max)),
                                       formula = 'EDI',
                                       vranking = 'droot',
                                       nrep = runs))

HS_node_ang_test<-lapply(trees_ang,
                          function(x)
                            rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                   x$tip.label),
                                        dis = dist(dplyr::select(data_angios, B_h)),
                                        formula = 'EDI',
                                        vranking = 'droot',
                                        nrep = runs))

LA_node_ang_test<-lapply(trees_ang,
                         function(x)
                           rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                  x$tip.label),
                                       dis = dist(dplyr::select(data_angios, LA)),
                                       formula = 'EDI',
                                       vranking = 'droot',
                                       nrep = runs))


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



# Contribucion por nodos gimnos
trees_gimn<-lapply(trees_sp_yule_10, function(x) 
  keep.tip(x, rownames(data_gimn)))


End_node_gimn_test<-lapply(trees_gimn,
                           function(x)
                             rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                   x$tip.label),
                                         dis = dist(dplyr::select(data_gimn, Endemic)),
                                         formula = 'EDI',
                                         vranking = 'droot',
                                         nrep = runs))


RA_node_gimn_test<-lapply(trees_gimn,
                           function(x)
                             rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                    x$tip.label),
                                         dis = dist(dplyr::select(data_gimn, RA_max)),
                                         formula = 'EDI',
                                         vranking = 'droot',
                                         nrep = runs))


HS_node_gimn_test<-lapply(trees_gimn,
                          function(x)
                            rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                   x$tip.label),
                                        dis = dist(dplyr::select(data_gimn, B_h)),
                                        formula = 'EDI',
                                        vranking = 'droot',
                                        nrep = runs))

LA_node_gimn_test<-lapply(trees_gimn,
                          function(x)
                            rtestdecdiv(phyl = x, vecab = setNames(rep(1, times = length(x$tip.label)),
                                                                   x$tip.label),
                                        dis = dist(dplyr::select(data_gimn, LA)),
                                        formula = 'EDI',
                                        vranking = 'droot',
                                        nrep = runs))

save(End_node_all_test, End_node_ang_test, End_node_gimn_test,
     RA_node_all_test, RA_node_ang_test, RA_node_gimn_test,
     HS_node_all_test, HS_node_ang_test, HS_node_gimn_test,
     LA_node_all_test, LA_node_ang_test, LA_node_gimn_test,
     file = 'rtestdecdiv.RData')