library(tidyverse)
library(reshape2)
library(openxlsx)
library(corrgram)

load('C:/Users/18172844S/Dropbox/DATA__LAB/__VEGETACION_HABITATS/INVENTARIOS/PIRINEOS/SIVIM/REGISTROS_INV_SIVIM_100621.RData')
load('~/Dropbox/DATA__LAB/__VEGETACION_HABITATS/INVENTARIOS/PIRINEOS/SIVIM/REGISTROS_INV_SIVIM_100621.RData')
red_list<-read.xlsx('C:/Users/18172844S/Dropbox/DATA__LAB/__FLORA/CATALOGACION/CAT_FLORAPYR/Catalogadas_PIRINEOS/FLORA VASCULAR AMENAZADA DE PIRINEOS - LISTA ROJA 2019.xlsx')
red_list<-read.xlsx('~/Dropbox/DATA__LAB/__FLORA/CATALOGACION/CAT_FLORAPYR/Catalogadas_PIRINEOS/FLORA VASCULAR AMENAZADA DE PIRINEOS - LISTA ROJA 2019.xlsx')

SIVIM_rare<-SIVIM
SIVIM_rare<-SIVIM_rare[SIVIM_rare$EUNIS_LEVEL_2_GROUPED != '',] #Quitamos habitats antropicos y sin definir
SIVIM_rare<-SIVIM_rare[!is.na(SIVIM_rare$EUNIS_LEVEL_2_GROUPED),]
SIVIM_rare<-SIVIM_rare[SIVIM_rare$NATIVE_PYR != '',] # Quitamos aloctonas
SIVIM_rare<-SIVIM_rare[!endsWith(SIVIM_rare$TAXON_REF_PYR, 'sp.'),] # Quitamos especies en duda
SIVIM_rare<-SIVIM_rare[!is.na(SIVIM_rare$ABUNDANCE_TRANS),] 
SIVIM_rare$ABUNDANCE_TRANS_MOD<-ifelse(SIVIM_rare$ABUNDANCE_TRANS == 0, 1,
                                       SIVIM_rare$ABUNDANCE_TRANS)
SIVIM_rare$EUNIS_LEVEL_2_GROUPED <-ifelse(SIVIM_rare$EUNIS_LEVEL_2_GROUPED == 'Oromediterranean shrublands',
       'Garrigue', SIVIM_rare$EUNIS_LEVEL_2_GROUPED )
unique(SIVIM_rare$EUNIS_LEVEL_2_GROUPED)

# Contamos en cuantas UTM 10 aparece cada planta (Regional abundance RA) y creamos un df con esos datos
a<-aggregate(UTM10~TAXON_REF_PYR, SIVIM_rare, table)
Rareza<-data.frame('sp' = a$TAXON_REF_PYR, 'RA' = apply(a, 1, function(x) length(x[[2]])))

# Lo mismo pero para los habitats (n_habs)
b<-aggregate(EUNIS_LEVEL_2_GROUPED~TAXON_REF_PYR, SIVIM_rare, table)
Rareza$n_habs<-apply(b, 1, function(x) length(x[[2]]))

# Abundancia media de cada especie (Local Abundance, LA) y su intervalo de confianza
se<-function(x) sd(x, na.rm = T)/(sqrt(length(x)))
Rareza$LA<-aggregate(ABUNDANCE_TRANS_MOD~TAXON_REF_PYR, SIVIM_rare, mean)[,2]
Rareza$LA.se<-aggregate(ABUNDANCE_TRANS_MOD~TAXON_REF_PYR, SIVIM_rare, se)[,2]
Rareza$LA.se<-ifelse(is.na(Rareza$LA.se), 0, Rareza$LA.se)

# DF con las especies (solo especies, no subespecies)
Endemic<-as.data.frame(unique(cbind(SIVIM_rare$TAXON_REF_PYR, SIVIM_rare$ENDEMIC_PYR)))
Endemic$Endemic<-ifelse(Endemic$V2 == "Endemic_sp", 'End', 'NoEnd')
Endemic$Endemic<-ifelse(is.na(Endemic$Endemic), 'NoEnd', Endemic$Endemic)
names(Endemic)[1]<-'sp'

# Comprobamos que no haya ninguna especie que aparezca en mas inventarios que 
# habitats o UTMs
raras<-as.data.frame(table(SIVIM_rare$TAXON_REF_PYR))
names(raras)[1]<-'sp'
utms<-left_join(raras, Rareza, by = 'sp')
if(nrow(utms[utms$Freq < utms$RA,]) == nrow(utms[utms$Freq < utms$habspe,]) & 
   nrow(utms[utms$Freq < utms$RA,]) == 0){rm(utms, raras, a)}

# Unimos endemicas, rango y especificidad
Rareza<-merge(Rareza, Endemic, by = 'sp')
Rareza<-Rareza[,-which(names(Rareza) == 'V2')]
Rareza$RA_p<-100*(Rareza$RA/max(Rareza$RA))

rm(b, Endemic)
names(Rareza)[1]<-'TAXON_REF_PYR'


# Calculamos indices de amplitud de nicho segun Feinsinger 1981 (https://dx.doi.org/10.2307/1936664). 
    # B = Indice de Levins 1968 (Eq 1)
    # B_norm = B normalizado entre 0 y 1 (Eq 2)
    # B_h = Indice de Hurlbert 1978 (Eq 3)
    # PS = Indice de Feisinger 1981 (Eq 4)

especies<-unique(SIVIM_rare$TAXON_REF_PYR)
habitats<-list()

a_i<-as.data.frame(
  table(
    unique(
      data.frame('INV_CODE' = SIVIM_rare$INV_CODE, 
                 'EUNIS_LEVEL_2_GROUPED' = SIVIM_rare$EUNIS_LEVEL_2_GROUPED))$EUNIS_LEVEL_2_GROUPED
  )
)

a_i$A<-sum(a_i$Freq)
names(a_i)<-c('EUNIS_LEVEL_2_GROUPED', 'a_i', 'A')
a_i$q_i<-a_i$a_i/a_i$A

for (i in 1:length(especies)) {
  
  temp<-as.data.frame(
    table(
      SIVIM_rare[SIVIM_rare$TAXON_REF_PYR == especies[i],]$EUNIS_LEVEL_2_GROUPED)
  )
  
  temp<- if(dim(temp)[1] == 14) {temp} else {
    rbind(temp,
          data.frame('Var1' = a_i$EUNIS_LEVEL_2_GROUPED[!a_i$EUNIS_LEVEL_2_GROUPED %in% temp$Var1],
                     'Freq' = 0))}
  
  temp$TAXON_REF_PYR<-especies[i]
  temp$total<-sum(temp$Freq)
  names(temp)<-c('EUNIS_LEVEL_2_GROUPED', 'x_i', 'TAXON_REF_PYR', 'X')
  temp$p_i<-temp$x_i/temp$X
  temp<-left_join(temp, a_i, 'EUNIS_LEVEL_2_GROUPED')
  temp$B<-1/sum(temp$p_i^2)
  temp$B_norm<-1/(dim(a_i)[1]*sum(temp$p_i^2))
  temp$B_h<-1/sum(temp$p_i^2/temp$q_i)
  temp$PS<-1-0.5*sum(abs(temp$p_i-temp$q_i))
  temp$PS_alt<-sum(apply(dplyr::select(temp, c(p_i, q_i)), MARGIN = 1, min))
  
  habitats[[i]]<-temp
  rm(temp)
}
habitats<-do.call('rbind', habitats)

Rareza<-left_join(Rareza, unique(data.frame('TAXON_REF_PYR' = habitats$TAXON_REF_PYR,
                                          'B' = habitats$B, 
                                          'B_norm' = habitats$B_norm,
                                          'B_h' = habitats$B_h, 
                                          'PS' = habitats$PS,
                                          'PS_alt' = habitats$PS_alt)))
rm(a_i)

names(Rareza)[1]<-'TAXON_REF_PYR_MOD'
ocurrencias_sps<-data.frame(table(SIVIM_rare$TAXON_REF_PYR))
names(ocurrencias_sps)<-c('TAXON_REF_PYR_MOD', 'n_ocurr')
Rareza<-left_join(Rareza, ocurrencias_sps) 

red_list<-red_list %>%
  dplyr::select(c(TAXON_REF_FLORAPYR, CAT_AMENAZA_PYR)) %>%
  rename(TAXON_REF_PYR_MOD = TAXON_REF_FLORAPYR)

Rareza<-left_join(Rareza, red_list, 'TAXON_REF_PYR_MOD')

Frecuencia_por_habitat<-SIVIM_rare %>% 
  dplyr::select(c(TAXON_REF_PYR, EUNIS_LEVEL_2_GROUPED)) %>%
  group_by(TAXON_REF_PYR) %>%
  count(EUNIS_LEVEL_2_GROUPED) %>%
  rename(TAXON_REF_PYR_MOD = TAXON_REF_PYR)

SIVIM_rare %>%
  dplyr::select(c(INV_CODE, EUNIS_LEVEL_2_GROUPED)) %>%
  distinct(.) %>%
  group_by(EUNIS_LEVEL_2_GROUPED) %>%
  summarise(N_Inventories = n()) %>%
  left_join(SIVIM_rare %>%
              dplyr::select(c(TAXON_REF_PYR, EUNIS_LEVEL_2_GROUPED)) %>%
              distinct(.) %>%
              group_by(EUNIS_LEVEL_2_GROUPED) %>%
              summarise(Richness = n()), 'EUNIS_LEVEL_2_GROUPED') %>%
  left_join(SIVIM_rare %>%
              dplyr::select(c(UTM10, EUNIS_LEVEL_2_GROUPED)) %>%
              distinct(.) %>%
              group_by(EUNIS_LEVEL_2_GROUPED) %>%
              summarise(N_UTM10 = n()),
            'EUNIS_LEVEL_2_GROUPED') %>%
  rename(EUNIS = EUNIS_LEVEL_2_GROUPED) %>%
  write.xlsx(file = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Table_1_Habitats.xlsx')

save(Rareza, Frecuencia_por_habitat, file = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData')
save(Rareza, Frecuencia_por_habitat, file = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData')
