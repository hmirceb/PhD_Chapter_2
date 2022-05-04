library(picante)
library(dplyr)
library(tidytree)
library(reshape2)
library(emmeans)
library(visreg)
library(DHARMa)
library(lme4)

load("/home/hector/Dropbox/Hector_Bego/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
load('~/Dropbox/Hector_Bego/Paper_microrefugia_SI/Data_analysis/Scripts/DatosMicroR.RData')
load('~/Dropbox/Tesis/PD/Resultados/PD_Pirineo.RData')
rm(tre, tab_v2, SIVIM_abund_invXsp, SIVIM_abund_invXsp_angios)

# Sacamos los nodos correspondientes a cada especie
ramas<-as_tibble(trees_sp_yule_10[[1]])
l<-list()
for (i in 1:dim(ramas[!is.na(ramas$label),])[1]) {
  l[[i]]<-ancestor(ramas, ramas[[i,2]])
}
# Y les asignamos el nombre de su especie
for (j in 1:length(l)) {
  l[[j]]<-cbind(l[[j]], ramas[[j,4]])
}

# Las juntamos y filtramos las especies que no estan en los inventarios
prueba<-do.call('rbind', l)
names(prueba)[5]<-'TAXON_REF_PYR_MOD'
cosa<-prueba[prueba$TAXON_REF_PYR_MOD %in% SIVIM$TAXON_REF_PYR_MOD,]

# Cogemos las UTM10, habitas y codigos de los inventarios y lo juntamos con los nodos
sps_data<-data.frame('TAXON_REF_PYR_MOD' = SIVIM$TAXON_REF_PYR_MOD, 
           'UTM10' = SIVIM$UTM10, 'EUNIS_LEVEL_2_GROUPED' = SIVIM$EUNIS_LEVEL_2_GROUPED, 'INV_CODE' = SIVIM$INV_CODE)
locura<-full_join(cosa, sps_data, by = 'TAXON_REF_PYR_MOD')
locura$node<-as.factor(locura$node)

# Calculamos en cuantas UMT10 y habitats aparece cada nodo
phylo_utm<-aggregate(UTM10~node, data = locura, function(x)length(unique(x)))
names(phylo_utm)[2]<-'N_UTM10'
phylo_hab<-aggregate(EUNIS_LEVEL_2_GROUPED~node, data = locura, function(x)length(unique(x)))
names(phylo_hab)[2]<-'N_Habs'
nos<-full_join(phylo_hab, phylo_utm, by = 'node')

# Clasificamos el rango y la especificidad
nos$N_Habs_cat<-ifelse(nos$N_Habs/max(nos$N_Habs) > 0.1, 'NoEsp_Phylo', 'Esp_Phylo')
nos$N_UTM10_cat<-ifelse(nos$N_UTM10/max(nos$N_UTM10) > 0.1, 'Wide_Phylo', 'Narrow_Phylo')
locura<-full_join(locura, nos, by = 'node')
locura<-locura[locura$EUNIS_LEVEL_2_GROUPED != '',] # Quitamos inventarios sin habitats
locura$Phylo_Rarity<-paste(locura$N_Habs_cat, locura$N_UTM10_cat)
locura$Phylo_Rarity<-ifelse(locura$Phylo_Rarity == 'NoEsp_Phylo Wide_Phylo', 'Common_Phylo', 'Rare_Phylo')

t<-merge(merge(dcast(locura, INV_CODE~N_UTM10_cat, length), 
               dcast(locura, INV_CODE~N_Habs_cat, length), by = 'INV_CODE'),
         dcast(locura, INV_CODE~Phylo_Rarity, length), by = 'INV_CODE')
t<-right_join(t, limpios, by = 'INV_CODE')
t<-right_join(t, PD_Pirineo, by = 'INV_CODE')
t<-t[complete.cases(t),]

rm(nos, phylo_hab, phylo_utm, ramas, sps_data, l, cosa, prueba, i, j)

mod_phylo<-glmer(cbind(Rare_Phylo, Common_Phylo)~EUNIS_LEVEL_2_GROUPED_CODE+(1|INV_CODE), 
           family = binomial(link = logit), data = t, control = glmerControl(optimizer = 'bobyqa'))
summary(mod_phylo)
visreg(mod_phylo, 'EUNIS_LEVEL_2_GROUPED_CODE', scale = 'response')
emmeans(mod_phylo, 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response')

mod1_phylo<-glmer(cbind(Narrow_Phylo, Wide_Phylo)~EUNIS_LEVEL_2_GROUPED_CODE+(1|INV_CODE), 
            family = binomial(link = logit), data = t, control = glmerControl(optimizer = 'bobyqa'))
summary(mod1_phylo)
visreg(mod1_phylo, 'EUNIS_LEVEL_2_GROUPED_CODE', scale = 'response')
emmeans(mod1_phylo, 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response')

mod2_phylo<-glm(cbind(Esp_Phylo, NoEsp_Phylo)~EUNIS_LEVEL_2_GROUPED_CODE, family = binomial(link = logit), data = t)
summary(mod2_phylo)
visreg(mod2_phylo, 'EUNIS_LEVEL_2_GROUPED_CODE', scale = 'response')
emmeans(mod2_phylo, 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response')

simulateResiduals(mod_phylo, plot = T)
simulateResiduals(mod1_phylo, plot = T)
simulateResiduals(mod2_phylo, plot = T)

save(mod_phylo, mod1_phylo, mod2_phylo, t, file = 'Modelos_Rareza_Filo.RData')

aggregate(t$Rare_Phylo/(t$Rare_Phylo+t$Common_Phylo), by = list(t$EUNIS_LEVEL_2_GROUPED_CODE), mean)
aggregate(t$Rare_Phylo/(t$Rare_Phylo+t$Common_Phylo), by = list(t$EUNIS_LEVEL_2_GROUPED_CODE), sd)
aggregate(t$Narrow_Phylo/(t$Narrow_Phylo+t$Wide_Phylo), by = list(t$EUNIS_LEVEL_2_GROUPED_CODE), mean)
aggregate(t$Narrow_Phylo/(t$Narrow_Phylo+t$Wide_Phylo), by = list(t$EUNIS_LEVEL_2_GROUPED_CODE), sd)
aggregate(t$Esp_Phylo/(t$Esp_Phylo+t$NoEsp_Phylo), by = list(t$EUNIS_LEVEL_2_GROUPED_CODE), mean)
aggregate(t$Esp_Phylo/(t$Esp_Phylo+t$NoEsp_Phylo), by = list(t$EUNIS_LEVEL_2_GROUPED_CODE), sd)

raras_phylo<-unique(data.frame('EUNIS_LEVEL_2_GROUPED' = locura$EUNIS_LEVEL_2_GROUPED, 
                         'Rare' = locura$Phylo_Rarity, 'sp' = locura$TAXON_REF_PYR_MOD))
conteo_phylo<-dcast(raras_phylo, EUNIS_LEVEL_2_GROUPED~Rare, length)
conteo_phylo$prop_rare<-conteo_phylo$Rare/(conteo_phylo$Rare+conteo_phylo$Common)