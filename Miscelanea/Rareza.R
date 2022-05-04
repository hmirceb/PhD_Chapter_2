load('~/Dropbox/Hector_Bego/Frontiers Ecol Evol/Data_analysis/Scripts/DatosMicroR.RData')
load('~/Dropbox/Tesis/PD/Resultados/PD_Pirineo.RData')
load("/home/hector/Dropbox/Hector_Bego/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")

library(dplyr)
library(emmeans)
library(visreg)
library(DHARMa)
library(lme4)
library(reshape2)
library(readxl)
library(ggplot2)
library(PhyloMeasures)
library(popbio)
library(plotly)
library(scales)
library(ggtree)

rm(SIVIM_abund_invXsp_angios, limpios, limpios.angios, tre, tab_v2)
    
    SIVIM_rare<-SIVIM; rm(SIVIM)
    SIVIM_rare<-SIVIM_rare[SIVIM_rare$EUNIS_LEVEL_2_GROUPED != '',] #Quitamos habitats antropicos y sin definir
    SIVIM_rare<-SIVIM_rare[SIVIM_rare$NATIVE_PYR != '',] # Quitamos aloctonas
    SIVIM_rare<-SIVIM_rare[!endsWith(SIVIM_rare$TAXON_REF_PYR_MOD, 'sp.'),] # Quitamos especies en duda
    
    # Contamos en cuantas UTM 10 aparece cada planta y creamos un df con esos datos
    a<-aggregate(UTM10~TAXON_REF_PYR_MOD, SIVIM_rare, table)
    range<-data.frame('sp' = a$TAXON_REF_PYR, 'rango' = apply(a, 1, function(x) length(x[[2]])))
    
    # Lo mismo pero para los habitats
    a<-aggregate(EUNIS_LEVEL_2_GROUPED~TAXON_REF_PYR_MOD, SIVIM_rare, table)
    range$habspe<-apply(a, 1, function(x) length(x[[2]]))
    
    # DF con las especies (solo especies, no subespecies)
    end<-as.data.frame(unique(cbind(SIVIM_rare$TAXON_REF_PYR_MOD, SIVIM_rare$ENDEMIC_PYR)))
    end$end<-ifelse(end$V2 == "Endemic_sp", 'End', 'NoEnd')
    names(end)[1]<-'sp'
    
    # Comprobamos que no haya ninguna especie que aparezca en mas inventarios que 
      # habitats o UTMs
    raras<-as.data.frame(table(SIVIM_rare$TAXON_REF_PYR_MOD))
    names(raras)[1]<-'sp'
    utms<-left_join(raras, range, by = 'sp')
    utms[utms$Freq < utms$rango,]
    utms[utms$Freq < utms$habspe,]
    if(nrow(utms[utms$Freq < utms$rango,]) == nrow(utms[utms$Freq < utms$habspe,]) & 
       nrow(utms[utms$Freq < utms$rango,]) == 0){rm(utms, raras, a)}
    
    # Unimos endemicas, rango y especificidad
    range<-merge(range, end, by = 'sp')
    range<-range[,-which(names(range) == 'V2')]
    range$rango_p<-100*(range$rango/max(range$rango))
   
    # Efecto de distintos umbrales en la captura de especies raras:
    np<-list()
    for (i in 1:100) {
      np[[i]]<-ifelse(range$rango_p <= i, 'Narrow', 'Wide')
    }
    np<-as.data.frame(do.call('rbind', lapply(np, FUN = function(x) table(x))))
    np[100, 2]<-0
    np$narrow_p<-np$Narrow/(np$Narrow+np$Wide)
    ep<-list()
    for (i in 1:length(unique(SIVIM_rare$EUNIS_LEVEL_2_GROUPED))) {
      ep[[i]]<-as.factor(ifelse(range$habspe <= i &
                                  range$rango > 1 , 'Esp', 'NoEsp'))
    }
    ep2<-list()
    for (i in 1:length(unique(SIVIM_rare$EUNIS_LEVEL_2_GROUPED))) {
      ep2[[i]]<-as.factor(ifelse(range$habspe <= i, 'Esp', 'NoEsp'))
    }
    ep<-as.data.frame(do.call('rbind', lapply(ep, FUN = function(x) table(x))))
    ep2<-as.data.frame(do.call('rbind', lapply(ep2, FUN = function(x) table(x))))
    ep2$NoEsp[15]<-0
    ep$esp_p<-ep$Esp/(ep$Esp+ep$NoEsp)
    ep2$esp_p<-ep2$Esp/(ep2$Esp+ep2$NoEsp)
    ep$a<-'Specialist > 1UTM'; ep$n<-100*(1:15/15); names(ep)<-c('Rare', 'Common', 'Prop', 'a', 'n')
    ep2$a<-'Specialist All'; ep2$n<-100*(1:15/15); names(ep2)<-c('Rare', 'Common', 'Prop', 'a', 'n')
    np$a<-'Narrow'; np$n<-1:100; names(np)<-c('Rare', 'Common', 'Prop', 'a', 'n')
    ggplot(as.data.frame(rbind(ep, ep2, np)) , aes(x = n, y = Prop, color = a))+geom_line(size = 1)+
      scale_x_continuous(breaks = seq(0, 100, 5))+
      geom_vline(aes(xintercept = 5))+geom_vline(aes(xintercept = min(ep$n)))
    # Si para asegurarnos de que las especialistas
          # estan bien no consideramos aquellas que solo esten en 1 UTM, entonces
          # hay 269 especies que quedan fuera y por eso la grafica no llega al 100%
    
    
    # Umbrales para considerar los rangos amplio o estrechos y especificos o no
      #Probamos con 10 y 5 % del rango maximo observado y con especies en 1 habitat
    range$rango_cat10<-as.factor(ifelse(range$rango_p > 10, 'Wide', 'Narrow'))
    range$rango_cat5<-as.factor(ifelse(range$rango_p > 5, 'Wide', 'Narrow'))
    
    # Especialistas si estan en un solo habitat y mas de una utm
    range$habspe_cat<-as.factor(ifelse(range$habspe == 1 &
                                         range$rango >1 , 'Esp', 'NoEsp'))
    range$Rare10<-as.factor(paste(range$rango_cat10, range$habspe_cat, range$end))
    range$Rare5<-as.factor(paste(range$rango_cat5, range$habspe_cat, range$end))
    
    # Determinamos las especies raras
    range$Rarity5<-as.factor(ifelse(range$Rare5 == 'Wide NoEsp NoEnd', 'Common', 'Rare'))
    range$Rarity10<-as.factor(ifelse(range$Rare10 == 'Wide NoEsp NoEnd', 'Common', 'Rare'))
    range$end<-as.factor(range$end)
    
    table(range$Rare5)
    table(range$Rare10)
    100*sum(table(range$Rare5)[1:6])/sum(table(range$Rare5)) ## 70% de las especies son raras
    100*sum(table(range$Rare5)[1:4])/sum(table(range$Rare5)[1:6]) ## 99% de las especies raras tienen narrow distribution
    table(range$Rarity5)/sum(table(range$Rarity5))
    table(range$Rarity10)/sum(table(range$Rarity10))
    
    # Ponemos la rareza de cada especie en SIVIM
    names(range)[1]<-'TAXON_REF_PYR_MOD'
    range<-range[!endsWith(range$TAXON_REF_PYR_MOD, 'sp.'),]
    SIVIM_rare<-left_join(SIVIM_rare, range, by = 'TAXON_REF_PYR_MOD')
    SIVIM_rare<-SIVIM_rare[SIVIM_rare$EUNIS_LEVEL_2_GROUPED != '',] # Quitamos inventarios sin habitat
    
    # Conteo de cada tipo de rareza por inventario
    b<-merge(merge(merge(dcast(SIVIM_rare, INV_CODE~rango_cat5, length), 
                   dcast(SIVIM_rare, INV_CODE~rango_cat10, length), by = 'INV_CODE'),
                   dcast(SIVIM_rare, INV_CODE~habspe_cat, length), by = 'INV_CODE'),
             dcast(SIVIM_rare, INV_CODE~end, length) , by = 'INV_CODE')
    names(b)[endsWith(names(b), '.x')]<-c('Narrow_5', 'Wide_5')
    names(b)[endsWith(names(b), '.y')]<-c('Narrow_10', 'Wide_10')
    
    # Le ponemos el habitat
    b<-left_join(b, unique(data.frame('INV_CODE' = SIVIM_rare$INV_CODE, 'Hab' = SIVIM_rare$EUNIS_LEVEL_2_GROUPED)), 'INV_CODE')
    
    # Contamos las raras por inventario y lo juntamos al resto
    z<-merge(dcast(SIVIM_rare, INV_CODE~Rarity5, length),
             dcast(SIVIM_rare, INV_CODE~Rarity5, length), by = 'INV_CODE')
    names(z)[endsWith(names(z), '.y')]<-c('Common_10', 'Rare_10')
    names(z)[endsWith(names(z), '.x')]<-c('Common_5', 'Rare_5')
    
    b<-merge(b, z, by = 'INV_CODE')
    sum(is.na(b))
    b<-right_join(b, PD_Pirineo, by = 'INV_CODE')
    
    # Reclasificamos la habitats para que salgan abreviados
    oldvalues <- unique(TablaCodigos$Descrip_Grouped)
    newvalues <- factor(unique(TablaCodigos$Code_Grouped))
    b$EUNIS_LEVEL_2_GROUPED_CODE <- newvalues[ match(b$Hab, oldvalues) ]
    
    b<-b[complete.cases(b),] # Nos quedamos con los inventarios que tengan TODOS los datos
    # si quitase las variable ambientales saldrian mas
  rm(end, z, oldvalues, newvalues)
  
  # Plot 3D de como influye cada tipo de rareza en los inventarios
  plot_ly(data = data.frame('Rango_restringido' = (b$Narrow_5/(b$Narrow_5+b$Wide_5)), 
                            'Especialización' = (b$Esp/(b$Esp+b$NoEsp)),
                            'Endemismo' = (b$End/(b$End+b$NoEnd)),
                            'hab' = (b$EUNIS_LEVEL_2_GROUPED_CODE)), 
          x = ~Rango_restringido, y = ~Especialización, z = ~Endemismo, color = ~hab,
          colors = hue_pal()(15), size = 5, alpha = 0.8)
  
  
  # Calculamos la proporcion de PD correspondiente a las especies raras
  matriz_raras<-ifelse(SIVIM_abund_invXsp[, names(SIVIM_abund_invXsp) 
                                          %in% range[range$Rarity10 == 'Rare',]$TAXON_REF_PYR_MOD] >0,
                       1, 0)
  pd_raras<-mean.list(lapply(trees_sp_yule_10, FUN = function(x) as.matrix(pd.query(x, matriz_raras, standardize = F))))
  
  b<-left_join(b, data.frame('INV_CODE' = SIVIM_abund_invXsp$INV_CODE, pd_raras), 'INV_CODE')
  b$pd_common<-b$pd.obs-b$pd_raras
  b$pd_prop<-(b$pd_raras/b$pd.obs)
  rm(pd_raras, matriz_raras)
  
  ### El OLRE (observation level random effect) es para quitar la sobredispersion.
    # En algunos no hace falta usarla y se puede usar un GLM normal
    # En las que no se puede usar es porque hay mucho valores muy bajos y no 
    # converge el algoritmo
  
  mod<-glmer(cbind(Narrow_5, Wide_5)~EUNIS_LEVEL_2_GROUPED_CODE + (1|INV_CODE), 
             family = binomial(link = logit), data = b, control = glmerControl(optimizer = 'bobyqa'))
  summary(mod)
  emmeans(mod, specs = 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response')
  visreg(mod, 'EUNIS_LEVEL_2_GROUPED_CODE', scale = 'response')
  plot(summary(emmeans(mod, specs = 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response')))
  
  mod1<-glm(cbind(Esp, NoEsp)~EUNIS_LEVEL_2_GROUPED_CODE, 
              family = binomial(link = logit), data = b)
  summary(mod1)
  emmeans(mod1, specs = 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response')
  visreg(mod1, 'EUNIS_LEVEL_2_GROUPED_CODE', scale = 'response')
  
  mod2<-glm(cbind(End, NoEnd)~EUNIS_LEVEL_2_GROUPED_CODE, 
              family = binomial(link = logit), data = b)
  summary(mod2)
  emmeans(mod2, specs = 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response')
  visreg(mod2, 'EUNIS_LEVEL_2_GROUPED_CODE', scale = 'response')
  
  mod3<-glmer(cbind(Rare_5, Common_5)~EUNIS_LEVEL_2_GROUPED_CODE + (1|INV_CODE), 
              family = binomial(link = logit), data = b, control = glmerControl(optimizer = 'bobyqa'))
  summary(mod3)
  emmeans(mod3, specs = 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response')
  visreg(mod3, 'EUNIS_LEVEL_2_GROUPED_CODE', scale = 'response')
  car::Anova(mod3)
  plot(summary(emmeans(mod3, specs = 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response')))
  
  mod5<-glmer(cbind(Rare_5, Common_5)~pd.obs.z + (1|EUNIS_LEVEL_2_GROUPED_CODE), 
              family = binomial(link = logit), data = b, control = glmerControl(optimizer = 'bobyqa'))
  summary(mod5)
  emmeans(mod5, specs = 'pd.obs.z', type = 'response')
  visreg(mod5, 'pd.obs.z', scale = 'response')

  mod6<-glm(cbind(pd_raras, pd_common)~EUNIS_LEVEL_2_GROUPED_CODE, 
              family = binomial(link = logit), data = b)
  summary(mod6)
  emmeans(mod6, specs = 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response')
  visreg(mod6, 'EUNIS_LEVEL_2_GROUPED_CODE', scale = 'response')
  
### 
simulateResiduals(mod, plot = T)
simulateResiduals(mod1, plot = T)
simulateResiduals(mod2, plot = T)
simulateResiduals(mod3, plot = T)
simulateResiduals(mod4, plot = T)
simulateResiduals(mod5, plot = T)

###
save(mod, mod1, mod2, mod3, mod4, mod5, b, file = 'ModelosRareza.RData')

aggregate(b$Rare_5/b$ntaxa, by = list(b$EUNIS_LEVEL_2_GROUPED_CODE), mean)
aggregate(b$Rare_5/b$ntaxa, by = list(b$EUNIS_LEVEL_2_GROUPED_CODE), sd)
aggregate(b$Narrow_5/b$ntaxa, by = list(b$EUNIS_LEVEL_2_GROUPED_CODE), mean)
aggregate(b$Narrow_5/b$ntaxa, by = list(b$EUNIS_LEVEL_2_GROUPED_CODE), sd)
aggregate(b$Esp/b$ntaxa, by = list(b$EUNIS_LEVEL_2_GROUPED_CODE), mean)
aggregate(b$Esp/b$ntaxa, by = list(b$EUNIS_LEVEL_2_GROUPED_CODE), sd)
aggregate(b$End/b$ntaxa, by = list(b$EUNIS_LEVEL_2_GROUPED_CODE), mean)
aggregate(b$End/b$ntaxa, by = list(b$EUNIS_LEVEL_2_GROUPED_CODE), sd)

raras<-unique(data.frame('EUNIS_LEVEL_2_GROUPED' = SIVIM_rare$EUNIS_LEVEL_2_GROUPED, 
                         'Rare' = SIVIM_rare$Rare5, 'sp' = SIVIM_rare$TAXON_REF_PYR_MOD, 
                         'Rarity' = SIVIM_rare$Rarity5))
conteo<-(merge(dcast(raras, EUNIS_LEVEL_2_GROUPED~Rare, length),
      dcast(raras, EUNIS_LEVEL_2_GROUPED~Rarity, length)))
conteo$prop_rare<-conteo$Rare/(conteo$Rare+conteo$Common)
View(conteo)
mean(conteo$prop_rare)
sd(conteo$prop_rare)
b$prop_rare<-100*b$Rare_5/b$ntaxa
aggregate(prop_rare~EUNIS_LEVEL_2_GROUPED_CODE, b, mean)

boxplot(prop_rare~EUNIS_LEVEL_2_GROUPED_CODE, b)
pairs(emmeans(mod3, specs = 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response'))
summary(mod1)$deviance/(summary(mod1)$null.deviance+summary(mod1)$deviance)
a<-summary(emmeans(mod3, specs = 'EUNIS_LEVEL_2_GROUPED_CODE', type = 'response'))
ggplot(a, aes(x = EUNIS_LEVEL_2_GROUPED_CODE, y = prob, color = EUNIS_LEVEL_2_GROUPED_CODE))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))+theme_bw()

plot(x = b$pd.obs.z, y = b$Rare_5/b$Common_5)
plot(x = b$ntaxa, y = b$Rare_5/(b$Rare_5+b$Common_5))

# Para ver que funcion de enlace usar y si conviene meter un termino cuadratico
lin<-glm(cbind(Rare_5, Common_5)~pd.obs.z, 
          family = binomial(link = logit), data = b)
quad<-glm(cbind(Rare_5, Common_5)~pd.obs.z+I(pd.obs.z^2), 
          family = binomial(link = logit), data = b)
logi<-glm(cbind(Rare_5, Common_5)~pd.obs.z+I(pd.obs.z^2), 
            family = binomial(link = logit), data = b)
probi<-glm(cbind(Rare_5, Common_5)~pd.obs.z+I(pd.obs.z^2), 
          family = binomial(link = probit), data = b)
clog<-glm(cbind(Rare_5, Common_5)~pd.obs.z+I(pd.obs.z^2),
          family = binomial(link = cloglog), data = b)
visreg(logi, scale = 'response')
visreg(probi, scale = 'response')
visreg(clog, scale = 'response')
AIC(lin, quad)
AIC(logi, probi, clog)

#### Arboles por rareza
arbol<-keep.tip(trees_sp_yule_10[[1]], range$TAXON_REF_PYR_MOD)
arbol<-groupOTU(arbol, arbol$tip.label %in% range[range$habspe_cat == 'Esp',]$TAXON_REF_PYR_MOD, 
                group_name = 'Especialista')
arbol<-groupOTU(arbol, arbol$tip.label %in% range[range$rango_cat5 == 'Narrow',]$TAXON_REF_PYR_MOD, 
                group_name = 'Rango restringido')
arbol<-groupOTU(arbol, arbol$tip.label %in% range[range$end == 'End',]$TAXON_REF_PYR_MOD, 
                group_name = 'Endemico')

arboles<-ggarrange(ggtree(arbol, layout = 'circular')+geom_tree(aes(color = Endemico))+
                     scale_color_manual(values = c('grey', 'red'), labels = c('NoEnd', 'End'))+
                     theme(legend.position = 'bottom',
                           text = element_text(size = 20)),
                   ggtree(arbol, layout = 'circular')+geom_tree(aes(color = `Rango restringido`))+
                     scale_color_manual(values = c('grey', 'blue'), labels = c('Wide', 'Narrow'))+
                     theme(legend.position = 'bottom',
                           text = element_text(size = 20)),
                   ggtree(arbol, layout = 'circular')+geom_tree(aes(color = Especialista))+
                     scale_color_manual(values = c('grey', 'green'), labels = c('NoEsp', 'Esp'))+
                     theme(legend.position = 'bottom',
                           text = element_text(size = 20)),
                   nrow = 1)
ggsave(arboles, width = 2*29.7, height = 2*21, units = 'cm', filename = 'Arboles.jpg')

##########################
##### Balloon plot ######
#########################
library(ggridges)
data<-data.frame('Rango_restringido' = (b$Narrow_5/(b$Narrow_5+b$Wide_5)), 
                 'Especialización' = (b$Esp/(b$Esp+b$NoEsp)),
                 'Endemismo' = (b$End/(b$End+b$NoEnd)),
                 'hab' = (b$EUNIS_LEVEL_2_GROUPED_CODE))

grid.arrange(ggplot(data, aes(x = Rango_restringido, y = hab))+
               geom_density_ridges(fill = 'red', alpha = 0.5)+scale_x_continuous(limits = c(0, 0.25)),
             ggplot(data, aes(x = Especialización, y = hab))+
               geom_density_ridges(fill = 'green', alpha = 0.5)+scale_x_continuous(limits = c(0, 0.25)),
             ggplot(data, aes(x = Endemismo, y = hab))+
               geom_density_ridges(fill = 'blue', alpha = 0.5)+scale_x_continuous(limits = c(0, 0.25)))

a<-merge(merge(aggregate(Rango_restringido~hab, data, mean),
               aggregate(Especialización~hab, data, mean), by = 'hab'),
         aggregate(Endemismo~hab, data, mean), by = 'hab')
b<-merge(merge(aggregate(Rango_restringido~hab, data, sd),
               aggregate(Especialización~hab, data, sd), by = 'hab'),
         aggregate(Endemismo~hab, data, sd), by = 'hab')

p<-gather(a, rareza, valor, Rango_restringido:Endemismo)
p2<-gather(b, rareza2, valor2, Rango_restringido:Endemismo)

ggplot(p)+geom_point(aes(x = hab, y = rareza2, size = valor2), data = p2)+
  geom_point(aes(x = hab, y = rareza, size = valor, fill = hab), shape = 21)+theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = 'none')+
  scale_size(range = c(5,20))
