library(treedataverse)
library(tidyverse)
library(phylosignal)
library(phylobase)
library(picante)
library(ggnewscale)
library(ggtreeExtra)
library(popbio)
library(phytools)
library(ggpubr)
library(ggridges)
library(lme4)
library(rstanarm)
library(insight)
library(gridExtra)

load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/FiloSignal_Rareza_scaled_All.RData")
load('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/FiloSignal_Rareza_scaled_All.RData')
load("C:/Users/18172844S/Dropbox/DATA__LAB/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
load("~/Dropbox/DATA__LAB/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
load("C:/Users/18172844S/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")
load("~/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")
rm(tab_v2, tre, w1)

data<-eiv %>% 
  dplyr::select(c(TAXON_REF_PYR_MOD, RA, RA_max, LA, n_habs, B, B_h, PS, Endemic, CAT_AMENAZA_PYR)) %>%
  mutate(RA_max = as.numeric(RA_max)) %>%
  mutate(Endemic = ifelse(Endemic == 'End', 1, 0)) %>%
  mutate(Threat = ifelse(CAT_AMENAZA_PYR %in% c('CR', 'VU', 'EN'), 'Thrt', 'NoThrt')) %>%
  distinct()


lipa_all_df<-as.data.frame(mean.list(lapply(lipa_all, function(x)
  cbind(x$lipa, x$p.value))))
names(lipa_all_df)<-paste('LIPA', names(lipa_all_df), sep = '_')
names(lipa_all_df)[5:8]<-paste(names(lipa_all_df)[5:8], 'pval', sep = '_')
lipa_all_df$TAXON_REF_PYR_MOD<-rownames(lipa_all_df)

data_all<-inner_join(data, lipa_all_df, 'TAXON_REF_PYR_MOD') %>%
  mutate(Endemic = as.factor(Endemic))

tree_figures<-lapply(trees_sp_yule_10, function(x) keep.tip(x, data_all$TAXON_REF_PYR_MOD))

runs<-999

####
eliminar_thrt<-data_all %>%
  filter(Threat == 'Thrt') %>%
  pull(TAXON_REF_PYR_MOD)

obs_total<-mean(sapply(tree_figures, function(x) sum(x$edge.length)))
obs_loss_thrt<-mean(sapply(tree_figures, function(x) (sum(x$edge.length)-
                                                        sum(drop.tip(x, eliminar_thrt)$edge.length))/sum(x$edge.length)))

rand_loss_thrt<-list()
for (i in 1:runs) {
  rand_loss_thrt[[i]]<-mean(sapply(tree_figures, function(x) (sum(x$edge.length)-
                                                                sum(drop.tip(x, sample(data_all$TAXON_REF_PYR_MOD, 
                                                                                       size = length(eliminar_thrt)))$edge.length))/sum(x$edge.length)))
}

###
eliminar_end<-data_all %>%
  filter(LIPA_Endemic > 0 & 
            LIPA_Endemic_pval <= 0.05 & Endemic == '1') %>%
  pull(TAXON_REF_PYR_MOD)

obs_total<-mean(sapply(tree_figures, function(x) sum(x$edge.length)))
obs_loss_end<-mean(sapply(tree_figures, function(x) (sum(x$edge.length)-
                                                       sum(drop.tip(x, eliminar_end)$edge.length))/sum(x$edge.length)))

rand_loss_end<-list()
for (i in 1:runs) {
  rand_loss_end[[i]]<-mean(sapply(tree_figures, function(x) (sum(x$edge.length)-
                                                               sum(drop.tip(x, sample(data_all$TAXON_REF_PYR_MOD, 
                                                                                      size = length(eliminar_end)))$edge.length))/sum(x$edge.length)))
}

###
eliminar_ra<-data_all %>%
  filter(LIPA_RA_max > 0 & 
            LIPA_RA_max_pval <= 0.05 & scale(RA_max) < 0) %>%
  pull(TAXON_REF_PYR_MOD)

obs_loss_ra<-mean(sapply(tree_figures, function(x) (sum(x$edge.length)-
                                                      sum(drop.tip(x, eliminar_ra)$edge.length))/sum(x$edge.length)))

rand_loss_ra<-list()
for (i in 1:runs) {
  rand_loss_ra[[i]]<-mean(sapply(tree_figures, function(x) (sum(x$edge.length)-
                                                              sum(drop.tip(x, sample(data_all$TAXON_REF_PYR_MOD, 
                                                                                     size = length(eliminar_ra)))$edge.length))/sum(x$edge.length)))
}
rank(c(obs_loss_ra, unlist(rand_loss_ra)))[1]

###
eliminar_hs<-data_all %>%
  filter(LIPA_B_h > 0 & 
            LIPA_B_h_pval <= 0.05 & scale(B_h) < 0) %>%
  pull(TAXON_REF_PYR_MOD)

obs_loss_hs<-mean(sapply(tree_figures, function(x) (sum(x$edge.length)-
                                                      sum(drop.tip(x, eliminar_hs)$edge.length))/sum(x$edge.length)))

rand_loss_hs<-list()
for (i in 1:runs) {
  rand_loss_hs[[i]]<-mean(sapply(tree_figures, function(x) (sum(x$edge.length)-
                                                              sum(drop.tip(x, sample(data_all$TAXON_REF_PYR_MOD, 
                                                                                     size = length(eliminar_hs)))$edge.length))/sum(x$edge.length)))
}
rank(c(obs_loss_hs, unlist(rand_loss_hs)))[1]

###
eliminar_la<-data_all %>%
  filter(LIPA_LA > 0 & LIPA_LA_pval <= 0.05 & scale(LA) < 0) %>%
  pull(TAXON_REF_PYR_MOD)

obs_loss_la<-mean(sapply(tree_figures, function(x) (sum(x$edge.length)-
                                                      sum(drop.tip(x, eliminar_la)$edge.length))/sum(x$edge.length)))

rand_loss_la<-list()
for (i in 1:runs) {
  rand_loss_la[[i]]<-mean(sapply(tree_figures, function(x) (sum(x$edge.length)-
                                                              sum(drop.tip(x, sample(data_all$TAXON_REF_PYR_MOD, 
                                                                                     size = length(eliminar_la)))$edge.length))/sum(x$edge.length)))
}

save(obs_loss_end, obs_loss_hs, obs_loss_la, obs_loss_ra, 
     rand_loss_end, rand_loss_hs, rand_loss_la, rand_loss_ra,
     rand_loss_thrt, obs_loss_thrt, 
     obs_total,
     file = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/PD_loss_all.RData')

load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/PD_loss_all.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/PD_loss_all.RData")

q<-data.frame('Value' = c(unlist(rand_loss_end), 
                          unlist(rand_loss_ra),
                          unlist(rand_loss_hs),
                          unlist(rand_loss_la)),
              'Rarity' = rep(c('Endemic',
                               'Regional geographic range',
                               'Habitat specialization',
                               'Local abundance'), each = 999)) %>%
  mutate(Rarity = factor(Rarity, levels = c( 
    'Endemic',
    'Regional geographic range',
    'Habitat specialization',
    'Local abundance'))) %>%
  group_by(Rarity) %>%
  summarise(quantile(Value, c(0.025, 0.975))) %>%
  ungroup() %>%
  mutate(Interval = rep(c('Low', 'Up'), times = 4)) %>%
  rename(q =`quantile(Value, c(0.025, 0.975))` )

datos<-data.frame('Value' = c(unlist(rand_loss_end), 
                              unlist(rand_loss_ra),
                              unlist(rand_loss_hs),
                              unlist(rand_loss_la),
                              unlist(rand_loss_thrt)),
                  'Rarity' = rep(c('Endemic',
                                   'Regional geographic range',
                                   'Habitat specialization',
                                   'Local abundance',
                                   'Red List'), each = 999)) %>% 
  mutate(Rarity = factor(Rarity, levels = c('Red List', 
                                            'Endemic',
                                            'Regional geographic range',
    'Habitat specialization',
    'Local abundance')))


aaa<-data.frame(
  'Line' = c(obs_loss_end, obs_loss_ra, obs_loss_hs, obs_loss_la, obs_loss_thrt),
  'Mean' = c(mean(unlist(rand_loss_end)), 
             mean(unlist(rand_loss_ra)),
             mean(unlist(rand_loss_hs)),
             mean(unlist(rand_loss_la)),
             mean(unlist(rand_loss_thrt))),
  'SD' = c(sd(unlist(rand_loss_end)), 
           sd(unlist(rand_loss_ra)),
           sd(unlist(rand_loss_hs)),
           sd(unlist(rand_loss_la)),
           sd(unlist(rand_loss_thrt))),
  'Rarity' = factor(c('Endemic', 'RGR', 'HS', 'LA', 'Red List')))
aaa$ses<-(aaa$Line-aaa$Mean)/aaa$SD

ses_loss_plot<-aaa %>% 
  mutate(Rarity = factor(Rarity, levels = c('LA', 'HS', 'RGR', 'Endemic', 'Red List'))) %>%
  ggplot(aes(y = Rarity))+
  geom_segment(aes(xend = ses, x = 0, yend = Rarity),
               linetype = 'dashed', size = 1)+
  theme_bw()+
  geom_point(aes(x = ses, color = Rarity), size = 5)+
  scale_color_manual(values = c('Red List' = 'mediumorchid',
                                'Endemic' = '#FF8C00', 
                                'RGR' = 'firebrick',
                                'HS' = 'steelblue',
                                'LA' = 'forestgreen'), 
                     name = "Rarity type:", guide = 'none')+
  theme(text = element_text(size = 20))+
  xlab('Standard effect size of PD loss')+
  ylab('Rarity type')+
  annotate(geom = 'rect', xmin = -1.96, 
           xmax = 1.96, ymin = 0, ymax = 6,
           alpha = 0.1)+
  ggtitle('B')

eiv$Endemic<-factor(eiv$Endemic, levels = c('NoEnd', 'End'))
eiv<-eiv %>% mutate(Threat = ifelse(CAT_AMENAZA_PYR %in% c('CR', 'VU', 'EN'), 'Thrt', 'NoThrt'))
eiv$Threat<-factor(eiv$Threat, levels = c('NoThrt', 'Thrt'))

m<-stan_glmer(Endemic~1+(1|Order/FAMILIA_WFO/GENERO_RP), family = binomial, data = eiv,
              iter = 4000)

m0<-stan_glmer(Threat~1+(1|Order/FAMILIA_WFO/GENERO_RP), family = binomial, data = eiv,
               iter = 4000)

m1<-stan_glmer(scale(as.numeric(RA_max))~1+(1|Order/FAMILIA_WFO/GENERO_RP), data = eiv, 
               iter = 4000)

m2<-stan_glmer(scale(B_h)~1+(1|Order/FAMILIA_WFO/GENERO_RP), data = eiv,
               iter = 4000)

m3<-stan_glmer(scale(LA)~1+(1|Order/FAMILIA_WFO/GENERO_RP), data = eiv,
               iter = 4000)


save(m ,m1, m2, m3, m0, file = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Modelos_taxo_all.RData')
save(m ,m1, m2, m3, m0, file = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Modelos_taxo_all.RData')

load('C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Modelos_taxo.RData')
load('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Modelos_taxo.RData')

ranef_families<-data.frame('Family' = rownames(ranef(m)$`FAMILIA_WFO:Order`),
                           'Beta Red List' = unlist(ranef(m3)$`FAMILIA_WFO:Order`),
                           'Beta Endemic' = unlist(ranef(m)$`FAMILIA_WFO:Order`), 
                           'Beta RGR' = unlist(ranef(m1)$`FAMILIA_WFO:Order`), 
                           'Beta HS' = unlist(ranef(m2)$`FAMILIA_WFO:Order`), 
                           'Beta LA' = unlist(ranef(m3)$`FAMILIA_WFO:Order`)
                           )
ranef_families$Family<-unlist(lapply(strsplit(ranef_families$Family, ':'), function(x) x[1]))

dat_var<-data.frame('Var' = c(unlist(get_variance(m0)[c(3, 6)]),
                              unlist(get_variance(m)[c(3, 6)]),
                              unlist(get_variance(m1)[c(3, 6)]),
                              unlist(get_variance(m2)[c(3, 6)]),
                              unlist(get_variance(m3)[c(3, 6)]))) 
var_prop<-dat_var %>%
  mutate(Rarity = rep(c('Endemic', 'RGR', 'HS', 'LA', 'Red List'), each = 4)) %>%
  mutate(Level = rep(c('Residual', 'Genus', 'Family', 'Order'), times = 5)) %>%
  group_by(Rarity) %>%
  mutate(Var_prop = Var/sum(Var)) %>%
  ungroup() %>%
  mutate(Rarity = factor(Rarity, levels = c('LA', 'HS', 'RGR', 'Endemic', 'Red List'))) %>%
  mutate(Level = factor(Level, levels = c('Residual', 'Genus', 'Family', 'Order'))) %>%
  ggplot(aes(x = Rarity, y = 100*Var_prop, fill = Level))+
  geom_bar(stat = 'identity', color = 'black')+
  theme_bw()+
  ylab('Variation explained (%)')+
  xlab('Rarity type')+
  ggtitle('C')+
  theme(text = element_text(size = 20),
        axis.title.y=element_blank(),
        legend.position = 'bottom')+
  scale_fill_manual(values=c( 'gray90', "#fb9a99", "#6a3d9a", "#b15928", ''))

rand_total<-data.frame('Var' = c(unlist(get_variance(m)[c(2, 3)]),
                                 unlist(get_variance(m1)[c(2,3 )]),
                                 unlist(get_variance(m2)[c(2, 3)]),
                                 unlist(get_variance(m3)[c(2, 3)]))) %>%
  mutate(Rarity = rep(c('Endemic', 'RGR', 'HS', 'LA'), each = 2)) %>%
  mutate(Level = rep(c('Random', 'Residual'), times = 4)) %>% 
  group_by(Rarity) %>%
  mutate(Var_prop = 100*Var/sum(Var)) %>%
  ungroup()


grid.arrange(ses_loss_plot, var_prop+coord_flip(),
             nrow = 1)
