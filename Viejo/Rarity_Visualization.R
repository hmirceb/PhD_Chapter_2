library(treedataverse)
library(tidyverse)
library(phylosignal)
library(phylobase)
library(picante)
library(ggnewscale)
library(ggtreeExtra)
library(popbio)
library(ggtree)
library(phytools)

load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/FiloSignal_Rareza_scaled.RData")
load('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/FiloSignal_Rareza_scaled.RData')
load("C:/Users/18172844S/Dropbox/DATA__LAB/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
load("~/Dropbox/DATA__LAB/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
load("C:/Users/18172844S/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")
load("~/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")

rm(tab_v2, tre, w1)

data<-eiv %>% 
  dplyr::select(c(TAXON_REF_PYR_MOD, RA, RA_max, LA, n_habs, B, B_h, PS, Endemic)) %>%
  filter(complete.cases(.))  %>%
  mutate(RA_max = as.numeric(RA_max)) %>%
  mutate(Endemic = ifelse(Endemic == 'End', 1, 0)) %>%
  left_join(EDs, 'TAXON_REF_PYR_MOD') 

lipa_all_df<-as.data.frame(mean.list(lapply(lipa_all, function(x)
  cbind(x$lipa, x$p.value))))
names(lipa_all_df)<-paste('LIPA', names(lipa_all_df), sep = '_')
names(lipa_all_df)[5:8]<-paste(names(lipa_all_df)[5:8], 'pval', sep = '_')
lipa_all_df$TAXON_REF_PYR_MOD<-rownames(lipa_all_df)

data_all<-inner_join(data, lipa_all_df, 'TAXON_REF_PYR_MOD') %>%
  mutate(Endemic = as.factor(Endemic))

tree_figure<-keep.tip(trees_sp_yule_10[[1]], data_all$TAXON_REF_PYR_MOD)

#### Arbol base todo ####

data_all$end_color<-ifelse(data_all$Endemic == 1 & data_all$LIPA_Endemic_pval <= 0.05, 
                           '1', 
                           ifelse(data_all$Endemic == 1 & data_all$LIPA_Endemic_pval > 0.05,
                                  '2', '3'))
q<-ggtree(tree_figure, layout = 'fan') %<+% data_all

##### Nodos todo #####

ra<-apply(do.call('cbind', RA_node_all), 1, mean)

RA_nodes<-as.data.frame(cbind(names(ra), ra, sapply(as.numeric(names(RA_node_all[[1]][,1])),
                                             function(x) phytools::nodeheight(tree = tree_figure, node = x))))
RA_nodes<-RA_nodes %>% mutate_all(as.numeric) %>% rename(var = ra)
RA_nodes_top<-RA_nodes[RA_nodes$var >= quantile(RA_nodes$var, seq(0, 1, 0.01))[100],]

hs<-apply(do.call('cbind', HS_node_all), 1, mean)

HS_nodes<-as.data.frame(cbind(names(hs), hs, sapply(as.numeric(names(HS_node_all[[1]][,1])),
                                             function(x) phytools::nodeheight(tree = tree_figure, node = x))))
HS_nodes<-HS_nodes %>% mutate_all(as.numeric) %>% rename(var = hs)
HS_nodes_top<-HS_nodes[HS_nodes$var >= quantile(HS_nodes$var, seq(0, 1, 0.01))[100],]

la<-apply(do.call('cbind', LA_node_all), 1, mean)

LA_nodes<-as.data.frame(cbind(names(la), la, sapply(as.numeric(names(LA_node_all[[1]][,1])),
                                             function(x) phytools::nodeheight(tree = tree_figure, node = x))))
LA_nodes<-LA_nodes %>% mutate_all(as.numeric) %>% rename(var = la)
LA_nodes_top<-LA_nodes[LA_nodes$var >= quantile(LA_nodes$var, seq(0, 1, 0.01))[100],]

end<-apply(do.call('cbind', End_node_all), 1, mean)

End_nodes<-as.data.frame(cbind(names(end), end, sapply(as.numeric(names(End_node_all[[1]][,1])),
                                               function(x) phytools::nodeheight(tree = tree_figure, node = x))))
End_nodes<-End_nodes %>% mutate_all(as.numeric) %>% rename(var = end)
End_nodes_top<-End_nodes[End_nodes$var >= quantile(End_nodes$var, seq(0, 1, 0.01))[100],]


#q<-q+geom_nodepoint(aes(subset = node %in% RA_nodes_top$V1), color = 'firebrick', size = 3, alpha = 0.5)+
  geom_nodepoint(aes(subset = node %in% LA_nodes_top$V1), color = 'forestgreen', size = 3, alpha = 0.5)+
  geom_nodepoint(aes(subset = node %in% HS_nodes_top$V1), color = 'steelblue', size = 3, alpha = 0.5)+
  geom_nodepoint(aes(subset = node %in% End_nodes_top$V1), color = '#FF8C00BF', size = 3, alpha = 0.5)
q


cols_endemic<-c("1" = "#FF8C00BF", "2" = "#979797", '3' = "#ffffff")
cols_RA<-c("TRUE" = "firebrick", "FALSE" = "#b2bed14D")
cols_HS<-c("TRUE" = "steelblue", "FALSE" = "#b2bed14D")
cols_LA<-c("TRUE" = "forestgreen", "FALSE" = "#b2bed14D")

p<-q + geom_tippoint(
  mapping=aes(color= end_color, fill = end_color, size = end_color,
              alpha = end_color, group = rev(end_color), x = x+12,
              stroke = 1),
  position="identity")+
  scale_color_manual(values = cols_endemic)+
  scale_fill_manual(values = cols_endemic)+
  scale_alpha_manual(values = c('1' = 1, '2' = 0.4, '3' = 0))+
  scale_size_manual(values = c('1' = 3, '2' = 2, '3' = 0))+
  theme(legend.position = 'none')

p1<-p+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(RA_max),
                color = LIPA_RA_max_pval < 0.05, fill = LIPA_RA_max_pval < 0.05),
  stat = 'identity', orientation = 'y', offset = 0.12
)+scale_color_manual(values = cols_RA, name = "RA LIPA sig")+
  scale_fill_manual(values = cols_RA, name = "RA LIPA sig")

p2<-p1+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(B_h),
                color = LIPA_B_h_pval < 0.05, fill = LIPA_B_h_pval < 0.05),
  stat = 'identity', orientation = 'y'
)+scale_color_manual(values = cols_HS, name = "HS LIPA sig")+
  scale_fill_manual(values = cols_HS, name = "HS LIPA sig")

p3<-p2+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(LA),
                color = LIPA_LA_pval < 0.05, fill = LIPA_LA_pval < 0.05),
  stat = 'identity', orientation = 'y'
)+scale_color_manual(values = cols_LA, name = "LA LIPA sig")+
  scale_fill_manual(values = c('Endemic' = '#FF8C00BF', 
                               'RGR' = 'firebrick',
                               'HS' = 'steelblue',
                               'LS' = 'forestgreen'), 
                    name = "Rarity type:")
p3

ggsave(p3, filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/arbol.jpeg', height = 29, width = 29, units = 'cm')
ggsave(p3, filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/arbol.jpeg', height = 29, width = 29, units = 'cm')


#### Arbol circular angiospermas ####


lipa_angios_df<-as.data.frame(mean.list(lapply(lipa_angios, function(x)
  cbind(x$lipa, x$p.value))))
names(lipa_angios_df)<-paste('LIPA', names(lipa_angios_df), sep = '_')
names(lipa_angios_df)[5:8]<-paste(names(lipa_angios_df)[5:8], 'pval', sep = '_')
lipa_angios_df$TAXON_REF_PYR_MOD<-rownames(lipa_angios_df)

data_angios<-inner_join(data, lipa_angios_df, 'TAXON_REF_PYR_MOD')

tree_figure_ang<-keep.tip(trees_sp_yule_10[[1]], data_angios$TAXON_REF_PYR_MOD)
data_angios$end_color<-ifelse(data_angios$Endemic == 1 & data_angios$LIPA_Endemic_pval <= 0.05, 
                           '1', 
                           ifelse(data_angios$Endemic == 1 & data_angios$LIPA_Endemic_pval > 0.05,
                                  '2', '3'))
q_ang<-ggtree(tree_figure_ang, layout = 'fan') %<+% data_angios


##### Nodos angios #####

ra_ang<-apply(do.call('cbind', RA_node_ang), 1, mean)

RA_nodes_ang<-as.data.frame(cbind(names(ra_ang), ra_ang, sapply(as.numeric(names(RA_node_ang[[1]][,1])),
                                                                function(x) phytools::nodeheight(tree = tree_figure_ang, node = x))))
RA_nodes_ang<-RA_nodes_ang %>% mutate_all(as.numeric) %>% rename(var = ra_ang)
RA_nodes_top_ang<-RA_nodes_ang[RA_nodes_ang$var >= quantile(RA_nodes_ang$var, seq(0, 1, 0.01))[100],]

hs_ang<-apply(do.call('cbind', HS_node_ang), 1, mean)

HS_nodes_ang<-as.data.frame(cbind(names(hs_ang), hs_ang, sapply(as.numeric(names(HS_node_ang[[1]][,1])),
                                                                function(x) phytools::nodeheight(tree = tree_figure_ang, node = x))))
HS_nodes_ang<-HS_nodes_ang %>% mutate_all(as.numeric) %>% rename(var = hs_ang)
HS_nodes_top_ang<-HS_nodes_ang[HS_nodes_ang$var >= quantile(HS_nodes_ang$var, seq(0, 1, 0.01))[100],]

la_ang<-apply(do.call('cbind', LA_node_ang), 1, mean)

LA_nodes_ang<-as.data.frame(cbind(names(la_ang), la_ang, sapply(as.numeric(names(LA_node_ang[[1]][,1])),
                                                                function(x) phytools::nodeheight(tree = tree_figure_ang, node = x))))
LA_nodes_ang<-LA_nodes_ang %>% mutate_all(as.numeric) %>% rename(var = la_ang)
LA_nodes_top_ang<-LA_nodes_ang[LA_nodes_ang$var >= quantile(LA_nodes_ang$var, seq(0, 1, 0.01))[100],]

end_ang<-apply(do.call('cbind', End_node_ang), 1, mean)

End_nodes_ang<-as.data.frame(cbind(names(end_ang), end_ang, sapply(as.numeric(names(End_node_ang[[1]][,1])),
                                                                   function(x) phytools::nodeheight(tree = tree_figure_ang, node = x))))
End_nodes_ang<-End_nodes_ang %>% mutate_all(as.numeric) %>% rename(var = end_ang)
End_nodes_top_ang<-End_nodes_ang[End_nodes_ang$var >= quantile(End_nodes_ang$var, seq(0, 1, 0.01))[100],]


q_ang<-q_ang+geom_nodepoint(aes(subset = node %in% RA_nodes_top_ang$V1), color = 'firebrick', size = 3, alpha = 0.5)+
  geom_nodepoint(aes(subset = node %in% LA_nodes_top_ang$V1), color = 'forestgreen', size = 3, alpha = 0.5)+
  geom_nodepoint(aes(subset = node %in% HS_nodes_top_ang$V1), color = 'steelblue', size = 3, alpha = 0.5)+
  geom_nodepoint(aes(subset = node %in% End_nodes_top_ang$V1), color = '#FF8C00BF', size = 3, alpha = 0.5)
q_ang

cols_endemic<-c("1" = "#FF8C00BF", "2" = "#979797", '3' = "#ffffff")
cols_RA<-c("TRUE" = "firebrick", "FALSE" = "#b2bed14D")
cols_HS<-c("TRUE" = "steelblue", "FALSE" = "#b2bed14D")
cols_LA<-c("TRUE" = "forestgreen", "FALSE" = "#b2bed14D")

p_ang<-q_ang + geom_tippoint(
  mapping=aes(color= end_color, fill = end_color, size = end_color,
              alpha = end_color, group = rev(end_color), x = x+12,
              stroke = 1),
  position="identity")+
  scale_color_manual(values = cols_endemic)+
  scale_fill_manual(values = cols_endemic)+
  scale_alpha_manual(values = c('1' = 1, '2' = 0.4, '3' = 0))+
  scale_size_manual(values = c('1' = 3, '2' = 2, '3' = 0))

p1_ang<-p_ang+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(RA_max),
                color = LIPA_RA_max_pval < 0.05, fill = LIPA_RA_max_pval < 0.05),
  stat = 'identity', orientation = 'y', offset = 0.12
)+scale_color_manual(values = cols_RA, name = "RA LIPA sig")+
  scale_fill_manual(values = cols_RA, name = "RA LIPA sig")

p2_ang<-p1_ang+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(B_h),
                color = LIPA_B_h_pval < 0.05, fill = LIPA_B_h_pval < 0.05),
  stat = 'identity', orientation = 'y'
)+scale_color_manual(values = cols_HS, name = "HS LIPA sig")+
  scale_fill_manual(values = cols_HS, name = "HS LIPA sig")

p3_ang<-p2_ang+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(LA),
                color = LIPA_LA_pval < 0.05, fill = LIPA_LA_pval < 0.05),
  stat = 'identity', orientation = 'y'
)+scale_color_manual(values = cols_LA, name = "LA LIPA sig")+
  scale_fill_manual(values = cols_LA, name = "LA LIPA sig")+
  theme(legend.position = 'none')
p3_ang


#### Arbol circular gimnos ####
lipa_gimn_df<-as.data.frame(mean.list(lapply(lipa_gimn, function(x)
  cbind(x$lipa, x$p.value))))
names(lipa_gimn_df)<-paste('LIPA', names(lipa_gimn_df), sep = '_')
names(lipa_gimn_df)[5:8]<-paste(names(lipa_gimn_df)[5:8], 'pval', sep = '_')
lipa_gimn_df$TAXON_REF_PYR_MOD<-rownames(lipa_gimn_df)

data_gimn<-inner_join(data, lipa_gimn_df, 'TAXON_REF_PYR_MOD')

tree_figure_gimn<-keep.tip(trees_sp_yule_10[[1]], data_gimn$TAXON_REF_PYR_MOD)
data_gimn$end_color<-ifelse(data_gimn$Endemic == 1 & data_gimn$LIPA_Endemic_pval <= 0.05, 
                              '1', 
                              ifelse(data_gimn$Endemic == 1 & data_gimn$LIPA_Endemic_pval > 0.05,
                                     '2', '3'))
q_gimn<-ggtree(tree_figure_gimn, layout = 'fan') %<+% data_gimn


##### Nodos gimnos #####

ra_gimn<-apply(do.call('cbind', RA_node_gimn), 1, mean)

RA_nodes_gimn<-as.data.frame(cbind(names(ra_gimn), ra_gimn, sapply(as.numeric(names(RA_node_gimn[[1]][,1])),
                                                                function(x) phytools::nodeheight(tree = tree_figure_gimn, node = x))))
RA_nodes_gimn<-RA_nodes_gimn %>% mutate_all(as.numeric) %>% rename(var = ra_gimn)
RA_nodes_top_gimn<-RA_nodes_gimn[RA_nodes_gimn$var >= quantile(RA_nodes_gimn$var, seq(0, 1, 0.05))[20],]

hs_gimn<-apply(do.call('cbind', HS_node_gimn), 1, mean)

HS_nodes_gimn<-as.data.frame(cbind(names(hs_gimn), hs_gimn, sapply(as.numeric(names(HS_node_gimn[[1]][,1])),
                                                                function(x) phytools::nodeheight(tree = tree_figure_gimn, node = x))))
HS_nodes_gimn<-HS_nodes_gimn %>% mutate_all(as.numeric) %>% rename(var = hs_gimn)
HS_nodes_top_gimn<-HS_nodes_gimn[HS_nodes_gimn$var >= quantile(HS_nodes_gimn$var, seq(0, 1, 0.05))[20],]

la_gimn<-apply(do.call('cbind', LA_node_gimn), 1, mean)

LA_nodes_gimn<-as.data.frame(cbind(names(la_gimn), la_gimn, sapply(as.numeric(names(LA_node_gimn[[1]][,1])),
                                                                function(x) phytools::nodeheight(tree = tree_figure_gimn, node = x))))
LA_nodes_gimn<-LA_nodes_gimn %>% mutate_all(as.numeric) %>% rename(var = la_gimn)
LA_nodes_top_gimn<-LA_nodes_gimn[LA_nodes_gimn$var >= quantile(LA_nodes_gimn$var, seq(0, 1, 0.05))[20],]

end_gimn<-apply(do.call('cbind', End_node_gimn), 1, mean)

End_nodes_gimn<-as.data.frame(cbind(names(end_gimn), end_gimn, sapply(as.numeric(names(End_node_gimn[[1]][,1])),
                                                                   function(x) phytools::nodeheight(tree = tree_figure_gimn, node = x))))
End_nodes_gimn<-End_nodes_gimn %>% mutate_all(as.numeric) %>% rename(var = end_gimn)
End_nodes_top_gimn<-End_nodes_gimn[End_nodes_gimn$var >= quantile(End_nodes_gimn$var, seq(0, 1, 0.05))[20],]


q_gimn<-q_gimn+geom_nodepoint(aes(subset = node %in% RA_nodes_top_gimn$V1), color = 'firebrick', size = 3, alpha = 0.5)+
  geom_nodepoint(aes(subset = node %in% LA_nodes_top_gimn$V1), color = 'forestgreen', size = 3, alpha = 0.5)+
  geom_nodepoint(aes(subset = node %in% HS_nodes_top_gimn$V1), color = 'steelblue', size = 3, alpha = 0.5)
q_gimn

cols_endemic<-c("1" = "#FF8C00BF", "2" = "#979797", '3' = "#ffffff")
cols_RA<-c("TRUE" = "firebrick", "FALSE" = "#b2bed14D")
cols_HS<-c("TRUE" = "steelblue", "FALSE" = "#b2bed14D")
cols_LA<-c("TRUE" = "forestgreen", "FALSE" = "#b2bed14D")

p_gimn<-q_gimn + geom_tippoint(
  mapping=aes(color= end_color, fill = end_color, size = end_color,
              alpha = end_color, group = rev(end_color), x = x+12,
              stroke = 1),
  position="identity")+
  scale_color_manual(values = cols_endemic)+
  scale_fill_manual(values = cols_endemic)+
  scale_alpha_manual(values = c('1' = 1, '2' = 0.4, '3' = 0))+
  scale_size_manual(values = c('1' = 3, '2' = 2, '3' = 0))

p1_gimn<-p_gimn+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(RA_max),
                color = LIPA_RA_max_pval < 0.05, fill = LIPA_RA_max_pval < 0.05),
  stat = 'identity', orientation = 'y', offset = 0.12
)+scale_color_manual(values = cols_RA, name = "RA LIPA sig")+
  scale_fill_manual(values = cols_RA, name = "RA LIPA sig")

p2_gimn<-p1_gimn+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(B_h),
                color = LIPA_B_h_pval < 0.05, fill = LIPA_B_h_pval < 0.05),
  stat = 'identity', orientation = 'y'
)+scale_color_manual(values = cols_HS, name = "HS LIPA sig")+
  scale_fill_manual(values = cols_HS, name = "HS LIPA sig")

p3_gimn<-p2_gimn+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(LA),
                color = LIPA_LA_pval < 0.05, fill = LIPA_LA_pval < 0.05),
  stat = 'identity', orientation = 'y'
)+scale_color_manual(values = cols_LA, name = "LA LIPA sig")+
  scale_fill_manual(values = cols_LA, name = "LA LIPA sig")+
  theme(legend.position = 'none')
p3_gimn

# Boxplot de la senyal filo ####
cols_boxplot<-c('Endemic' = "#FF8C00BF", 'RGR' = 'firebrick', 'HS' = 'steelblue', 'LA' = 'forestgreen')

signal<-rbind(do.call('rbind', lapply(signal_all, function(x) 
  data.frame('rarity' = rownames(x$stat),
             'lambda' = x$stat$Lambda, 
             'p_val' = x$pvalue$Lambda, 
             'group' = 'All'))), 
  do.call('rbind', lapply(signal_angios, function(x) 
    data.frame('rarity' = rownames(x$stat),
               'lambda' = x$stat$Lambda, 
               'p_val' = x$pvalue$Lambda, 
               'group' = 'Angiosperms'))),
  do.call('rbind', lapply(signal_gimn, function(x) 
    data.frame('rarity' = rownames(x$stat),
               'lambda' = x$stat$Lambda, 
               'p_val' = x$pvalue$Lambda, 
               'group' = 'Gymnosperms & Ferns'))))
dat<-signal %>% 
  filter(rarity != 'Endemic') %>%
  bind_rows(data.frame('rarity'= 'Endemic', 
                       'lambda' = 1-c(unlist(lapply(Endemic_PhyloD, function(x) x$DEstimate)),
                                      unlist(lapply(Endemic_PhyloD_ang, function(x) x$DEstimate))),
                       'p_val' = 0, 
                       'group' = rep(c('All', 'Angiosperms'), each = 10))) %>%
  mutate(rarity=recode(rarity, 
                       RA_max = 'RGR',
                       B_h = 'HS')) %>%
  mutate(rarity = factor(rarity, levels = c('Endemic', 'RGR', 'HS', 'LA')))
dat_text<-dat %>%
  aggregate(p_val~rarity+group, ., mean) %>%
  mutate(yloc = max(dat$lambda))

boxplot_signal<-ggplot(dat, aes(x = rarity, y = lambda, fill = rarity))+geom_boxplot()+
  facet_wrap(~group)+theme_bw()+ylab('Signal')+xlab('Rarity type')+
  theme(text = element_text(size = 20),
        strip.text = element_text(face = 'bold'))+
  scale_fill_manual(values = cols_boxplot, name = 'Rarity type')+
  geom_text(data = dat_text, aes(x = rarity, y = yloc+0.05, label = 
                  ifelse(p_val <= 0.05, '***', '')), size = 8)
boxplot_signal
ggsave(boxplot_signal, 
       filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Boxplot_signal.jpeg',
       height = 21, width = 29.7, units = 'cm')

#### ####
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/PD_loss.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/PD_loss.RData")

aaa<-data.frame(
  'Line' = c(obs_loss_end, obs_loss_ra, obs_loss_hs, obs_loss_la),
  'Mean' = c(mean(unlist(rand_loss_end)), 
             mean(unlist(rand_loss_ra)),
             mean(unlist(rand_loss_hs)),
             mean(unlist(rand_loss_la))),
  'SD' = c(sd(unlist(rand_loss_end)), 
           sd(unlist(rand_loss_ra)),
           sd(unlist(rand_loss_hs)),
           sd(unlist(rand_loss_la))),
  'Rarity' = factor(c('Endemic', 'RGR', 'HS', 'LA')))
aaa$ses<-(aaa$Line-aaa$Mean)/aaa$SD

ses_loss_plot<-aaa %>% 
  mutate(Rarity = factor(Rarity, levels = c('LA', 'HS', 'RGR', 'Endemic'))) %>%
  ggplot(aes(y = Rarity))+
  geom_segment(aes(xend = ses, x = 0, yend = Rarity),
               linetype = 'dashed', size = 1)+
  theme_bw()+
  geom_point(aes(x = ses, color = Rarity), size = 5)+
  scale_color_manual(values = c('Endemic' = '#FF8C00', 
                                'RGR' = 'firebrick',
                                'HS' = 'steelblue',
                                'LA' = 'forestgreen'), 
                     name = "Rarity type:", guide = 'none')+
  theme(text = element_text(size = 20))+
  xlab('Standard effect size of PD loss')+
  ylab('Rarity type')+
  annotate(geom = 'rect', xmin = -1.96, 
           xmax = 1.96, ymin = 0, ymax = 5,
           alpha = 0.1)+
  ggtitle('B')

load('C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Modelos_taxo.RData')
load('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Modelos_taxo.RData')

ranef_families<-data.frame('Family' = rownames(ranef(m)$`FAMILIA_WFO:Order`),
                           'Beta Endemic' = unlist(ranef(m)$`FAMILIA_WFO:Order`), 
                           'Beta RGR' = unlist(ranef(m1)$`FAMILIA_WFO:Order`), 
                           'Beta HS' = unlist(ranef(m2)$`FAMILIA_WFO:Order`), 
                           'Beta LA' = unlist(ranef(m3)$`FAMILIA_WFO:Order`))
ranef_families$Family<-unlist(lapply(strsplit(ranef_families$Family, ':'), function(x) x[1]))

dat_var<-data.frame('Var' = c(unlist(get_variance(m)[c(3, 6)]),
                              unlist(get_variance(m1)[c(3, 6)]),
                              unlist(get_variance(m2)[c(3, 6)]),
                              unlist(get_variance(m3)[c(3, 6)]))) 
var_prop<-dat_var %>%
  mutate(Rarity = rep(c('Endemic', 'RGR', 'HS', 'LA'), each = 4)) %>%
  mutate(Level = rep(c('Residual', 'Genus', 'Family', 'Order'), times = 4)) %>%
  group_by(Rarity) %>%
  mutate(Var_prop = Var/sum(Var)) %>%
  ungroup() %>%
  mutate(Rarity = factor(Rarity, levels = c('LA', 'HS', 'RGR', 'Endemic'))) %>%
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
  scale_fill_manual(values=c( 'gray90', "#009E73", "#D55E00", "mediumorchid4"))


#### Analisis importancia de nodos por distancia a la raiz ####
#### Todo ####

smooth_nodes_all<-End_nodes %>%
  bind_rows(RA_nodes, HS_nodes, LA_nodes) %>%
  mutate(rarity = rep(c('Pyrenean Endemic', 'Regional Geographic Range', 'Habitat Specialization', 'Local Abundance'), each = dim(End_nodes)[1])) %>%
  rename(node = V1, height = V3)  %>%
  group_by(rarity) %>%
  mutate(var_scaled = var/sum(var)) %>%
  ungroup() %>%
  ggplot(aes(x = height, y = var_scaled, color = rarity, fill = rarity))+
  theme_bw()+
  scale_color_manual(name = 'Rarity type:', values = c('Pyrenean Endemic' = '#FF8C00BF', 'Regional Geographic Range' = 'firebrick', 'Habitat Specialization' = 'steelblue', 'Local Abundance' = 'forestgreen'))+
  scale_fill_manual(name = 'Rarity type:', values = c('Pyrenean Endemic' = '#FF8C00BF', 'Regional Geographic Range' = 'firebrick', 'Habitat Specialization' = 'steelblue', 'Local Abundance' = 'forestgreen'))+
  geom_smooth(method = 'gam', alpha = 0.08, formula = y~s(x, bs = 'ts'))+
  ylab('Node relative importance')+
  xlab('Node height')+
  theme(text = element_text(size = 20),
        legend.position = 'bottom')+
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.7)+
  geom_point()
smooth_nodes_all

points_nodes_all<-End_nodes %>%
  bind_rows(RA_nodes, HS_nodes, LA_nodes) %>%
  mutate(rarity = rep(c('Pyrenean Endemic', 'Regional Geographic Range', 'Habitat Specialization', 'Local Abundance'), each = dim(End_nodes)[1])) %>%
  rename(node = V1, height = V3)  %>%
  group_by(rarity) %>%
  mutate(var_scaled = var/sum(var)) %>%
  ungroup() %>%
  ggplot(aes(x = height, y = var_scaled, color = rarity, fill = rarity))+
  theme_bw()+
  scale_color_manual(name = 'Rarity type:', values = c('Pyrenean Endemic' = '#FF8C00BF', 'Regional Geographic Range' = 'firebrick', 'Habitat Specialization' = 'steelblue', 'Local Abundance' = 'forestgreen'))+
  scale_fill_manual(name = 'Rarity type:', values = c('Pyrenean Endemic' = '#FF8C00BF', 'Regional Geographic Range' = 'firebrick', 'Habitat Specialization' = 'steelblue', 'Local Abundance' = 'forestgreen'))+
  geom_point()+
  ylab('Node relative importance')+
  xlab('Node height')+
  theme(text = element_text(size = 20))+
  ggtitle('B')


n_bins<-20
nodes_top_count_all<-bind_rows(End_nodes %>% mutate(t_bins = cut(V3, 
                                                             breaks = 10, 
                                                             include.lowest=TRUE)) %>%
                             group_by(t_bins) %>%
                             summarize(n_sig = sum(V1 %in% End_nodes_top$V1)),
                           
                           RA_nodes %>% mutate(t_bins = cut(V3, 
                                                            breaks = 10, 
                                                            include.lowest=TRUE)) %>%
                             group_by(t_bins) %>%
                             summarize(n_sig =sum(V1 %in% RA_nodes_top$V1)),
                           
                           HS_nodes %>% mutate(t_bins = cut(V3, 
                                                            breaks = 10, 
                                                            include.lowest=TRUE)) %>%
                             group_by(t_bins) %>%
                             summarize(n_sig =sum(V1 %in% HS_nodes_top$V1)),
                           
                           LA_nodes %>% mutate(t_bins = cut(V3, 
                                                            breaks = 10, 
                                                            include.lowest=TRUE)) %>%
                             group_by(t_bins) %>%
                             summarize(n_sig =sum(V1 %in% LA_nodes_top$V1))) %>%
  mutate(rarity = rep(c('Endemic', 'RGR', 'HS', 'LA'), each = 10)) %>%
  left_join(
    bind_rows(End_nodes %>% mutate(t_bins = cut(V3, 
                                                breaks = 10, 
                                                include.lowest=TRUE)) %>%
                group_by(t_bins) %>%
                summarize(n_total = length(V1)),
              
              RA_nodes %>% mutate(t_bins = cut(V3, 
                                               breaks = 10, 
                                               include.lowest=TRUE)) %>%
                group_by(t_bins) %>%
                summarize(n_total = length(V1)),
              
              HS_nodes %>% mutate(t_bins = cut(V3, 
                                               breaks = 10, 
                                               include.lowest=TRUE)) %>%
                group_by(t_bins) %>%
                summarize(n_total = length(V1)),
              
              LA_nodes %>% mutate(t_bins = cut(V3, 
                                               breaks = 10, 
                                               include.lowest=TRUE)) %>%
                group_by(t_bins) %>%
                summarize(n_total = length(V1))) %>%
      mutate(rarity = rep(c('Endemic', 'RGR', 'HS', 'LA'), each = 10)),
    c('t_bins', 'rarity')
  ) %>%
  mutate(n_sig_prop = round(100*(n_sig/n_total), 1))

# Porcentaje de nodos top en cada time slice
perc_nodes_bar_all<-nodes_top_count_all %>%
  ggplot(aes(x = t_bins, y = n_sig_prop, fill = rarity))+
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black')+
  scale_fill_manual(name = 'Rarity type:', values = c('Endemic' = '#FF8C00BF', 'RGR' = 'firebrick', 'HS' = 'steelblue', 'LA' = 'forestgreen'))+
  theme_bw()+
  ylab('% of significant nodes')+
  xlab('Node height (binned)')+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15))+
  ggtitle('B')

# Boxplot de importancia de valores por time slice  de todo 
boxplot_nodes_all<-End_nodes %>%
  bind_rows(RA_nodes, HS_nodes, LA_nodes) %>%
  mutate(rarity = rep(c('Pyrenean Endemic', 'Regional Geographic Range', 'Habitat Specialization', 'Local Abundance'), each = dim(End_nodes)[1])) %>%
  rename(node = V1, height = V3)  %>%
  group_by(rarity) %>%
  mutate(var_scaled = var/sum(var)) %>%
  ungroup() %>%
  mutate(t_bins = cut(height, 
                      breaks = 10, 
                      include.lowest=TRUE)) %>%
  ggplot(aes(x = t_bins, y = var_scaled, groups = t_bins, fill = rarity))+
  geom_boxplot()+
  scale_fill_manual(name = 'Rarity type:', values = c('Pyrenean Endemic' = '#FF8C00BF', 'Regional Geographic Range' = 'firebrick', 'Habitat Specialization' = 'steelblue', 'Local Abundance' = 'forestgreen'))+
  theme_bw()+
  ylab('Node relative importance')+
  xlab('Node height (binned)')+
  ggtitle('B')+
  scale_x_discrete(labels = c("[0, 42.1]", "[42.1, 84.2]", "[84.2, 126]", "[126, 168]", "[168, 210]", "[210, 253]",    
                              "[253, 295]", "[337, 379]", "[295, 337]", "[379, 421]"))+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 10))

ggarrange(smooth_nodes_all, points_nodes_all, ncol = 2, common.legend = T, legend = 'bottom')


##### Angios #####

ra_ang<-apply(do.call('cbind', RA_node_ang), 1, mean)
RA_nodes_ang<-as.data.frame(cbind(names(ra_ang), ra_ang, sapply(as.numeric(names(RA_node_ang[[1]][,1])),
                                                                function(x) phytools::nodeheight(tree = tree_figure_ang, node = x))))
RA_nodes_ang<-RA_nodes_ang %>% mutate_all(as.numeric) %>% rename(var = ra_ang)
RA_nodes_top_ang<-RA_nodes[RA_nodes_ang$var >= quantile(RA_nodes_ang$var, seq(0, 1, 0.05))[20],]

hs_ang<-apply(do.call('cbind', HS_node_ang), 1, mean)
HS_nodes_ang<-as.data.frame(cbind(names(hs_ang), hs_ang, sapply(as.numeric(names(HS_node_ang[[1]][,1])),
                                                                function(x) phytools::nodeheight(tree = tree_figure_ang, node = x))))
HS_nodes_ang<-HS_nodes_ang %>% mutate_all(as.numeric) %>% rename(var = hs_ang)
HS_nodes_top_ang<-HS_nodes[HS_nodes_ang$var >= quantile(HS_nodes_ang$var, seq(0, 1, 0.05))[20],]

la_ang<-apply(do.call('cbind', LA_node_ang), 1, mean)
LA_nodes_ang<-as.data.frame(cbind(names(la_ang), la_ang, sapply(as.numeric(names(LA_node_ang[[1]][,1])),
                                                                function(x) phytools::nodeheight(tree = tree_figure_ang, node = x))))
LA_nodes_ang<-LA_nodes_ang %>% mutate_all(as.numeric) %>% rename(var = la_ang)
LA_nodes_top_ang<-LA_nodes_ang[LA_nodes_ang$var >= quantile(LA_nodes_ang$var, seq(0, 1, 0.05))[20],]

end_ang<-apply(do.call('cbind', End_node_ang), 1, mean)
End_nodes_ang<-as.data.frame(cbind(names(end_ang), end_ang, sapply(as.numeric(names(End_node_ang[[1]][,1])),
                                                                   function(x) phytools::nodeheight(tree = tree_figure_ang, node = x))))
End_nodes_ang<-End_nodes_ang %>% mutate_all(as.numeric) %>% rename(var = end_ang)
End_nodes_top_ang<-End_nodes_ang[End_nodes_ang$var >= quantile(End_nodes_ang$var, seq(0, 1, 0.05))[20],]


smooth_nodes_angios<-End_nodes_ang %>%
  bind_rows(RA_nodes_ang, HS_nodes_ang, LA_nodes_ang) %>%
  mutate(rarity = rep(c('Pyrenean Endemic', 'Regional Geographic Range', 'Habitat Specialization', 'Local Abundance'), each = dim(End_nodes_ang)[1])) %>%
  rename(node = V1, height = V3)  %>%
  group_by(rarity) %>%
  mutate(var_scaled = var/sum(var)) %>%
  ungroup() %>%
  ggplot(aes(x = height, y = var_scaled, color = rarity, fill = rarity))+
  theme_bw()+
  scale_color_manual(name = 'Rarity type:', values = c('Pyrenean Endemic' = '#FF8C00BF', 'Regional Geographic Range' = 'firebrick', 'Habitat Specialization' = 'steelblue', 'Local Abundance' = 'forestgreen'))+
  scale_fill_manual(name = 'Rarity type:', values = c('Pyrenean Endemic' = '#FF8C00BF', 'Regional Geographic Range' = 'firebrick', 'Habitat Specialization' = 'steelblue', 'Local Abundance' = 'forestgreen'))+
  geom_smooth(method = 'gam', alpha = 0.08, formula = y~s(x, bs = 'ts'))+
  ylab('Node relative importance')+
  xlab('Node height')+
  theme(text = element_text(size = 20),
        legend.position = 'bottom')+
  geom_point()

boxplot_nodes_angios<-End_nodes_ang %>%
  bind_rows(RA_nodes_ang, HS_nodes_ang, LA_nodes_ang) %>%
  mutate(rarity = rep(c('Endemic', 'RGR', 'HS', 'LA'), each = dim(End_nodes_ang)[1])) %>%
  rename(node = V1, height = V3)  %>%
  group_by(rarity) %>%
  mutate(var_scaled = var/sum(var)) %>%
  ungroup() %>%
  mutate(t_bins = cut(height, 
                      breaks = 10, 
                      include.lowest=TRUE)) %>%
  ggplot(aes(x = t_bins, y = var_scaled, groups = t_bins, fill = rarity))+
  geom_boxplot()+
  scale_fill_manual(name = 'Rarity', values = c('Endemic' = '#FF8C00BF', 'RGR' = 'firebrick', 'HS' = 'steelblue', 'LA' = 'forestgreen'))+
  theme_bw()


##### Gimnos #####
ra_gimn<-apply(do.call('cbind', RA_node_gimn), 1, mean)
RA_nodes_gimn<-as.data.frame(cbind(names(ra_gimn), ra_gimn, sapply(as.numeric(names(RA_node_gimn[[1]][,1])),
                                                                   function(x) phytools::nodeheight(tree = tree_figure_gimn, node = x))))
RA_nodes_gimn<-RA_nodes_gimn %>% mutate_all(as.numeric) %>% rename(var = ra_gimn)
RA_nodes_top_gimn<-RA_nodes[RA_nodes_gimn$var >= quantile(RA_nodes_gimn$var, seq(0, 1, 0.05))[20],]

hs_gimn<-apply(do.call('cbind', HS_node_gimn), 1, mean)
HS_nodes_gimn<-as.data.frame(cbind(names(hs_gimn), hs_gimn, sapply(as.numeric(names(HS_node_gimn[[1]][,1])),
                                                                   function(x) phytools::nodeheight(tree = tree_figure_gimn, node = x))))
HS_nodes_gimn<-HS_nodes_gimn %>% mutate_all(as.numeric) %>% rename(var = hs_gimn)
HS_nodes_top_gimn<-HS_nodes[HS_nodes_gimn$var >= quantile(HS_nodes_gimn$var, seq(0, 1, 0.05))[20],]

la_gimn<-apply(do.call('cbind', LA_node_gimn), 1, mean)
LA_nodes_gimn<-as.data.frame(cbind(names(la_gimn), la_gimn, sapply(as.numeric(names(LA_node_gimn[[1]][,1])),
                                                                   function(x) phytools::nodeheight(tree = tree_figure_gimn, node = x))))
LA_nodes_gimn<-LA_nodes_gimn %>% mutate_all(as.numeric) %>% rename(var = la_gimn)
LA_nodes_top_gimn<-LA_nodes_gimn[LA_nodes_gimn$var >= quantile(LA_nodes_gimn$var, seq(0, 1, 0.05))[20],]

end_gimn<-apply(do.call('cbind', End_node_gimn), 1, mean)
End_nodes_gimn<-as.data.frame(cbind(names(end_gimn), end_gimn, sapply(as.numeric(names(End_node_gimn[[1]][,1])),
                                                                      function(x) phytools::nodeheight(tree = tree_figure_gimn, node = x))))
End_nodes_gimn<-End_nodes_gimn %>% mutate_all(as.numeric) %>% rename(var = end_gimn)
End_nodes_top_gimn<-End_nodes_gimn[End_nodes_gimn$var >= quantile(End_nodes_gimn$var, seq(0, 1, 0.05))[20],]


smooth_nodes_gimn<-End_nodes_gimn %>%
  bind_rows(RA_nodes_gimn, HS_nodes_gimn, LA_nodes_gimn) %>%
  mutate(rarity = rep(c('Pyrenean Endemic', 'Regional Geographic Range', 'Habitat Specialization', 'Local Abundance'), each = dim(End_nodes_gimn)[1])) %>%
  rename(node = V1, height = V3)  %>%
  group_by(rarity) %>%
  mutate(var_scaled = var/sum(var)) %>%
  ungroup() %>%
  ggplot(aes(x = height, y = var_scaled, color = rarity, fill = rarity))+
  theme_bw()+
  scale_color_manual(name = 'Rarity type:', values = c('Pyrenean Endemic' = '#FF8C00BF', 'Regional Geographic Range' = 'firebrick', 'Habitat Specialization' = 'steelblue', 'Local Abundance' = 'forestgreen'))+
  scale_fill_manual(name = 'Rarity type:', values = c('Pyrenean Endemic' = '#FF8C00BF', 'Regional Geographic Range' = 'firebrick', 'Habitat Specialization' = 'steelblue', 'Local Abundance' = 'forestgreen'))+
  geom_smooth(method = 'gam', alpha = 0.08, formula = y~s(x, bs = 'ts'))+
  ylab('Node relative importance')+
  xlab('Node height')+
  theme(text = element_text(size = 20),
        legend.position = 'bottom')+
  geom_point()

ggpubr::ggarrange(smooth_nodes_angios, smooth_nodes_gimn)


##### Arbol circular con las familias que tienen especies con LIPA significativos #####


data<-eiv %>% 
  dplyr::select(c(TAXON_REF_PYR_MOD, RA, RA_max, LA, n_habs, B, B_h, PS, Endemic)) %>%
  filter(complete.cases(.))  %>%
  mutate(RA_max = as.numeric(RA_max)) %>%
  mutate(Endemic = ifelse(Endemic == 'End', 1, 0))

lipa_all_df<-as.data.frame(mean.list(lapply(lipa_all, function(x)
  cbind(x$lipa, x$p.value))))
names(lipa_all_df)<-paste('LIPA', names(lipa_all_df), sep = '_')
names(lipa_all_df)[5:8]<-paste(names(lipa_all_df)[5:8], 'pval', sep = '_')
lipa_all_df$TAXON_REF_PYR_MOD<-rownames(lipa_all_df)

data_all<-inner_join(data, lipa_all_df, 'TAXON_REF_PYR_MOD') %>%
  mutate(Endemic = as.factor(Endemic))

tree_figure<-keep.tip(trees_sp_yule_10[[1]], data_all$TAXON_REF_PYR_MOD)

#### Arbol base todo ####

data_all$end_color<-ifelse(data_all$Endemic == 1 & data_all$LIPA_Endemic_pval <= 0.05, 
                           '1', 
                           ifelse(data_all$Endemic == 1 & data_all$LIPA_Endemic_pval > 0.05,
                                  '2', '3'))
q<-ggtree(tree_figure, layout = 'fan') %<+% data_all


cols_endemic<-c("1" = "#FF8C00BF", "2" = "#979797", '3' = "#ffffff")
cols_RA<-c("TRUE" = "firebrick", "FALSE" = "#b2bed14D")
cols_HS<-c("TRUE" = "steelblue", "FALSE" = "#b2bed14D")
cols_LA<-c("TRUE" = "forestgreen", "FALSE" = "#b2bed14D")

p<-q + geom_tippoint(
  mapping=aes(color= end_color, fill = end_color, size = end_color,
              alpha = end_color, group = rev(end_color), x = x+12,
              stroke = 1),
  position="identity")+
  scale_color_manual(values = cols_endemic)+
  scale_fill_manual(values = cols_endemic)+
  scale_alpha_manual(values = c('1' = 1, '2' = 0.4, '3' = 0))+
  scale_size_manual(values = c('1' = 3, '2' = 2, '3' = 0))+
  theme(legend.position = 'none')

p1<-p+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(RA_max),
                color = LIPA_RA_max_pval < 0.05, fill = LIPA_RA_max_pval < 0.05),
  stat = 'identity', orientation = 'y', offset = 0.12
)+scale_color_manual(values = cols_RA, name = "RA LIPA sig")+
  scale_fill_manual(values = cols_RA, name = "RA LIPA sig")

p2<-p1+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(B_h),
                color = LIPA_B_h_pval < 0.05, fill = LIPA_B_h_pval < 0.05),
  stat = 'identity', orientation = 'y'
)+scale_color_manual(values = cols_HS, name = "HS LIPA sig")+
  scale_fill_manual(values = cols_HS, name = "HS LIPA sig")

p3<-p2+new_scale_fill()+new_scale_color()+geom_fruit(
  geom = geom_bar,
  mapping = aes(y = TAXON_REF_PYR_MOD, x = scale(LA),
                color = LIPA_LA_pval < 0.05, fill = LIPA_LA_pval < 0.05),
  stat = 'identity', orientation = 'y'
)+scale_color_manual(values = cols_LA, name = "LA LIPA sig")+
  scale_fill_manual(values = c('Endemic' = '#FF8C00BF', 
                               'RGR' = 'firebrick',
                               'HS' = 'steelblue',
                               'LS' = 'forestgreen'), 
                    name = "Rarity type:")


data_fams<-fp_fams %>%
  rename(TAXON_REF_PYR_MOD = TAXON_RP) %>%
  right_join(data_all, 'TAXON_REF_PYR_MOD') %>%
  select(-c(B, PS, LIPA_Endemic, LIPA_RA_max, LIPA_LA, LIPA_B_h))


familias_end<-data_fams %>%
  filter(LIPA_Endemic_pval <= 0.05 & Endemic != 0) %>%
  pull(FAMILIA_WFO)

nodes_end<-list()
for (i in 1:length(familias_end)) {
  sps_temp<-data_fams[data_fams$FAMILIA_WFO == familias_end[i],]$TAXON_REF_PYR_MOD
  nodes_end[[i]]<-as.numeric(findMRCA(tips = sps_temp, tree = tree_figure))
}
nodes_end[sapply(nodes_end, function(x) length(x) == 0)] <- NA
node_families_end<-data.frame('fam' = familias_end, 'mrca' = unlist(nodes_end))
node_families_end<-unique(node_families_end[complete.cases(node_families_end$mrca),])
node_families_end$color<-'#FF8C00BF'
node_families_end$rarity<-'end'


familias_hs<-data_fams %>%
  filter(LIPA_B_h_pval <= 0.05 & scale(B_h) < 0) %>%
  pull(FAMILIA_WFO)

nodes_hs<-list()
for (i in 1:length(familias_hs)) {
  sps_temp<-data_fams[data_fams$FAMILIA_WFO == familias_hs[i],]$TAXON_REF_PYR_MOD
  nodes_hs[[i]]<-as.numeric(findMRCA(tips = sps_temp, tree = tree_figure))
}
nodes_hs[sapply(nodes_hs, function(x) length(x) == 0)] <- NA
node_families_hs<-data.frame('fam' = familias_hs, 'mrca' = unlist(nodes_hs))
node_families_hs<-unique(node_families_hs[complete.cases(node_families_hs$mrca),])
node_families_hs$color<-'steelblue'
node_families_hs$rarity<-'hs'

familias_ra<-data_fams %>%
  filter(LIPA_RA_max_pval <= 0.05 & scale(RA_max) < 0) %>%
  pull(FAMILIA_WFO)

nodes_ra<-list()
for (i in 1:length(familias_ra)) {
  sps_temp<-data_fams[data_fams$FAMILIA_WFO == familias_ra[i],]$TAXON_REF_PYR_MOD
  nodes_ra[[i]]<-as.numeric(findMRCA(tips = sps_temp, tree = tree_figure))
}
nodes_ra[sapply(nodes_ra, function(x) length(x) == 0)] <- NA
node_families_ra<-data.frame('fam' = familias_ra, 'mrca' = unlist(nodes_ra))
node_families_ra<-unique(node_families_ra[complete.cases(node_families_ra$mrca),])
node_families_ra$color<-'firebrick'
node_families_ra$rarity<-'ra'

familias_la<-data_fams %>%
  filter(LIPA_LA_pval <= 0.05 & scale(LA) < 0) %>%
  pull(FAMILIA_WFO)

nodes_la<-list()
for (i in 1:length(familias_la)) {
  sps_temp<-data_fams[data_fams$FAMILIA_WFO == familias_la[i],]$TAXON_REF_PYR_MOD
  nodes_la[[i]]<-as.numeric(findMRCA(tips = sps_temp, tree = tree_figure))
}
nodes_la[sapply(nodes_la, function(x) length(x) == 0)] <- NA
node_families_la<-data.frame('fam' = familias_la, 'mrca' = unlist(nodes_la))
node_families_la<-unique(node_families_la[complete.cases(node_families_la$mrca),])
node_families_la$color<-'forestgreen'
node_families_la$rarity<-'la'

nodos<-bind_rows(node_families_end, node_families_ra, node_families_hs, node_families_la)
nodos$fam[duplicated(nodos$fam)]<-NA

aa<-p3
for (i in 1:dim(nodos[nodos$rarity == 'end',])[1]) {
  n<-nodos[nodos$rarity == 'end',]
  input<-ifelse(is.na(n[i,1]), 
                paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  n[i,1], "'", 
                      ", geom = 'text'" ,', alpha = 0, align = F, offset.text = 100, barsize = 1.5, offset = 320)',
                      sep = ''), 
                paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  n[i,1], "'", 
                      ", geom = 'text'" ,", angle = 'auto'" ,', align = F, offset.text = 100, barsize = 1.5, offset = 320)',
                      sep = ''))
  input2<-paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  '', "'",
                ', color = ', "'", n[i,3],"'", ', align = T, offset.text = 100, barsize = 1.5, offset = 320)',
                sep = '')
  aa<-aa+eval(parse(text = input))+
    eval(parse(text = input2))
}

for (i in 1:dim(nodos[nodos$rarity == 'ra',])[1]) {
  n<-nodos[nodos$rarity == 'ra',]
  input<-ifelse(is.na(n[i,1]), 
                paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  n[i,1], "'", 
                      ", geom = 'text'" ,', alpha = 1, align = F, offset.text = 100, barsize = 1.5, offset = 340)',
                      sep = ''), 
                paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  n[i,1], "'", 
                      ", geom = 'text'" ,", angle = 'auto'" ,', align = F, offset.text = 100, barsize = 1.5, offset = 340)',
                      sep = ''))
  input2<-paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  '', "'",
                ', color = ', "'", n[i,3],"'", ', align = T, offset.text = 100, barsize = 1.5, offset = 340)',
                sep = '')
  aa<-aa+eval(parse(text = input))+
    eval(parse(text = input2))
}


for (i in 1:dim(nodos[nodos$rarity == 'hs',])[1]) {
  n<-nodos[nodos$rarity == 'hs',]
  input<-ifelse(is.na(n[i,1]), 
                paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  n[i,1], "'", 
                      ", geom = 'text'" ,', alpha = 0, align = F, offset.text = 100, barsize = 1.5, offset = 360)',
                      sep = ''), 
                paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  n[i,1], "'", 
                      ", geom = 'text'" ,", angle = 'auto'" ,', align = F, offset.text = 100, barsize = 1.5, offset = 360)',
                      sep = ''))
  input2<-paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  '', "'",
                ', color = ', "'", n[i,3],"'", ', align = T, offset.text = 100, barsize = 1.5, offset = 360)',
                sep = '')
  aa<-aa+eval(parse(text = input))+
    eval(parse(text = input2))
}

for (i in 1:dim(nodos[nodos$rarity == 'la',])[1]) {
  n<-nodos[nodos$rarity == 'la',]
  input<-ifelse(is.na(n[i,1]), 
                paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  n[i,1], "'", 
                      ", geom = 'text'" ,', alpha = 0, align = F, offset.text = 100, barsize = 1.5, offset = 380)',
                      sep = ''), 
                paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  n[i,1], "'", 
                      ", geom = 'text'" ,", angle = -20" ,', align = F, offset.text = 100, barsize = 1.5, offset = 380)',
                      sep = ''))
  input2<-paste('geom_cladelabel(node = ', n[i,2], ', label = ', "'",  '', "'",
                ', color = ', "'", n[i,3],"'", ', align = T, offset.text = 100, barsize = 1.5, offset = 380)',
                sep = '')
  aa<-aa+eval(parse(text = input))+
    eval(parse(text = input2))
}

aa
ggsave(aa, filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/arbol_fams.jpeg', height = 29, width = 29, units = 'cm')
ggsave(aa, filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/arbol_fams.jpeg', height = 29, width = 29, units = 'cm')


a<-eiv %>%
  dplyr::select(c(TAXON_REF_PYR_MOD, FAMILIA_WFO, Order)) %>%
  left_join(data_all, 'TAXON_REF_PYR_MOD') %>%
  filter(scale(B_h) < 0 & LIPA_B_h_pval < 0.05) %>%
  mutate(Rarity = 'HS') %>%
  bind_rows(eiv %>%
              dplyr::select(c(TAXON_REF_PYR_MOD, FAMILIA_WFO, Order)) %>%
              left_join(data_all, 'TAXON_REF_PYR_MOD') %>%
              filter(scale(RA_max) < 0 & LIPA_RA_max_pval < 0.05) %>%
              mutate(Rarity = 'RGR')) %>%
  bind_rows(eiv %>%
              dplyr::select(c(TAXON_REF_PYR_MOD, FAMILIA_WFO, Order)) %>%
              left_join(data_all, 'TAXON_REF_PYR_MOD') %>%
              filter(scale(LA) < 0 & LIPA_LA_pval < 0.05) %>%
              mutate(Rarity = 'LA')) %>%
  bind_rows(eiv %>%
              dplyr::select(c(TAXON_REF_PYR_MOD, FAMILIA_WFO, Order)) %>%
              left_join(data_all, 'TAXON_REF_PYR_MOD') %>%
              filter(Endemic == 1 & LIPA_Endemic_pval < 0.05) %>%
              mutate(Rarity = 'Endemic')) %>%
  distinct() %>%
  dplyr::select(-c(RA, n_habs, B, PS, LIPA_Endemic, LIPA_RA_max,
            LIPA_LA, LIPA_B_h, LIPA_Endemic_pval, LIPA_RA_max_pval,
            LIPA_LA_pval,LIPA_B_h_pval)) %>%
  left_join(Frecuencia_por_habitat %>%
              group_by(TAXON_REF_PYR_MOD) %>%
              filter(n == max(n)) %>% 
              filter(row_number()==1) %>%
              distinct(), 'TAXON_REF_PYR_MOD') %>%
  rename(Species = TAXON_REF_PYR_MOD, Family = FAMILIA_WFO,
         RGR = RA_max, HS = B_h, `Main habitat` = EUNIS_LEVEL_2_GROUPED) %>% 
  dplyr::select(-n) %>% 
  dplyr::select(c(Species, Family, Order, Endemic, RGR, HS, LA, Rarity, `Main habitat`))
a %>% view

a %>%
  dplyr::select(c(Species, Rarity)) %>%
  mutate(Value = Rarity) %>%
  pivot_wider(values_from = Value, id_cols = Species, values_fill = NA, names_from = Rarity) %>%
  replace(is.na(.), '') %>%
  mutate(Rarity = paste(Endemic, RGR, HS, LA)) %>%
  mutate(Rarity = str_trim(Rarity)) %>%
  dplyr::select(c(Species, Rarity)) %>%
  left_join(a %>% dplyr::select(-Rarity), 'Species') %>%
  dplyr::select(c(Species, Family, Order, Rarity, Endemic, RGR, HS, LA, Rarity, `Main habitat`)) %>%
  distinct() %>%
  xlsx::write.xlsx(file = 'LIPA_significant.xlsx')

load('C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Modelos_taxo.RData')

ranef_families<-data.frame('Family' = rownames(ranef(m)$`FAMILIA_WFO:Order`),
                           'Beta Endemic' = unlist(ranef(m)$`FAMILIA_WFO:Order`), 
                           'Beta RGR' = unlist(ranef(m1)$`FAMILIA_WFO:Order`), 
                           'Beta HS' = unlist(ranef(m2)$`FAMILIA_WFO:Order`), 
                           'Beta LA' = unlist(ranef(m3)$`FAMILIA_WFO:Order`))
ranef_families$Family<-unlist(lapply(strsplit(ranef_families$Family, ':'), function(x) x[1]))


eiv %>%
  group_by(FAMILIA_WFO) %>%
  summarise(RGR = round(mean(as.numeric(RA_max)), 2), RGR_sd = round(sd(as.numeric(RA_max)), 2),
            HS = round(mean(B_h), 2), HS_sd = round(sd(B_h), 2),
            LAbu = round(mean(LA), 2), LA_sd = round(sd(LA), 2),
            Endemic = sum(Endemic == 'End'),
            N = n()) %>%
  mutate(RGR_ = paste(RGR, ' ', '(', RGR_sd, ')', sep = ''),
         HS_ = paste(HS, ' ', '(', HS_sd, ')', sep = ''),
         LA_ = paste(LAbu, ' ', '(', LA_sd, ')', sep = '')) %>%
  rename(Family = FAMILIA_WFO) %>%
  left_join(ranef_families, 'Family') %>%
  dplyr::select(c(Family, N, Endemic, RGR_, HS_, LA_, Beta.Endemic, Beta.RGR, Beta.HS, Beta.LA)) %>%
  rename(RGR = RGR_, HS = HS_, LA = LA_) %>%
  xlsx::write.xlsx(file = '~/Familia_mean.xlsx')
  