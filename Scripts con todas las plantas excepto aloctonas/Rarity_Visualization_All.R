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
library(insight)
library(rstanarm)
library(gridExtra)
library(phytools)

load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Valores_Rareza.RData")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/FiloSignal_Rareza_scaled_All.RData")
load('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/FiloSignal_Rareza_scaled_All.RData')
load("C:/Users/18172844S/Dropbox/DATA__LAB/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
load("~/Dropbox/DATA__LAB/PHYLO/de_Hector/ArbolesPirineo_10_12_2020.RData")
load("C:/Users/18172844S/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")
load("~/Dropbox/Tesis/Angiospermas FLORAPYR/Resultados/FLORAPYR_Grupos_WFO.RData")

rm(tab_v2, tre, w1, 
   RA_node_all, RA_node_ang, RA_node_gimn,
   End_node_all, End_node_ang, End_node_gimn,
   LA_node_all, LA_node_ang, LA_node_gimn,
   HS_node_all, HS_node_ang, HS_node_gimn)

data<-eiv %>% 
  dplyr::select(c(TAXON_REF_PYR_MOD, RA_max, LA, B_h, Endemic, CAT_AMENAZA_PYR)) %>%
  mutate(CAT_AMENAZA_PYR = ifelse(is.na(CAT_AMENAZA_PYR), 'No',
                                  CAT_AMENAZA_PYR)) %>%
  dplyr::rename(Red_list = CAT_AMENAZA_PYR) %>%
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
#### Figure 2 ####
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
  theme(legend.position = 'none',
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

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

ggsave(p3, filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Figure_2_Tree.jpeg', 
       height = 29, width = 29, units = 'cm')
ggsave(p3, filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Figure_2_Tree.jpeg', 
       height = 29, width = 29, units = 'cm')


#### Figure 1 ####
# Boxplot de la senyal filo ####
cols_boxplot<-c('Endemic' = "#FF8C00BF", 'RGR' = 'firebrick', 'HS' = 'steelblue', 'LA' = 'forestgreen')

signal<-do.call('rbind', lapply(signal_all, function(x) 
  data.frame('rarity' = rownames(x$stat),
             'lambda' = x$stat$Lambda, 
             'p_val' = x$pvalue$Lambda, 
             'group' = 'All')))
dat_signal<-signal %>%
  mutate(rarity=dplyr::recode(rarity, 
                       RA_max = 'RGR',
                       B_h = 'HS')) %>%
  mutate(rarity = factor(rarity, levels = c('Endemic', 'RGR', 'HS', 'LA'))) %>%
  dplyr::rename(Rarity = rarity)
dat_text<-dat_signal %>%
  aggregate(p_val~Rarity+group, ., mean) %>%
  mutate(yloc = max(dat_signal$lambda))

boxplot_signal<-ggplot(dat_signal, aes(x = Rarity, y = lambda, fill = Rarity))+geom_boxplot()+
  theme_classic()+ylab('Phylogenetic signal')+xlab('Rarity type')+
  theme(text = element_text(size = 20),
        strip.text = element_text(face = 'bold'),
        legend.position = 'none')+
  scale_fill_manual(values = cols_boxplot, name = 'Rarity type')+
  ggtitle('(a)')

boxplot_signal

#### PD loss ####
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/PD_loss_all.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/PD_loss_all.RData")

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
  mutate(Rarity = factor(Rarity, levels = c('Endemic', 'RGR', 'HS', 'LA'))) %>%
  ggplot(aes(x = Rarity))+
  geom_segment(aes(yend = ses, y = 0, xend = Rarity),
               linetype = 'dashed', size = 0.7)+
  theme_classic()+
  geom_point(aes(y = ses, color = Rarity), size = 3)+
  scale_color_manual(values = c('Endemic' = '#FF8C00', 
                                'RGR' = 'firebrick',
                                'HS' = 'steelblue',
                                'LA' = 'forestgreen'), 
                     name = "Rarity type:", guide = 'none')+
  theme(text = element_text(size = 20))+
  ylab('Standard effect size of PD loss')+
  xlab('Rarity type')+
  annotate(geom = 'rect', ymin = -1.96, 
           ymax = 1.96, xmin = 0, xmax = 5,
           alpha = 0.1)+
  ggtitle('(b)')

ggsave(grid.arrange(boxplot_signal, ses_loss_plot, ncol = 2), 
       filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Figure_1_Boxplot_And_Loss.jpeg',
       scale = 1.5, dpi = 300)

ggsave(grid.arrange(boxplot_signal, ses_loss_plot, ncol = 2), 
       filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Figure_1_Boxplot_And_Loss.jpeg',
       scale = 1.5, dpi = 300)

#### Figure 3 ####
#### Variation partition ####
load('C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Modelos_taxo_all.RData')
load('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Resultados/Modelos_taxo_all.RData')

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
var_prop_data<-dat_var %>%
  mutate(Rarity = rep(c('Endemic', 'RGR', 'HS', 'LA'), each = 4)) %>%
  mutate(Level = rep(c('Residual', 'Genus', 'Family', 'Order'), times = 4)) %>%
  group_by(Rarity) %>%
  dplyr::group_by(Rarity) %>%
  dplyr::summarise(Var = Var, 
                   Var_prop = Var/sum(Var),
                   Rarity = Rarity,
                   Level = Level) %>%
  ungroup() %>%
  mutate(Rarity = factor(Rarity, levels = c('LA', 'HS', 'RGR', 'Endemic'))) %>%
  mutate(Level = factor(Level, levels = c('Residual', 'Genus', 'Family', 'Order')))


var_prop<-ggplot(var_prop_data, 
                 aes(x = Rarity, 
                     y = 100*Var_prop, 
                     fill = Level))+
  geom_bar(stat = 'identity', color = 'black')+
  theme_classic()+
  ylab('Variation explained (%)')+
  xlab('Rarity type')+
  ggtitle('(a)')+
  theme(text = element_text(size = 20),
        legend.position = 'bottom')+
  scale_fill_manual(values=c( 'gray90', "#009E73", "#D55E00", "mediumorchid4"))

### Table outputs ####
#a<-eiv %>%
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
  dplyr::select(-c(LIPA_Endemic, LIPA_RA_max,
                   LIPA_LA, LIPA_B_h, LIPA_Endemic_pval, LIPA_RA_max_pval,
                   LIPA_LA_pval,LIPA_B_h_pval)) %>%
  left_join(Frecuencia_por_habitat %>%
              group_by(TAXON_REF_PYR_MOD) %>%
              filter(n == max(n)) %>% 
              filter(row_number()==1) %>%
              distinct(), 'TAXON_REF_PYR_MOD') %>%
  dplyr::rename(Species = TAXON_REF_PYR_MOD, Family = FAMILIA_WFO,
         RGR = RA_max, HS = B_h, `Main habitat` = EUNIS_LEVEL_2_GROUPED) %>% 
  dplyr::select(-n) %>% 
  dplyr::select(c(Species, Family, Order, Endemic, RGR, HS, LA, Rarity, Red_list, `Main habitat`))
a %>% view

#a %>%
  dplyr::select(c(Species, Rarity)) %>%
  mutate(Value = Rarity) %>%
  pivot_wider(values_from = Value, id_cols = Species, values_fill = NA, names_from = Rarity) %>%
  replace(is.na(.), '') %>%
  mutate(Rarity = paste(Endemic, RGR, HS, LA)) %>%
  mutate(Rarity = str_trim(Rarity)) %>%
  dplyr::select(c(Species, Rarity)) %>%
  left_join(a %>% dplyr::select(-Rarity), 'Species') %>%
  dplyr::rename(`Red list status` = Red_list) %>%
  dplyr::select(c(Species, Family, Order, Rarity, `Red list status`, Endemic, RGR, HS, LA, Rarity, `Main habitat`)) %>%
  distinct() %>%
  xlsx::write.xlsx(file = 'LIPA_significant.xlsx')


#eiv %>%
  group_by(FAMILIA_WFO) %>%
  summarise(RGR = round(mean(as.numeric(RA_max)), 2), RGR_sd = round(sd(as.numeric(RA_max)), 2),
            HS = round(mean(B_h), 2), HS_sd = round(sd(B_h), 2),
            LAbu = round(mean(LA), 2), LA_sd = round(sd(LA), 2),
            Endemic = sum(Endemic == 'End'),
            N = n()) %>%
  mutate(RGR_ = paste(RGR, ' ', '(', RGR_sd, ')', sep = ''),
         HS_ = paste(HS, ' ', '(', HS_sd, ')', sep = ''),
         LA_ = paste(LAbu, ' ', '(', LA_sd, ')', sep = '')) %>%
  dplyr::rename(Family = FAMILIA_WFO) %>%
  left_join(ranef_families, 'Family') %>%
  dplyr::select(c(Family, N, Endemic, RGR_, HS_, LA_, Beta.Endemic, Beta.RGR, Beta.HS, Beta.LA)) %>%
  dplyr::rename(RGR = RGR_, HS = HS_, LA = LA_) %>%
  xlsx::write.xlsx(file = '~/Familia_mean.xlsx')


#### PCA filo ####

# funcion para las etiquetas circulares
rshift = function(r, theta, a=0.1, b=0.2) { 
  r + a + b*abs(cos(theta))}

pca_list_loadings<-(mean.list(lapply(pca_filo_trees, function(x) abs(x$Evec))))*
  (pca_filo_trees[[3]]$Evec/abs(pca_filo_trees[[3]]$Evec))

pca_list_loadings<-as.data.frame(pca_list_loadings)
pca_list_loadings$Rarity<-rownames(pca_list_loadings)
pca_list_loadings<-pca_list_loadings %>%
  mutate(r = sqrt(PC1^2 + PC2^2),
         theta = atan2(PC2,PC1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta)) %>%
  mutate(Rarity = dplyr::recode(Rarity, 'RA_max' = 'RGR', 'B_h' = 'HS'))

lambda_pca<-data.frame('lamba' = mean(unlist(lapply(pca_filo_trees, function(x) x$lambda))),
                       'sd' =  sd(unlist(lapply(pca_filo_trees, function(x) x$lambda))))

pca_list_scores<-(mean.list(lapply(pca_filo_trees, function(x) abs(x$S)))*
                    (pca_filo_trees[[3]]$S/abs(pca_filo_trees[[3]]$S)))
pca_list_scores<-as.data.frame(pca_list_scores)
pca_list_scores$TAXON_REF_PYR_MOD<-rownames(pca_list_scores)

pca_list_scores<-pca_list_scores %>%
  mutate(Endemic = ifelse(TAXON_REF_PYR_MOD %in% 
                            data_all[data_all$Endemic == '1',]$TAXON_REF_PYR_MOD, 
                          'End', 'NoEnd')) %>%
  mutate(across(c(PC1, PC2), as.numeric))  %>%
  left_join(., data_all %>% 
              dplyr::select(c(TAXON_REF_PYR_MOD,
                              Red_list)), 'TAXON_REF_PYR_MOD') %>%
  mutate(Red_list_simple = dplyr::recode(Red_list, 'NT' = 'No', 
                                         'CR' = 'Si',
                                         'VU' = 'Si',
                                         'EN' = 'Si')) %>%
  mutate(Red_list_simple = dplyr::recode(Red_list_simple, 'Si' = 'Threatened')) %>%
  mutate(Color = paste(Endemic, Red_list_simple))

# Proporcion de varianza explicada por cada componente (sale de los eigenvalues (Eval))
pc_variance<-mean.list(
  lapply(pca_filo_trees,
         function(x) as.matrix(100*round(diag(x$Eval)/sum(diag(x$Eval)), 3))))

# Funcion para calcular el rango sin signo
unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)), 
                                abs(max(x, na.rm = TRUE)))

# Ratio de conversion entre las flechas de los eigenvectors y los scores (para que las flechas queden bien)
ratio<-max(unsigned.range(pca_list_loadings$PC1)/unsigned.range(pca_list_scores$PC1), 
           unsigned.range(pca_list_loadings$PC2)/unsigned.range(pca_list_scores$PC2))
phylo_pca<-pca_list_scores %>%
  arrange(Red_list_simple) %>%
  ggplot(aes(x = PC1, y = PC2))+
  geom_vline(xintercept = 0, alpha = 0.5, linetype = 'dashed')+
  geom_hline(yintercept = 0, alpha = 0.5, linetype = 'dashed')+
  geom_point(aes(alpha = ifelse(Red_list_simple != 'No' |
                                  Endemic == 'End', 1, 0.3),
                 fill = Red_list_simple,
                 color = Color,
                 shape = Endemic,
                 size = Color))+
  geom_segment(data = pca_list_loadings,
               inherit.aes = F, 
               aes(x = 0, y = 0, xend = (0.9*PC1)/ratio, yend = (0.9*PC2)/ratio), #convertimos por el ratio
               arrow = arrow(length = unit(0.2, "cm")),
               color = '#C30E2D')+
  geom_text(data = pca_list_loadings,
            inherit.aes = F, 
            aes(x = PC1/ratio, y = PC2/ratio, label = Rarity),
            color = '#C30E2D', 
            fontface = 'bold')+
  theme_classic()+
  scale_fill_manual(values = c('No' = 'black',
                               'Threatened' = 'red'),
                    name = 'Red List: ')+
  scale_color_manual(values = c('End No' = 'white',
                                'End Threatened' = 'black',
                                'NoEnd No' = 'black',
                                'NoEnd Threatened' = 'black'),
                     name = 'Red List: ',
                     guide = 'none')+
  scale_size_manual(values = c('End No' = 4,
                                'End Threatened' = 4,
                                'NoEnd No' = 3,
                                'NoEnd Threatened' = 4),
                     name = 'Red List: ',
                     guide = 'none')+
  scale_shape_manual(values = c('NoEnd' = 21, 'End' = 22), 
                     guide = 'none')+
  scale_alpha_continuous(guide = 'none')+
  ggtitle('(b)')+
  theme(text = element_text(size = 20),
        legend.position = 'bottom')+
  guides(color = 'none',
         size = 'none',
         fill = guide_legend(override.aes = list(size = 7,
                                                 shape = 21, 
                                                 color = 'black')))+
  coord_fixed(xlim = range(pca_list_scores$PC1), ylim = range(pca_list_scores$PC2))+
  xlab(paste('PC1 ', '(', pc_variance[1], '% explained variance)', sep = ''))+
  ylab(paste('PC2 ', '(', pc_variance[2], '% explained variance)', sep = ''))
phylo_pca


### Figure 3 ####
Multipanel<-grid.arrange(var_prop+coord_flip(),
                         phylo_pca,
                         nrow = 1, ncol = 2)
ggsave(Multipanel,
       filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Fig_3_Partition_PCA.jpeg',
       dpi = 300,
       scale = 1.8)
ggsave(Multipanel,
       filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Fig_3_Partition_PCA.jpeg', 
       dpi = 300,
       scale = 1.8)

#### Supp mat Figure 1 ####
cols_boxplot<-c('Endemic' = "#FF8C00BF", 'RGR' = 'firebrick', 'HS' = 'steelblue', 'LA' = 'forestgreen')

signal_supp<-rbind(do.call('rbind', lapply(signal_angios, function(x) 
  data.frame('rarity' = rownames(x$stat),
             'lambda' = x$stat$Lambda, 
             'p_val' = x$pvalue$Lambda, 
             'group' = 'Angiosperms'))),
  do.call('rbind', lapply(signal_gimn, function(x) 
    data.frame('rarity' = rownames(x$stat),
               'lambda' = x$stat$Lambda, 
               'p_val' = x$pvalue$Lambda, 
               'group' = 'Gimnosperms & Pteridophytes'))))

dat_signal_supp<-signal_supp %>%
  mutate(rarity=dplyr::recode(rarity, 
                              RA_max = 'RGR',
                              B_h = 'HS')) %>%
  mutate(rarity = factor(rarity, levels = c('LA', 'HS', 'RGR', 'Endemic'))) %>%
  dplyr::rename(Rarity = rarity) %>%
  filter(!is.na(p_val))
dat_text<-dat_signal_supp %>%
  aggregate(p_val~Rarity+group, ., mean) %>%
  mutate(yloc = max(dat_signal_supp$lambda))

boxplot_signal_supp<-ggplot(dat_signal_supp, aes(x = Rarity, y = lambda, fill = Rarity))+
  geom_boxplot()+
  theme_classic()+
  ylab('Signal')+
  xlab('Rarity type')+
  geom_text(data = dat_text, aes(x = Rarity, y = yloc,
                                 label = ifelse(p_val < 0.05,'***', '')),
            inherit.aes = F, size = 8)+
  theme(text = element_text(size = 20),
        strip.text = element_text(face = 'bold'),
        legend.position = 'none')+
  facet_wrap(~group)+
  scale_fill_manual(values = cols_boxplot, name = 'Rarity type')

ggsave(boxplot_signal_supp, scale = 1.5,
       filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Figure_S1_Boxplot_signal.jpeg')

ggsave(boxplot_signal_supp, scale = 1.5,
       filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Figure_S1_Boxplot_signal.jpeg')


#### TIFFs ####
ggsave(grid.arrange(boxplot_signal, ses_loss_plot, ncol = 2), 
       filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Figure_1_Boxplot_And_Loss.tiff',
       scale = 1.5, dpi = 300, compression = 'jpeg')

ggsave(grid.arrange(boxplot_signal, ses_loss_plot, ncol = 2), 
       filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Figure_1_Boxplot_And_Loss.tiff',
       scale = 1.5, dpi = 300, compression = 'jpeg')


ggsave(p3, filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Figure_2_Tree.tiff', 
       height = 29, width = 29, units = 'cm',
       compression = 'jpeg')
ggsave(p3, filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Figure_2_Tree.tiff', 
       height = 29, width = 29, units = 'cm',
       compression = 'jpeg')


ggsave(Multipanel,
       filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Fig_3_Partition_PCA.tiff',
       dpi = 300,
       scale = 1.8, compression = 'jpeg')
ggsave(Multipanel,
       filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 2 - Senal Filogenetica en la rareza del Pirineo/Figures and tables/Fig_3_Partition_PCA.tiff', 
       dpi = 300,
       scale = 1.8, compression = 'jpeg')