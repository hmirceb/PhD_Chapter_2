freqs<-as.data.frame(table(SIVIM$TAXON_REF_PYR_MOD))
freqs<-freqs[freqs$Freq >= 10,]
freqs<-freqs[sample(1:dim(freqs)[1], 500),]

runs<-10
sites<-10
theta<-list()
temp<-list()
SIVIM$pres<-1
for (i in 1:length(freqs$Var1)) {
  for (q in 1:runs) {
    invs.temp<-sample(unique(SIVIM[SIVIM$TAXON_REF_PYR_MOD == freqs$Var1[i],]$INV_CODE),  sites)
    SIVIM_temp<-SIVIM[SIVIM$INV_CODE %in% invs.temp,] %>% 
      select(INV_CODE, TAXON_REF_PYR_MOD, pres) %>%
      pivot_wider(names_from = INV_CODE, 
                  values_from = pres, values_fill = 0)
    temp[[q]]<-beta.multi(SIVIM_temp[,-1])$beta.SIM
  }
  theta[[i]]<-do.call('mean', temp)
  
}
aa<-data.frame('TAXON_REF_PYR_MOD' = freqs$Var1, 
           'theta' = unlist(theta))
