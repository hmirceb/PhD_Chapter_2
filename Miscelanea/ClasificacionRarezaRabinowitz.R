florapyr<-read.csv('C:/Users/18172844S/Dropbox/Tesis/Datos/FLORAPYR/FLORAPYR.csv', sep = ';')

# # 
florapyr<-florapyr[florapyr$AUTOCT_SSP == 'x',]

# # Eliminamos las especies que no tienen informacion sobre el tamaño poblacional
florapyr<-florapyr[florapyr$POPSIZE != '',]

# # Vamos a crear tres variables, una para el rango de distribucion, otra para especificidad de habitat y otra para el tamaño poblacional
rabinowitz<-data.frame(florapyr$TAXON_PYR_SP)

# # Rango de distribucion. Amplio (Wide) o estrecho (Narrow). Se basa en la informacion sobre la corologia de cada especie. Consideramos que todas aquellas especies
# # no clasificadas como Cosmopolitas o de amplia distribucion (AD) se considera que tienen un rango Amplio
rabinowitz$range<-as.factor(ifelse(florapyr$CORO_PYR_SP == 'AD', 'Wide', 'Narrow'))

# # Especificidad de habitat. Amplio (Broad) o restringido (Restricted). En funcion de los requerimientos edaficos de cada planta. Las
# # plantas indiferentes al sustrato se consideran de rango Amplio.
rabinowitz$specific<-as.factor(ifelse(florapyr$EDAF_SP == 'Indiferent', 'Broad', 'Restricted'))

# # Tamaño poblacional. Pequeño (Small) o grande (Large). Se juntan las categorias Medium y Large como una sola. Unicamente las plantas
# # con poblaciones pequeñas se consideran pequeñas.
rabinowitz$popsize<-as.factor(ifelse(florapyr$POPSIZE == 'Small', 'Small', 'Large'))

# # Clasificamos en una de las 8 categorias de rareza de Rabinowitz basandonos en las variables anteriores
rabinowitz$rareness<-as.factor(paste(substr(rabinowitz$range, 1, 1), substr(rabinowitz$specific, 1, 1), substr(rabinowitz$popsize, 1, 1)))

table(rabinowitz$rareness)/dim(rabinowitz)[1]
table(rabinowitz$range)/dim(rabinowitz)[1]
table(rabinowitz$specific)/dim(rabinowitz)[1]
table(rabinowitz$popsize)/dim(rabinowitz)[1]
