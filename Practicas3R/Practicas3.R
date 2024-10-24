#cargamos las librerias

library(ape)
library(phangorn)
library(phytools)
#cambiamos el archivo fasta a otro compatible con phanrgon
fraxatin <- read.phyDat(file = "fraxatin_aligned.fasta", 
                       format = "FASTA", type = "AA")
#cambiamos el objeto a .AAbin
matrizdist <- as.AAbin(fraxatin)
#calcula la matriz de distancias por pares
matrizdist <- dist.aa(matrizdist)

matrizdist

#Arboles de distancias
#metodo UPGMA
  arbolUPGMA <- upgma(matrizdist)
  plot(arbolUPGMA, type= "fan", cex=0.7, edge.width=2, edge.color="green", font=3)
#método NJ, unión de vecinos
  arbolNJ <- nj(matrizdist)
  plot(arbolNJ, type= "p", edge.lty=5, cex=0.8, edge.width=3, edge.color="purple", font=3)#probando distintos tipos de plots;)
  plot(arbolNJ, type= "c", cex=0.8, edge.width=3, edge.color="purple", font=3, bg="black")
#hacer gráfico con phytools
  plotTree(arbolNJ)
  plotTree(arbolNJ, ftype="b", fsize=0.8, offset=1, color="gold", lwd=2)
  plotTree(ladderize(arbolNJ))#sirve para cambiar el orden de los grupos, en este caso por orden alfabetico
  
#guardamos el tree
  write.tree(arbolNJ, file = "arbolNJ.nex")
  read.tree(file = "arbolNJ.nex") #leemos el archivo
  
#para enrraizar el arbol
  arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
  plot(arbolNJraiz)
  arbolUPGMAraiz <-root(arbolUPGMA, outgroup = "Ornitorrinco", r=TRUE)
  plot(arbolUPGMAraiz)  
  
#ver dos trees a la vez
  layout(matrix(c(1,2)), height=c(10,10))
  par(mar=c(1,1,1,1))
  plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex=0.4)
  plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex=0.4)
  
  
  
#PARSIMONIA
  parsimony(arbolUPGMAraiz, fraxatin)#Cuenta el número de pasos del arbol ([1] 313)
  parsimony(arbolUPGMA, fraxatin)#en el de sin raiz es lo mismo
  
  mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin)#así se obtiene el arbol con paxima parsimonia(menos pasos)
  mejorNJ <- optim.parsimony(arbolNJraiz, fraxatin)
  #CUIDADO, mas de 10 taxones necesita unos requerimientos computacionales muy altos!!
  
  #para poder hacer árboles con maxima parsimonia sin necesitar muchos requerimientos, podemos hacer que en vez de leer todo el genoma vea solo un porcentaje
  fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE)
  fraxatin_parsimonia  #4 árboles filogenéticos
  #los enrraizamos
  fraxatin_parsimoniaR <- root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
  plot(fraxatin_parsimoniaR, cex = 0.6, color="brown")#salen los 4 árboles
  
  #creamos el arbol consenso
  estrictode100 <- consensus(fraxatin_parsimoniaR, p = 1)
  plot(estrictode100, cex = .6)
  estrictode30 <- consensus(fraxatin_parsimoniaR, p = 0.3)#podemos cambiar los parametros segun como de estrictos queramos ser
  plot(estrictode30, cex = .6)


#BOOTSTRAP
  #El bootstrap permite crear seudoreplicas, con cada seudoréplica creamos un arbol consenso, proporcionando valores de confianza por cada clado de un arbol según la proporción de ese mismo clado.
  arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs = 10)#la función prachet hace que tarde menos
  #generamos un consenso
  estricto60 <- consensus(arbolesbootstrap, p = 0.6)
  plot(estricto60, cex = .6)
  
  
#MÉTODOS PROBABILISTICOS, modelo evolutivo
  #máxima verosimilitud, probabilidad de obtener un dato según un modelo, 3 PASOS
  arbolazar <- rtree(n = 11, tip.label = names(fraxatin))#como necesita muchos recursos, creamos un arbol de 11 ramas aleatorio como base
  plot(arbolazar, cex = .5)
  #lo enraizamos y lo escalamos
  arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
  plot(ladderize(arbolazarR), cex = .5); add.scale.bar()
  #Calculamos la verosimilitud con pml
  ajustado <- pml(arbolazarR, fraxatin)
  ajustado
  #model: Mk loglikelihood: -4080.529  unconstrained loglikelihood: -1479.871 Rate matrix:
  
  #Ahora vamos a encontrar un árbol que optimice la verosimilitud usando la matriz de sustitución de Dayhoff
  ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")
  
  ajustadoconDay$tree#para ver los resultados del árbol oculto
  ajustadoconDayraíz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")#lo enraizamos
  plot(ladderize(ajustadoconDayraíz), cex = .5); add.scale.bar()
  
  #podemos hacer lo mismo con otros modelos de sustitución distintos (blosum62)
  ajustadoconBlo <- optim.pml(object = ajustado, model = "Blosum62", rearrangement = "ratchet")
  #JTT
  ajustadoconJTT <- optim.pml(object = ajustado, model = "JTT", rearrangement = "ratchet")  
  
  #ahora podemos comparar los modelos
  AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)
  ##                df (grados de libertad)     AIC ( más bajo, mejor)
  ## ajustadoconDay 19                          5188.886 
  ## ajustadoconBlo 19                          5197.264
  ## ajustadoconJTT 19                          5061.388 (el mejor)
  
  #el mejor es el JTT, es el que utilizaremos
 mejorarbol <- optim.pml(object = ajustadoconDay, model = "JTT",  rearrangement = "ratchet")
 mejorarbol
 mejorarbolR <- root(mejorarbol$tree, outgroup = "Ornitorrinco")#enraizamos
 plot(ladderize(mejorarbolR), cex = 0.5); add.scale.bar()
 #ya tendriamos nuestro mejor árbol.