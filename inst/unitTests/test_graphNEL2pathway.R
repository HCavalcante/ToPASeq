
test_graphNEL2pathway<-function(){
g1<-pathways("hsapiens","kegg")[[1]]


g2<-as(g1,"graphNEL")
g3<-graphNEL2Pathway(g2)

checkTrue(length(nodes(g1))==length(nodes(g3)))
checkTrue(nrow(edges(g1))==nrow(edges(g3)))

set.seed(1237)
V<-letters[1:10]
from <- sample(V,7)
 to <- sample(V,7)
 ft <- cbind(from, to)
 ft <- rbind(ft,c("a","b"))
 g <- ftM2graphNEL(ft)
gg<-graphNEL2Pathway(g)

checkTrue(length(nodes(g))==length(nodes(gg)))
checkTrue(length(unlist(edges(g)))==nrow(edges(gg)))
checkTrue(nlevels(edges(gg)[,4])==1)

pathways<-pathways("hsapiens","kegg")[[1]]

G2<-graphNEL2Pathway(pathwayGraph(pathways))
checkTrue(identical(sort(nodes(pathways)),sort( nodes(G2)) ))
E.orig<-edges(pathways)
E.orig<-E.orig[order(E.orig[,1], E.orig[,2], E.orig[,4]),]
rownames(E.orig)<-seq_len(nrow(E.orig))
checkTrue(identical(E.orig[,1:2], edges(G2)[,1:2]))
#checkTrue(identical(as.character(E.orig[,3]), as.character(edges(G2)[,3])))
checkTrue(identical(as.character(E.orig[,4]), as.character(edges(G2)[,4])))

set.seed(123)
rg <- randomEGraph(LETTERS[1:20], edges = 30)

p<-graphNEL2Pathway(rg)
checkEquals(length(names(which(degree(rg)>0))), length(nodes(p)))
checkTrue(identical( edgeNames(rg), paste(edges(p)[,1], edges(p)[,2], sep="~") ))


}
