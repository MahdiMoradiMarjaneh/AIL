# Returns the prob of QTL genotype | Marker genotype for F2 design
# Routine called by MapAIL_main.R

r1<-0.5*(1-exp(-2*(d-d1)))
r2<-0.5*(1-exp(-2*(d2-d)))
r12<-0.5*(1-exp(-2*(d2-d1)))
# D&S (1005) AIL adjustments
r1<-(1-(1-r1)^(n.AIL-2)*(1-2*r1))/2
r2<-(1-(1-r2)^(n.AIL-2)*(1-2*r2))/2
r12<-(1-(1-r12)^(n.AIL-2)*(1-2*r12))/2

n1<-1-r1
n2<-1-r2

pset.QQ<-matrix(0,nrow=n,ncol=9)
pset.Qq<-matrix(0,nrow=n,ncol=9)
pset.qq<-matrix(0,nrow=n,ncol=9)

pset.QQ[,1]<-n1^2*n2^2     / (1-r12)^2
pset.Qq[,1]<-2*r1*r2*n1*n2 / (1-r12)^2
pset.qq[,1]<-r1^2*r2^2     / (1-r12)^2

pset.QQ[,2]<-2*n1^2*r2*n2              / (2*r12*(1-r12))
pset.Qq[,2]<-2*(r1*n1*n2^2+r1*n1*r2^2) / (2*r12*(1-r12))
pset.qq[,2]<-2*r1^2*r2*n2              / (2*r12*(1-r12))

pset.QQ[,3]<-n1^2*r2^2     / r12^2
pset.Qq[,3]<-2*r1*n1*r2*n2 / r12^2	
pset.qq[,3]<-r1^2*n2^2     / r12^2

pset.QQ[,4]<-2*r1*n1*n2^2              / (2*r12*(1-r12))
pset.Qq[,4]<-2*(n1^2*r2*n2+r1^2*r2*n2) / (2*r12*(1-r12))
pset.qq[,4]<-2*r1*n1*r2^2              / (2*r12*(1-r12))

pset.QQ[,5]<-4*r1*n1*r2*n2                               / (2*((1-r12)^2+r12^2))
pset.Qq[,5]<-2*(n1^2*n2^2+r1^2*r2^2+n1^2*r2^2+r1^2*n2^2) / (2*((1-r12)^2+r12^2))
pset.qq[,5]<-4*r1*n1*r2*n2                               / (2*((1-r12)^2+r12^2))

pset.QQ[,6]<-2*r1*n1*r2^2              / (2*r12*(1-r12))
pset.Qq[,6]<-2*(n1^2*r2*n2+r1^2*r2*n2) / (2*r12*(1-r12))
pset.qq[,6]<-2*r1*n1*n2^2              / (2*r12*(1-r12))

pset.QQ[,7]<-r1^2*n2^2     / r12^2
pset.Qq[,7]<-2*r1*n1*r2*n2 / r12^2
pset.qq[,7]<-n1^2*r2^2     / r12^2

pset.QQ[,8]<-2*r1^2*r2*n2              / (2*r12*(1-r12))
pset.Qq[,8]<-2*(r1*n1*n2^2+r1*n1*r2^2) / (2*r12*(1-r12))
pset.qq[,8]<-2*n1^2*r2*n2              / (2*r12*(1-r12))

pset.QQ[,9]<-r1^2*r2^2     / (1-r12)^2
pset.Qq[,9]<-2*r1*n1*r2*n2 / (1-r12)^2
pset.qq[,9]<-n1^2*n2^2     / (1-r12)^2
