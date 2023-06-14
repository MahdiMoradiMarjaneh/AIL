# Interval map for Mouse QTL analysis
# AIL design for normally distrubuted trait
# Peter Thomson, University of Sydney
###################################################################

pset.AIL<-"MapAIL_pset.R"
fitq.AIL<-"MapAIL_fitq.R"

###############################################################
# Input data required:                                        #
# y = trait value                                             #
# dmark = map positions of markers (Morgan)                   #
# M = matrix of genotypes = no. of non-QSi5 alleles (0, 1, 2) #
# Q, M refer to QSi5, q. m refer to non-QSi5                  #
###############################################################

Traits<-read.csv("Phenotype file.csv")
Markers<-read.csv("MMU1 genotype.csv")

dmark1<-c(5.4549,7.704,10.308,11.3472,13.337,15.457,18.3493,20.7457,23.626,25.5443,26.6713,
   29.333,32.307,33.9014,37.996,39.511,41.6201,43.9953,45.4767,47.764,49.258,50.344,51.606,
   53.1661,55.798,58.0582,60.524,62.3284,63.32,64.7135)/100 # positions as Morgan
dmark<-dmark1

n.AIL<-14  # Number of AIL generations
epsd<-10^-6
ptrace<-0
Dmark<-c(dmark,-Inf,Inf)
delta.d<-0.0025 # Fit model every 0.25 cM
dmin<-floor(100*min(dmark))/100
dmax<-ceiling(100*max(dmark))/100
d<-dmin

# merge trait and genotype file: include all phenotypes
Merged<-merge(Traits,Markers,by="Mouse.ID",sort=T,all.x=T)

y<-Merged$FOW ## update this
# y<-Merged$FVL ## update this
ybar<-mean(y,na.rm=T)
y[is.na(y)]<-ybar
Age<-Merged$Age # get copies of variables to be adjusted for
Sex<-Merged$Sex
Weight<-Merged$Weight
n<-length(y)

# Recode to M = #QSi5 alleles (0,1,2)
Mtemp<-Merged[,-(1:21)] # Delete columns that are not markers
M<-1*Mtemp
M<-ifelse(Mtemp=="",NA,(Mtemp=="Q")*2+(Mtemp=="H")*1+(Mtemp=="S")*0)
nmark<-ncol(M)
MM<-cbind(M,0,0)

# Initial fit
y.fit<-lm(y~1) ## Update (no adjustment)
# y.fit<-lm(y~Sex) ## Update (adjust for Sex only)
# y.fit<-lm(y~Sex+Age) ## Update (adjust for Sex and Age)
g<-c(0,0)
b.lm<-c(coef(y.fit),g)
s.lm<-summary(y.fit)$sigma

# Duplication for mixture

q1<-rep(1,n)
q0<-rep(0,n)
Q.star<-rbind(cbind(q0,q0),cbind(q1,q0),cbind(q0,q1))
# Reoplicate variabkes that need to be adjusted for
Sex.star<-rep(Sex,3)
Age.star<-rep(Age,3)
Weight.star<-rep(Weight,3)
y.star<-rep(y,3)

# Interval map

dset<-as.numeric(names(table(c(seq(dmin,dmax,delta.d),dmark))))
map.fit<-NULL

for (i in 1:length(dset))
{
   d<-dset[i]

#  Determine flanking markers at a map position d for each animal

#  list of observed markers (1,...,nmark) if observed, NA otherwise
   omark<-matrix(rep(1:nmark,n),byrow=T,nrow=n)*(M*0+1)

   m10<-max((1:nmark)[d>=dmark])
   m10<-ifelse(is.na(m10),-Inf,m10)
   omarkl<-ifelse(omark<=m10&!is.na(omark),omark,-Inf)
   m1<-apply(omarkl,1,max)
   m1<-ifelse(m1==-Inf,nmark+1,m1)

   m20<-min((1:nmark)[d<=dmark])
   m20<-ifelse(is.na(m20),Inf,m20)
   omarku<-ifelse(omark>=m20&!is.na(omark),omark,Inf)
   m2<-apply(omarku,1,min)
   m2<-ifelse(m2==Inf,nmark+2,m2)

   d1<-Dmark[m1]
   d2<-Dmark[m2]
   r12<-0.5*(1-exp(-2*(d2-d1)))
   r12<-(1-(1-r12)^(n.AIL-2)*(1-2*r12))/2  # D&S (1005) AIL adjustment

#  Determine set of conditional QTL probabilities given marker information

   source(pset.AIL)

#  Flanking marker group	Flanking marker genotypes
#  1 = M1M1 M2M2		00
#  2 = M1M1 M2m2		01
#  3 = M1M1 m2m2		02
#  4 = M1m1 M2M2		10
#  5 = M1m1 M2m2		11
#  6 = M1m1 m2m2		12
#  7 = m1m1 M2M2		20
#  8 = m1m1 M2m2		21
#  9 = m1m1 m2m2		22

   mg<-3*diag(MM[,m1])+diag(MM[,m2])+1
#  p = Pr{QQ | M)
   p.QQ<-diag(pset.QQ[,mg])
   p.Qq<-diag(pset.Qq[,mg])
   p.qq<-diag(pset.qq[,mg])
   source(fitq.AIL) # Update MapAIL.R foe terms to adjyust for
   print(c(d,b.lm,s.lm,logL,iter))
   map.fit<-rbind(map.fit,c(d,b.lm,s.lm,logL,iter))

}

# fit a no QTL model = initial fit

y.fit<-lm(y~1) ## update this
# y.fit<-lm(y~Sex) ## update this
# y.fit<-lm(y~Sex+Age) ## update this
b.lm<-coef(y.fit)
s.lm<-summary(y.fit)$sigma
mu<-matrix(1,n,1) %*% b.lm # update this
# mu<-model.matrix(~Sex) %*% b.lm
# mu<-model.matrix(~Sex+Age) %*% b.lm
f<-dnorm(y,mu,s.lm)
logL<-sum(log(f))
map.fit<-rbind(map.fit,c(NA,b.lm,0,0,s.lm,logL,1))
LRT.col <- 6 # column number with LRT, 6: no adj, 7: 1 adj, 8: 2 ad.j etc
map.fit[,LRT.col]<-2*(map.fit[,LRT.col]-logL) # 2 adjustments (Sex+Age)

colnames(map.fit)<-c("d","Intercept","Qq","QQ","sd","LRT","Iter")
# colnames(map.fit)<-c("d","Intercept","SexM","Qq","QQ","sd","LRT","Iter")
# colnames(map.fit)<-c("d","Intercept","SexM","Age","Qq","QQ","sd","LRT","Iter")

# Graphical output of interval map

plot(100*map.fit[,1],map.fit[,LRT.col]/(2*log(10)),type="l",xlab="Map Position (cM)",
     ylab="LOD Score")
abline(v=100*dmark,lty=2)

write.csv(map.fit,"MMU1_FOW_no adj.csv",row.names=F) ## update this

