# Returns the fitted model at a map position

conv1<-1
conv2<-1
logL<--10^99
iter<-1

while (conv1>10^-8 & conv2>10^-4)
{

#  linear regression for y

   b.lm.old<-b.lm
   logL.old<-logL
   mu.lm.star<-model.matrix(~Q.star) %*% b.lm  ## update this
   mu.QQ.lm<-mu.lm.star[1:n]
   mu.Qq.lm<-mu.lm.star[(n+1):(2*n)]
   mu.qq.lm<-mu.lm.star[(2*n+1):(3*n)]
   f.QQ<-dnorm(y,mu.QQ.lm,s.lm)
   f.Qq<-dnorm(y,mu.Qq.lm,s.lm)
   f.qq<-dnorm(y,mu.qq.lm,s.lm)

   tau.QQ<-p.QQ*f.QQ/(p.QQ*f.QQ+p.Qq*f.Qq+p.qq*f.qq)
   tau.Qq<-p.Qq*f.Qq/(p.QQ*f.QQ+p.Qq*f.Qq+p.qq*f.qq)
   tau.qq<-p.qq*f.qq/(p.QQ*f.QQ+p.Qq*f.Qq+p.qq*f.qq)

   w.star<-c(tau.QQ,tau.Qq,tau.qq)

   fit.lm<-lm(y.star~Q.star,weights=w.star,na.action=na.omit) ## update this
#   fit.lm<-lm(y.star~Sex.star+Q.star,weights=w.star,na.action=na.omit) ## update this
#   fit.lm<-lm(y.star~Sex.star+Age.star+Q.star,weights=w.star,na.action=na.omit) ## update this
   b.lm<-coef(fit.lm)
   s.lm<-sqrt(sum(w.star*(resid(fit.lm))^2)/(sum(w.star)))

   logL<-sum(log(p.QQ*f.QQ+p.Qq*f.Qq+p.qq*f.qq))
   conv1<-sum((b.lm/b.lm.old-1)^2)
   conv2<-(logL-logL.old)
   if(ptrace==1) print(c(b.lm,s.lm,conv1,conv2,logL))
   iter<-iter+1
}
