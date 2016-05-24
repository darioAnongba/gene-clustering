make.block=function(a,b) outer(a,b,"+")

scalarprod=function(x,y) mean(x*y)


model1 = function(M, BIC=T, sigma2)
{
  n = nrow(M)
  m = ncol(M)
  
  temporalMean = apply(M,2,mean)
  
  make.block(rep(0,n), temporalMean)
}

residueModel1 = function(M, BIC=T, sigma2)
{
  n = nrow(M)
  m = ncol(M)
  
  temporalMean = apply(M,2,mean)
  
  residue = sum((M - make.block(rep(0,n), temporalMean))^2)
  
  #if(BIC) residue = residue + sigma2 * m*log(n*m)
  if(BIC) residue = residue + sigma2 * (m + 1)
  
  residue
}

model2 = function (M, BIC=T, sigma2) {
  n = nrow(M)
  m = ncol(M)
  
  temporalMean = apply(M,2,mean)
  
  sin = sin(2*pi/24*seq(0,46,by=2))
  cos = cos(2*pi/24*seq(0,46,by=2))
  
  a.s=apply(M - matrix(temporalMean,n,m, byrow=T), 2, scalarprod, y=sin)*2
  a.c=apply(M - matrix(temporalMean,n,m, byrow=T), 2, scalarprod, y=cos)*2
  
  alpha = mean(a.s)
  beta = mean(a.c)
  
  make.block(alpha*sin + beta*cos, temporalMean)
}

residueModel2 = function(M, BIC=T, sigma2)
{
  n = nrow(M)
  m = ncol(M)
  
  temporalMean = apply(M, 2, mean)
  
  sin = sin(2*pi/24*seq(0,46,by=2))
  cos = cos(2*pi/24*seq(0,46,by=2))
  
  a.s = apply(M - matrix(temporalMean, n, m, byrow=T), 2, scalarprod, y=sin)*2
  a.c = apply(M - matrix(temporalMean, n, m, byrow=T), 2, scalarprod, y=cos)*2
  
  alpha = mean(a.s)
  beta = mean(a.c)
  
  residue = sum((M - make.block(alpha*sin + beta*cos, temporalMean))^2)
  
  #if(BIC) residue = residue + sigma2 * (m + 2)*log(n*m)
  if(BIC) residue = residue + sigma2 * (m + 2 + 1)
  
  c(residue, alpha, beta)
}

models = function(j, M, min.size, max.size, BIC=T, s2)
{
  m = ncol(M)
  M = as.matrix(M[,j:m])
  m = ncol(M)
  
  type = 1
  r1minusr2 = 0
  
  if(m < min.size | m > max.size) residues = Inf
  else {
    # MODEL 1
    res.1 = residueModel1(M, BIC, s2)
    
    # MODEL 2
    result2 = residueModel2(M, BIC, s2)
    res.2 = result2[1]

    r  = c(res.1, res.2)
    imin = which.min(r)
    residues = r[imin]
    type = imin
    r1minusr2 = res.1 - res.2
  }
  
  c(residues, type, r1minusr2)
}

partitioning = function(M, min.size=2, max.size=30, BIC=T, sigma2=0)
{
  if(is.null(sigma2) & BIC) sigma2=var(as.vector(M))
  n = nrow(M)
  m = ncol(M)
  
  scores	= rep(0,m)
  change	= rep(0,m)
  type	= rep(0,m)
  
  r1minusr2 = rep(0, m)

  for (k in 1:m)
  {
    if(k == 1) s = c(0)
    if(k > 1)  s = c(0, scores[1:(k-1)])
    ss = sapply(1:k, models, M=as.matrix(M[,1:k]), min.size=min.size, max.size=max.size, BIC=BIC, s2=sigma2)
    s.tmp = s + ss[1,]
    jmin = which.min(s.tmp)
    scores[k] = s.tmp[jmin]
    change[k] = jmin
    type[k] = ss[2,jmin]
    
    r1minusr2[k] = ss[3, jmin]
  }
  
  #backtrack
  bk=change[m]
  bk.t=type[m]
  while(bk[1]>1)
  {
    tmp = bk[1]-1
    bk  = c(change[tmp],bk)
    bk.t = c(type[tmp],bk.t)
  }
  
  sizes=diff(c(bk, ncol(M)+1))
  scores=scores/(n*m)
  
  list(sizes=sizes,
       block.type=bk.t,
       start.pos=bk,
       changes=change,
       scores=scores,
       types=type,
       r1r2=r1minusr2)
}
