make.block=function(a,b) outer(a,b,"+")
scalarprod=function(x,y) mean(x*y)

findCuttofs <- function(chr, percentage, nGenes, nTimes)
{
  filename = paste("../partitions/", chr, "/partitions_no_penalty.Rda", sep = "")
  partitions <- readRDS(filename)
  
  r1r2 = partitions$res1 - partitions$res2
  
  cuttof <- quantile(r1r2, 1 - percentage/100)
  
  sigma2 = cuttof / (2 * log(nGenes*nTimes))
  
  c(cuttof, sigma2)
}

model1 = function(M)
{
  n = nrow(M)
  
  temporalMean = apply(M, 2, mean)
  
  make.block(rep(0,n), temporalMean)
}

residueModel1 = function(M, BIC, sigma2, nDataPoints)
{
  m = ncol(M)
  
  residue = sum((M - model1(M))^2)
  
  if(BIC) score = residue + sigma2 * (m + 1) * log(nDataPoints)
  else score = residue
  
  c(score, residue)
}

model2 = function(M) {
  n = nrow(M)
  m = ncol(M)
  
  temporalMean = apply(M, 2, mean)
  
  sin = sin(2*pi/24*seq(0,46,by=2))
  cos = cos(2*pi/24*seq(0,46,by=2))
  
  a.s = apply(M - matrix(temporalMean,n,m, byrow=T), 2, scalarprod, y=sin)*2
  a.c = apply(M - matrix(temporalMean,n,m, byrow=T), 2, scalarprod, y=cos)*2
  
  alpha = mean(a.s)
  beta = mean(a.c)
  
  make.block(alpha*sin + beta*cos, temporalMean)
}

residueModel2 = function(M, BIC, sigma2, nDataPoints)
{
  n = nrow(M)
  m = ncol(M)
  
  temporalMean = apply(M, 2, mean)
  
  sin = sin(2*pi/24*seq(0,46,by=2))
  cos = cos(2*pi/24*seq(0,46,by=2))
  
  a.s = apply(M - matrix(temporalMean,n,m, byrow=T), 2, scalarprod, y=sin)*2
  a.c = apply(M - matrix(temporalMean,n,m, byrow=T), 2, scalarprod, y=cos)*2
  
  alpha = mean(a.s)
  beta = mean(a.c)
  
  residue = sum((M - make.block(alpha*sin + beta*cos, temporalMean))^2)
  
  if(BIC) score = residue + sigma2 * (m + 2 + 1)*log(nDataPoints)
  else score = residue
  
  c(score, residue, alpha, beta)
}

models = function(j, M, min.size, max.size, BIC, s2, nDataPoints)
{
  m = ncol(M)
  M = as.matrix(M[, j:m])
  m = ncol(M)
  
  type = 1
  
  alpha = 0
  beta = 0
  
  score.1 = Inf
  score.2 = Inf
  
  res.1 = 0
  res.2 = 0
  
  if(m < min.size | m > max.size) score = Inf
  else {
    # MODEL 1
    result1 = residueModel1(M, BIC, s2, nDataPoints)
    score.1 = result1[1]
    res.1 = result1[2]
    
    # MODEL 2
    result2 = residueModel2(M, BIC, s2, nDataPoints)
    score.2 = result2[1]
    res.2 = result2[2]
    alpha = result2[3]
    beta = result2[4]
    
    r = c(score.1, score.2)
    imin = which.min(r)
    score = r[imin]
    type = imin
  }
  
  c(score, type, res.1, res.2, alpha, beta)
}

partitioning = function(M, chr, BIC, sigma2, min.size=1, max.size=100, percentage=100)
{
  n = nrow(M)
  m = ncol(M)
  
  scores	= rep(0,m)
  change	= rep(0,m)
  type	= rep(0,m)
  
  res.1 = rep(0, m)
  res.2 = rep(0, m)
  
  alpha = rep(0, m)
  beta = rep(0, m)
  
  sigma = sqrt(sigma2)
  
  for (k in 1:m)
  {
    if(k == 1) s = c(0)
    if(k > 1)  s = c(0, scores[1:(k-1)])
    
    ss = sapply(1:k,
                models,
                M=as.matrix(M[,1:k]),
                min.size=min.size,
                max.size=max.size,
                BIC=BIC,
                s2=sigma2,
                nDataPoints=n*m)
    
    s.tmp = s + ss[1,]
    
    jmin = which.min(s.tmp)
    
    scores[k] = s.tmp[jmin]
    change[k] = jmin
    type[k] = ss[2,jmin]
    res.1[k] = ss[3, jmin]
    res.2[k] = ss[4, jmin]
    alpha[k] = ss[5, jmin]
    beta[k] = ss[6, jmin]
  }
  
  #backtrack
  bk = change[m]
  bk.t = type[m]
  p.alpha = alpha[m]
  p.beta = beta[m]
  while(bk[1]>1)
  {
    tmp = bk[1]-1
    bk  = c(change[tmp], bk)
    bk.t = c(type[tmp], bk.t)
    p.alpha = c(alpha[tmp], p.alpha)
    p.beta = c(beta[tmp], p.beta)
  }
  
  sizes = diff(c(bk, ncol(M)+1))
  scores = scores / (n*m)
  
  rPercentage = round(length(which(type == 2)) / length(type) * 100, 2)
  
  list(sizes = sizes,
       block.types= bk.t,
       types = type,
       scores = scores,
       res1 = res.1,
       res2 = res.2,
       alphas = p.alpha,
       betas = p.beta,
       cuttof = cuttof,
       percentage.prior = percentage,
       percentage.resulted = rPercentage,
       sigma = sigma,
       max.size = max.size,
       chromosome = chr,
       nGenes = m
  )
}