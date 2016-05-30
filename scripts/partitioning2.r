make.block=function(a,b) outer(a,b,"+")
scalarprod=function(x,y) mean(x*y)

findCuttofs <- function(chr, cuttof, nGenes, nTimes)
{
  filename = paste("../partitions/", chr, "/partitions_sigma_0.Rda", sep = "")
  partitions <- readRDS(filename)
  
  r1r2 = partitions$alpha - partitions$beta
  
  percentage = length(which(r1r2 > cuttof)) / length(r1r2)
  
  sigma2 = cuttof / (2 * log(nGenes*nTimes))
  
  c(percentage, sigma2)
}

model1 = function(M)
{
  n = nrow(M)

  temporalMean = apply(M, 2, mean)
  
  make.block(rep(0,n), temporalMean)
}

residueModel1 = function(M, BIC=T, sigma2, nDataPoints)
{
  n = nrow(M)
  m = ncol(M)
  
  temporalMean = apply(M, 2, mean)
  
  residue = sum((M - model1(M))^2)
  
  # if(BIC) residue = residue + sigma2 * (m + 1)
  if(BIC) residue = residue + sigma2 * (m + 1) * log(nDataPoints)
  
  residue
}

model2 = function(M) {
  n = nrow(M)
  m = ncol(M)
  
  temporalMean = apply(M, 2, mean)
  
  sin = sin(2*pi/24*seq(0,46,by=2))
  cos = cos(2*pi/24*seq(0,46,by=2))
  
  a.s=apply(M - matrix(temporalMean,n,m, byrow=T), 2, scalarprod, y=sin)*2
  a.c=apply(M - matrix(temporalMean,n,m, byrow=T), 2, scalarprod, y=cos)*2
  
  alpha = mean(a.s)
  beta = mean(a.c)
  
  make.block(alpha*sin + beta*cos, temporalMean)
}

residueModel2 = function(M, BIC=T, sigma2, nDataPoints)
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
  
  residue = sum((M - model2(M))^2)

  # if(BIC) residue = residue + sigma2 * (m + 2 + 1)
  if(BIC) residue = residue + sigma2 * (m + 2 + 1)*log(nDataPoints)
  
  c(residue, alpha, beta)
}

models = function(j, M, min.size, max.size, BIC=T, s2, nDataPoints)
{
  m = ncol(M)
	M = as.matrix(M[, j:m])
	m = ncol(M)

	type = 1
	
	alpha = 0
	beta = 0
	
	res.1 = 0
	res.2 = 0
	
	if(m < min.size | m > max.size) residue = Inf
	else {
		# MODEL 1
		res.1 = residueModel1(M, BIC, s2, nDataPoints)

		# MODEL 2
		result2 = residueModel2(M, BIC, s2, nDataPoints)
		res.2 = result2[1]
		alpha = result2[2]
		beta = result2[3]
		
		r  = c(res.1, res.2)
		imin = which.min(r)
		residue = r[imin]
		type = imin
	}
	
	c(residue, type, res.1, res.2, alpha, beta)
}

partitioning = function(M, min.size, max.size, BIC, cuttof)
{
  print(paste('Computing the partitions of', chr, 'for cuttof:', cuttof))
  
  n = nrow(M)
  m = ncol(M)
  
  sigmaAndPercentage = findCuttofs(chr = chr, cuttof = cuttof, nGenes = m, nTimes = n)
  
  percentage = round(sigmaAndPercentage[1] * 100, 2)
  sigma2 = sigmaAndPercentage[2]
  sigma = sqrt(sigma2)
  
	print(paste('Number of genes :', m))
	print(paste('Maximum partition size :', max.size))
	print(paste('Percentage of blocks with circadian model: ', percentage, '%', sep = ''))
	print(paste('Value of sigma :', sigma))
	
	scores	= rep(0,m)
	change	= rep(0,m)
	type	= rep(0,m)
	
	res.1 = rep(0, m)
	res.2 = rep(0, m)
	
	alpha = rep(0, m)
	beta = rep(0, m)

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
	while(bk[1]>1)
	{
		tmp = bk[1]-1
		bk  = c(change[tmp],bk)
		bk.t = c(type[tmp],bk.t)
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
	     alphas = alpha,
	     betas = beta,
	     cuttof = cuttof,
	     percentage.prior = percentage,
	     percentage.resulted = rPercentage,
	     sigma = sigma,
	     max.size = max.size
	     )
}