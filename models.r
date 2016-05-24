make.block=function(a,b) outer(a,b,"+")
 
scalarprod = function(x,y) mean(x*y)

models = function(x, BIC, s2)
{
  #For all models
  b = mean(x)
  subs = (x - b)
  m = length(x)
  
  # MODEL 1
  res.1 = sum(subs^2)
  if(BIC) res.1 = res.1 + s2 * log(m)

  # MODEL 2
  s = sin(2*pi/24*seq(0,46,by=2))
  c = cos(2*pi/24*seq(0,46,by=2))
  
  a.s = scalarprod(subs, s) * 2
  a.c = scalarprod(subs, c) * 2
  
  alpha = mean(a.s)
  beta = mean(a.c)
  
  res.2 = sum((x - make.block(alpha*s + beta*c, b))^2)
  if(BIC) res.2 = res.2 + s2 * 2*log(m)
  
  c(res.1, res.2)
}

fitting = function(M, BIC=F, s2=1)
{
  scores = apply(M, 2, models, BIC, s2)
  
  colnames(scores) = colnames(M)
  rownames(scores) = c('Model 1', 'Model 2')
  
  scores
}

applyModelsHeatmap = function(x, BIC, s2) {
  scores = models(x, BIC, s2)
  
  if(scores[1] < scores[2]) 0 else 1
}

modelsHeatmap = function(M, BIC=T, s2)
{
  apply(M, 2, applyModelsHeatmap, BIC, s2)
}