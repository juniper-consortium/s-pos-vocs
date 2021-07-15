# mixture distribution combination

pmixture = function(q, cdfFn, ...) {
  params = rlang::list2(...)
  cdfs = purrr::pmap(params,function(...) return(function(x) cdfFn(x,...)))
  tmp = sapply(cdfs, function(z) z(q))
  if(is.array(tmp)) {
    tmp = rowMeans(tmp)
  } else {
    tmp = mean(tmp)
  }
  return(tmp)
}

# rough test
# pmixture(c(-1), pnorm, mean=c(0,2,4), sd=c(1,2,3))



qmixture = function(p, bounds, cdfFn, ... ) {
  quantGenFn = function(q) return(function(x) pmixture(x, cdfFn, ...)-q)
  estimates = sapply(p, function(x) tryCatch(uniroot(quantGenFn(x), interval=bounds)$root,error = function(e) NaN))
  return(tibble(
    quantile = p,
    estimate = estimates,
    method = "ensemble"
  ))
}

# 
# qmixture(c(0.025,0.5,0.975), c(-1000,1000), pnorm, mean=c(0,2,4), sd=c(1,2,3))

#plot(pmixture(seq(-5,5,length.out = 1001), pnorm, mean=c(0,2,4), sd=c(1,2,3)))
