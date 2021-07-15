#
.md5obj = function(obj) {
  as.character(openssl::md5(serialize(obj, connection = NULL)))
}

.cache <- new.env(parent=emptyenv())

#' Cache an expression to memory and to disk.
#'
#' @param expr - the expression that needs to be cached
#' @param inputs - a list(...) of dependencies that must be the same if any 
#' @param ... - no unnamed parameters allowed after inputs,
#' @param nocache - defeat the caching and force a refresh (defaults to option("nocache"))
#' @param cache - the directory for the caching (defaults to option("cache.dir"))
#' @param stale - the number fo days before a result is stale (defaults to option("cache.stale")) 
#'
#' @return the result of expr or a cached version
cached = function (
  expr,
  inputs = list(),
  ...,
  nocache = getOption("nocache", default=FALSE),
  cache = getOption("cache.dir", default=tempdir()),
  stale = getOption("cache.stale", default=Inf)
  )
{
  code = deparse(substitute(expr))
  
  inputs2 = rlang::list2(...)
  if (length(inputs2) > 0) stop("any dependencies must be explictily stated as an input in inputs=list(...)")
  
  md5code = .md5obj(code)
  
  if(!stringr::str_ends(cache,"/"))
    cache = paste0(cache,"/")

  dir.create(cache, recursive = TRUE, showWarnings = FALSE)

  md5params = .md5obj(inputs)

  cacheName = paste(md5code,md5params,sep = "-")
  
  # Check memory cache first
  if (nocache) suppressWarnings(rm(list=cacheName,envir = .cache,inherits = FALSE))
  
  if (exists(cacheName,where = .cache)) {
    obj = get(cacheName, envir=.cache)
    if (attr(obj,".cached-at") < Sys.Date()-stale+1) {
      rm(cacheName,where = .cache)  
    } else {
      message("already loaded: ",cacheName)
      return(obj)
    }
  } 
  
  path = paste0(cache,cacheName,".rda")

  if (nocache) unlink(path)
  if (file.exists(path)) {
    mtime = as.Date(file.info(path)$mtime)
    if (mtime < Sys.Date()-stale+1) unlink(path)
  }

  if (file.exists(path)) {
    message("using cached item: ",path)
    obj = readRDS(path)
    assign(cacheName, obj, envir=.cache)
  } else {
    message("caching item: ",path)
    obj = expr
    attr(obj,".cached-at") <- Sys.Date()
    saveRDS(obj, path)
    assign(cacheName, obj, envir=.cache)
  }
  
  return(obj)
}


