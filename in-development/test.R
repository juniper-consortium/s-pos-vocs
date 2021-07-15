# just check the caching ....
here::i_am("test.R")
source(here::here("utils-cache.R"))
options(cache.dir = here::here("working"))

tmp=cached({mtcars}, inputs=c(1, 2, 3))
tmp2=cached(ggplot2::diamonds, 1, 2, 3)


devtools::load_all("~/Git/uk-covid-datatools/")
ukcovidtools::reload("~/Git/uk-covid-datatools/config.yml")
dpc$spim$getLatestRawFile(dpc$spim$filter$voc351, here::here("input"))
