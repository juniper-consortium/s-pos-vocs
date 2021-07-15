# OLD CODE ----
# reinfection


```{r}

#' Interpret S gene status according to various cut off values
#' function to help interpret S gene CT values in context of N gene and ORF gene to give S gene status. 
#' With the defaults this produces the same result as the sgtf_30 column in the source SGTF line list
#' Defaults are S:30,ORF:30,N:30,Control:Inf
#'
#' @param sGeneLineList - a dataframe includeing 
#' @param S_CT - S gene detected when P2CH3CQ <= this value
#' @param ORF1ab_CT - ORF1ab gene detected when P2CH1CQ <= this value
#' @param N_CT - N gene detected when P2CH2CQ <= this value
#' @param Control_CT - control sample is positive when P2CH4CQ <= this value
#'
#' @return - the same dataframe with additional columns including "sGene" and "result"
#'
#' @examples coxData = coxData %>% interpretSGene()
interpretSGene = function(sGeneLineList, S_CT = 40, ORF1ab_CT = 40, N_CT = 40, Control_CT = Inf) {
  sGeneLineList %>% 
    mutate(
      ORF1ab_CT_threshold = ORF1ab_CT,
      N_CT_threshold = N_CT,
      S_CT_threshold = S_CT,
      S_pos = P2CH3CQ > 0 & P2CH3CQ <= S_CT,
      S_undetect = P2CH3CQ == 0,
      N_pos = P2CH2CQ > 0 & P2CH2CQ <= N_CT,
      ORF1ab_pos = P2CH1CQ > 0 & P2CH1CQ <= ORF1ab_CT,
      Control_pos = P2CH4CQ > 0 & P2CH4CQ <= Control_CT,
      sGene = case_when(
        is.na(P2CH1CQ) ~ "Unknown",
        S_pos & N_pos & ORF1ab_pos & Control_pos ~ "Positive",
        S_undetect & N_pos & ORF1ab_pos & Control_pos ~ "Negative",
        TRUE ~ "Equivocal"
      )
    ) %>% mutate(
      result = case_when(
        is.na(P2CH1CQ) ~ "Unknown",
        !Control_pos ~ "No control",
        TRUE ~ paste0(ifelse(S_pos,"S+","S-"),ifelse(N_pos,"N+","N-"),ifelse(ORF1ab_pos,"ORF+","ORF-")))
    )
}


# sgll - the s gene line list
# ll2 - the cases line list (which essentially has only the first occurrence of a case)

diff = sgll %>% left_join(ll2 %>% select("FINALID",first_specimen_date = "specimen_date"), by=c("FINALID"))
diff2 = diff %>% 
  # look for positive test results in the S gene list which are much later than original cases
  # these are possible reinfections
  filter(specimen_date > first_specimen_date + 56) %>% 
  group_by(FINALID) %>% 
  # get rid of any further repeat testing etc.
  filter(specimen_date == min(specimen_date)) %>%
  ungroup() %>%
  interpretSGene(40,40,40)




reinfect = sgll %>% 
  group_by(specimen_date) %>% 
  # baseline is number of S gene results on a given day.
  # Is this a good baseline?
  summarise(baseline = n()) %>%
  left_join(
    diff2 %>% 
      # hold onto sgtf status for sub
      group_by(specimen_date,sGene) %>% 
      summarise(cases = n()), 
    by="specimen_date")

# combine s gene status for initial plot.
reinfectAll = reinfect %>% 
  group_by(specimen_date,baseline) %>% 
  summarise(cases = sum(cases)) %>%
  filter(!is.na(cases)) %>%
  mutate(binom::binom.confint(cases, baseline, conf.level = 0.95, methods = "wilson"))

ggplot(reinfectAll,aes(x=specimen_date, y=mean*100))+geom_errorbar(aes(ymin=lower*100,ymax=upper*100),colour="grey50")+geom_point()+ylab("percent Pillar 2 tests that could be reinfections")+xlab("follow up specimen date")

# repeat keeping the s gene breakdown. Equivocal category is difficult to interpret here but more likey to be S pos.
reinfect2 = reinfect %>% 
  filter(!is.na(cases)) %>%
  mutate(binom::binom.confint(cases, baseline, conf.level = 0.95, methods = "wilson"))

ggplot(reinfect2,aes(x=specimen_date, y=mean*100))+geom_errorbar(aes(ymin=lower*100,ymax=upper*100),colour="grey50")+geom_point()+ylab("percent Pillar 2 tests that could be reinfections")+facet_wrap(vars(sGene))+xlab("follow up specimen date")

# Feel sick. Lie down.

delay = diff2 %>% mutate(delay = as.numeric(specimen_date-first_specimen_date)) %>% group_by(specimen_date) %>% summarise(
  lower = quantile(delay,0.25),
  upper = quantile(delay,0.75),
  median = quantile(delay,0.5)
)
ggplot(delay,aes(x=specimen_date, y=median, ymin=lower, ymax=upper))+geom_errorbar(colour="grey50")+geom_point()+ylab("median delay (+IQR)")+xlab("follow up specimen date")

delayAll = diff2 %>% mutate(delay = as.numeric(specimen_date-first_specimen_date)) %>% group_by(specimen_date,sGene) %>% summarise(
  lower = quantile(delay,0.25),
  upper = quantile(delay,0.75),
  median = quantile(delay,0.5)
)
ggplot(delayAll,aes(x=specimen_date, y=median, ymin=lower, ymax=upper))+geom_errorbar(colour="grey50")+geom_point()+ylab("median delay (+IQR)")+xlab("follow up specimen date")+facet_wrap(vars(sGene))

```

##

```{r}

reinfectLTLA = sgll %>% 
  inner_join(ll2 %>% select(FINALID,LTLA_code), by="FINALID") %>%
  # baseline is number of S gene results on a given day.
  # Is this a good baseline?
  group_by(LTLA_code,specimen_date) %>% 
  summarise(baseline = n()) %>%
  left_join(
    diff2 %>% 
      inner_join(ll2 %>% select(FINALID,LTLA_code), by="FINALID") %>%
      # hold onto sgtf status for sub
      group_by(LTLA_code,sGene,specimen_date) %>% 
      summarise(cases = n()), 
    by=c("LTLA_code","specimen_date"))

reinfectLTLA56 = reinfectLTLA %>% filter(!is.na(sGene) & specimen_date > Sys.Date()-56) %>% ungroup() %>% tidyr::complete(LTLA_code,specimen_date,sGene,fill=list(cases=0,baseline=NA)) %>% mutate(prop = ifelse(cases==0,0,cases/baseline))
#reinfect56 = reinfectAll %>% filter(specimen_date > Sys.Date()-56) %>% mutate(cases = ifelse(is.na(cases),0,cases), prop = cases/baseline)
reinfect56 = reinfect %>% filter(specimen_date > Sys.Date()-56) %>% ungroup() %>% tidyr::complete(specimen_date,sGene,fill=list(cases=0,baseline=NA)) %>% mutate(prop = ifelse(cases==0,0,cases/baseline))

reinfectLTLAP = reinfectLTLA56  %>% group_by(LTLA_code,sGene) %>% group_modify(function(d,g,...) {
  
  if(length(d$baseline) > 1) {
    tryCatch({
      tmp = t.test(
        x = d$prop,
        y = reinfect56 %>% filter(sGene == g$sGene) %>% pull(prop),
        paired=FALSE)
    },error=browser)
    return(tidy(tmp))
  } else {return(tibble())}
  
})


map = dpc$geog$getMap("LAD19")
# combine s gene status for initial plot.
# reinfectLTLA = reinfectLTLA %>% 
#   mutate(cases = ifelse(is.na(cases),0,cases)) %>%
#   mutate(binom::binom.confint(cases, baseline, conf.level = 0.95, methods = "wilson"))

map2 = map %>% inner_join(reinfectLTLAP, by=c("code"="LTLA_code"))

#ggplot(map2,aes(fill=mean*100))+geom_sf()+scale_fill_gradient(limits=c(0,10), oob=scales::squish)
ggplot()+
  geom_sf(data = map %>% filter(code %>% stringr::str_starts("E")),fill="black")+
  geom_sf(data = map2,mapping = aes(fill=estimate*100,alpha=1-p.value))+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint = 0,limits=c(-2,2),oob=scales::squish)+
  scale_alpha_identity()+facet_wrap(vars(sGene))

```


```{r}
# plotWithLabels = function(ltlaCounts, xlim = c(limits[2],limits[4]), ylim=c(limits[1],limits[3]), limits = c(NA,NA,NA,NA), labelsPerPlot = 4, nrow=2) {
#   
#   ltlaCounts = ltlaCounts %>% filter(category %in% c("unsequenced S+","B.1.351","B.1.617.2"))
#   ladMap = dpc$geog$getMap("LAD19")
#   ltlas = ladMap %>% as_tibble() %>% select(code,name) %>% distinct() %>%
#     mutate(label = code %>% stringr::str_replace("E0(.)0+(.+)","\\1.\\2"))
#   categories = ltlaCounts %>% select(category) %>% distinct() %>% arrange(category)
# 
#   ltlaCounts2 = ltlas %>% select(-name) %>% inner_join(categories, by=character()) %>%
#     left_join(ltlaCounts %>% select(-name), by=c("code","category")) %>% 
#     mutate(count = ifelse(is.na(count),0,count))
#   
#   mapVoc = ltlaCounts2 %>% left_join(ladMap, by=c("code"), suffix=c("",".map")) %>% 
#     filter(code %>% stringr::str_starts("E")) %>% sf::st_as_sf()
#   
#   labels = mapVoc %>% 
#     sf::st_centroid() %>%
#     mutate(x=sf::st_coordinates(.)[,"X"], y=sf::st_coordinates(.)[,"Y"]) %>% 
#     filter(x > ifna(xlim[1],-Inf) & x < ifna(xlim[2],Inf) ) %>%
#     filter(y > ifna(ylim[1],-Inf) & y < ifna(ylim[2],Inf) ) %>% 
#     as_tibble() %>%
#     group_by(category) %>% 
#     arrange(desc(count)) %>% 
#     filter(row_number() <= labelsPerPlot) %>%
#     ungroup()
#   
#   p1 = lapply(categories$category, function(cat) {
#     breaks = mapVoc %>% filter(category==cat & count>0) %>% pull(count) %>% quantile(c(0.5,0.75,0.9,0.95,0.975))
#     breaks = sort(unique(round(c(1,unname(breaks[breaks>1])))))
#     plotlabels = c("0",ifelse(c(0,breaks)+1 == c(breaks,Inf), as.character(breaks) , paste0( as.character(c(0,breaks)+1),c(rep("-",length(breaks)),"+"),c(as.character(breaks),""))))
#     tmp = mapVoc %>% filter(category==cat) %>% mutate(count2=cut(count,breaks = c(-Inf,0,breaks,Inf),labels = plotlabels,include.lowest = FALSE))
#     #browser()
#     ggplot()+
#     #geom_sf(data=dpc$geog$getMap("LAD19") ,size=0.05)+
#       geom_sf(data=tmp,mapping=aes(fill=count2),size=0.1,colour="white")+
#       scale_fill_viridis_d(option = "inferno", name="")+
#       standardPrintOutput::defaultMapLayout()+
#       standardPrintOutput::narrower()+
#       standardPrintOutput::smallLegend(textSize = 6)+ #spaceLegend = 0.8)+
#       facet_wrap(vars(category))+
#       coord_sf(xlim = xlim,ylim= ylim)+
#       ggrepel::geom_label_repel(data = labels %>% filter(category==cat), aes(x=x,y=y,label=label),min.segment.length = 0,segment.colour = "cyan",fill="#80808080",colour="cyan",size=1.5)
#   })
#   
#   
#   #print(labs %>% mutate(l = paste0(label,": ",name)) %>% pull(l) %>% paste0(collapse = "; "))
#   p4 = standardPrintOutput::simpleFigureTable(labels %>% as_tibble() %>% select(label,name) %>% distinct() %>% arrange(label),pts = 5,unwrapped = TRUE)
#   #browser()
#   return(patchwork::wrap_plots(c(p1,list(p4)),nrow=nrow)+patchwork::plot_annotation(title = paste0("cases between ",earliestDate," and ",lastDate))) #, widths = c(3,3,3,2)))#+patchwork::plot_layout(ncol = 2))
# }
```





## TODO----
## OLD CODE ----

# sgllEra = loader$getSGeneEras()


# take the latest values of the S pos for each person
# this could generate dups
# sgPosUniq = sgll %>% filter(sgtf_under30CT == "0") %>% 
#   group_by(FINALID) %>% 
#   filter(specimen_date==min(specimen_date,na.rm=TRUE))

# sgPosUniq = sgll %>% filter(sgtf_under30CT == "0") %>% 
#   group_by(FINALID) %>% 
#   arrange(desc(specimen_date)) %>%
#   filter(row_number()==1) %>% select(FINALID, specimen_date)


# isolationDelay = 
#   bind_rows(
#     ctas %>% filter(category=="case") %>%
#       select(finalid = genomic_finalid, date = genomic_specimen_date, isolation_date = start_of_isolation_date),
#     ctas %>% filter(category=="case") %>%
#       select(finalid = sgtf_finalid, date = sgtf_specimen_date, isolation_date = start_of_isolation_date)
#   ) %>% 
#   filter(!is.na(finalid) & !is.na(date) & !is.na(isolation_date)) %>%
#   distinct() %>% mutate(isolation_delay = as.integer(isolation_date - date)) %>%
#   filter(isolation_delay > -14 & isolation_delay < 21) %>%
#   mutate(
#     isolation_cat = case_when(
#       isolation_delay <= -7 ~ "very low risk",
#       isolation_delay <= -2 ~ "low risk",
#       isolation_delay <= 0 ~ "asymptomatic spread",
#       isolation_delay <= 4 ~ "high risk",
#       TRUE ~ "ineffective",
#     ) %>% ordered(c("very low risk","low risk","asymptomatic spread","high risk","ineffective"))
#   )

# TODO: number of traced

# potentialReinfectn = sgll %>%
#   inner_join(ll2 %>% select(FINALID,first_specimen_date = specimen_date), by="FINALID") %>%
#   mutate(reinfect56 = specimen_date > first_specimen_date+56) %>%
#   select(finalid=FINALID, date=specimen_date, reinfect56)

# tmp2 = tmp %>% inner_join(ll2 %>% select(-asymptomatic_indicator),by=c("finalid"="FINALID"),suffix=c("",".ll")) %>% 
#   left_join(sgPosUniq, by=c("finalid"="FINALID"),suffix=c("",".sgene")) %>% # creates duplicates
#   inner_join(lsoaXY, by = c("LSOA_code"="code"), suffix = c("",".cent")) %>%
#   left_join(traced, by="finalid") %>%
#   mutate(ctas_status = ifelse(is.na(ctas_status),"unknown",ctas_status)) %>%
#   left_join(asymptomatic, by=c("finalid","date")) %>%
#   left_join(potentialReinfectn, by=c("finalid","date")) %>%
#   left_join(isolationDelay, by=c("finalid","date"))

# mapSGTFtoSgene = function(df) {
#   if ("pillar" %in% names(df)) {
#     df %>% mutate(
#       sGene = case_when(
#         pillar != "Pillar 2" ~ "unknown",
#         is.na(sgtf_under30CT) ~ "equivocal",
#         sgtf_under30CT == 1 ~ "negative",
#         sgtf_under30CT == 0 ~ "positive"
#       ) %>% ordered(c("positive","negative","equivocal","unknown"))
#     )
#   } else {
#     df %>% mutate(
#       sGene = case_when(
#         is.na(sgtf_under30CT) ~ "equivocal",
#         sgtf_under30CT == 1 ~ "negative",
#         sgtf_under30CT == 0 ~ "positive"
#       ) %>% ordered(c("positive","negative","equivocal"))
#     )
#   }
# }




# Bolton cluster

```{r}

# sgeneIncidence = sgll %>% select(FINALID,specimen_date,sgtf_under30CT) %>% left_join(ll2 %>% select(-specimen_date), by=c("FINALID"))
# sgeneIncidence2 = sgeneIncidence %>% filter(LTLA_code %in% areasOfConcern$`Bolton etc`)
# #sgeneIncidence2 %>% group_by(FINALID) %>% filter(n() > 1) %>% arrange(specimen_date) %>% filter(lag(specimen_date)<specimen_date-56) %>% View()
# sgeneIncidence3 = sgeneIncidence2 %>% group_by(FINALID) %>% filter(specimen_date == min(specimen_date)) %>% ungroup()
sgSymptIncid = sgllLoc %>% filter(LTLA_code %in% areasOfConcern$`NW Cluster`) %>% dpc$spim$getLineListIncidence(ageBreaks = c(6,18,45,55,75),subgroup = sgtf_under30CT,filterExpr = asymptomatic_indicator=="N",codeTypes = "LAD") 
#sgSymptIncid = sgllLoc %>% filter(LTLA_code %in% areasOfConcern$`Bolton etc`) %>% dpc$spim$getLineListIncidence(ageBreaks = c(6,18,45,55,75),subgroup = sgtf_under30CT,codeTypes = "LAD") 
tmp = sgSymptIncid %>% dpc$demog$findDemographics()
tmp2 = tmp %>% mutate(code="Bolton etc", name="Bolton etc") %>% ukcovidtools::covidStandardDateGrouping() %>% summarise(value=sum(value), population = sum(population))
tmp2 = tmp2 %>% mutate(subgroup=case_when(
  subgroup==0 ~ "S+",
  subgroup==1 ~ "S-",
  TRUE ~ "equivocal"
))
(tsp$plotIncidenceQuantiles(tmp2 %>% filter(date>"2021-03-01" & ageCat !="unknown"),colour = subgroup,denominatorExpr = population/1000000)+facet_wrap(vars(ageCat))+ylab("daily incidence per 1M")) %>%
  standardPrintOutput::saveThirdPageFigure(here::here("output/BoltonIncidencePer1M"))
tmp3 = tmp2 %>% tsp$aggregateAge() %>% tsp$estimateRt(window = 14,quick = FALSE)
tsp$plotRt(tmp3,colour = subgroup,dates =  "2021-03-01", window=14, rtlim = c(0,3)) %>%
  standardPrintOutput::saveThirdPageFigure(here::here("output/BoltonRt"))

ggplot(
  sgllLoc %>% #filter(LTLA_code %in% areasOfConcern$`NW Cluster`) %>% 
    filter(specimen_date>"2021-04-01") %>% mutate(sgene=case_when(
      sgtf_under30CT==0 ~ "S+",
      sgtf_under30CT==1 ~ "S-",
      TRUE ~ "equivocal"
    )) ,aes(x=age, colour=sgene))+geom_histogram(binwidth = 1)+facet_wrap(vars(sgene), scales = "free")

```

# HERE NEEDS MOVING TO DATA

```{r}

ltlaSgeneCountsAll = sgllLoc %>% dpc$spim$getLineListIncidence(ageBreaks = NULL,codeTypes = "LAD", subgroup = sGene)

ltlaToAoc = function(l) l %>% 
  dpc$demog$findDemographics() %>% 
  inner_join(aocDf, by = c("code"="ltlaCode")) %>% 
  select(-code,-name,-name.original) %>%
  rename(name = area) %>%
  mutate(
    name = as.character(name),
    code = name,
    codeType = "AOC"
  ) %>%
  ukcovidtools::covidStandardDateGrouping() %>%
  summarise(value = sum(value), population = sum(population)) %>%
  ungroup()

aocSgeneCounts = bind_rows(
  ltlaToAoc(ltlaSgeneCounts),
  symptomaticPillar2CtrySgeneCounts %>% dpc$demog$findDemographics()
)

aocSgeneCounts %>% filter(subgroup!="equivocal") %>% readr::write_csv(here::here("output/symptomatic-pillar2-cases-by-aoc-and-sgene.csv"))

aocSgeneCountsAll = bind_rows(
  ltlaToAoc(ltlaSgeneCountsAll),
  allPillar2CtrySgeneCounts %>% dpc$demog$findDemographics()
)

aocSgeneCountsAll %>% filter(subgroup!="equivocal") %>% readr::write_csv(here::here("output/all-pillar2-cases-by-aoc-and-sgene.csv"))

aocDf %>% readr::write_csv(here::here("output/ltla-to-aoc.csv"))



tsp$plotIncidenceQuantiles(aocSgeneCounts, colour=subgroup, dates="2021-03-01")+facet_wrap(vars(name))+scale_y_continuous(trans="log1p", breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000))
tsp$plotIncidenceQuantiles(aocSgeneCountsAll, colour=subgroup, dates="2021-03-01")+facet_wrap(vars(name))+scale_y_continuous(trans="log1p", breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000))







# aocSgeneCountsLogGrowth = aocSgeneCounts %>% tsp$logIncidenceStats(smoothingWindow = 21,growthRateWindow = 7, leftSided=TRUE)
# (tsp$plotGrowthRate(aocSgeneCountsLogGrowth %>% 
#               filter(name %in% c("NW Cluster","Bedford etc","Sefton & Lpl")) %>% 
#               filter(subgroup!="equivocal"),colour = subgroup,dates="2021-03-01",rlim=c(-0.25,0.25))+facet_wrap(vars(name))+standardPrintOutput::watermark(lab="provisional")) %>%   
#   standardPrintOutput::saveThirdPageFigure(here::here("output/sGeneAOCHighGrowthLog"))
# 
# 
# aocSgeneCountsRtAll = aocSgeneCountsAll %>% tsp$estimateRt(window = 14)
# (tsp$plotRt(aocSgeneCountsRt %>% filter(subgroup!="equivocal"),colour = subgroup,dates="2021-03-01",rtlim=c(0,3.5))+facet_wrap(vars(name))+standardPrintOutput::watermark(lab="provisional")) %>%   
#   standardPrintOutput::saveHalfPageFigure(here::here("output/sGeneAllAOCRt"))


#ggplot(ltlaSgeneCounts, aes(x=date,y=value))
# ltlasSgeneCounts is asymptomatic == "N"
#ltlaSgeneCounts2 = ltlaSgeneCounts %>% filter(subgroup != "equivocal") %>% tsp$logIncidenceStats(growthRateWindow = 7,nocache=TRUE)
ltlaSgeneCounts2 = ltlaSgeneCounts %>% filter(subgroup != "equivocal") %>% tsp$estimateGrowthRate(window = 14, growthRateWindow = 7)
# ltlaSgeneCounts2 = ltlaSgeneCounts2 %>% rename(
#   Growth.windowed.value = Growth.windowed.poisson,
#   Growth.windowed.SE.value = Growth.windowed.SE.poisson
#   ) #etc.

#ltlaSgeneCounts2 = ltlaSgeneCounts %>% filter(subgroup != "equivocal") %>% tsp$estimateGrowthRate(window = 14, growthRateWindow = 7)
```

```{r}
tmp = ltlaSgeneCounts2 %>% select(name,code,date,subgroup,starts_with("Growth"))

tmp2 = tmp %>% 
  filter(subgroup=="S+") %>% inner_join(
    tmp %>% filter(subgroup=="S-"),
    by = c("code","name","date"),
    suffix = c(".pos",".neg")
  ) %>% 
  mutate(
    alpha = 1-sqrt(as.numeric(maxDate-date)/as.numeric(maxDate-minDate))
  ) %>%
  mutate(
    # assume independence between S- and S+ growth rates and that CIs are normally distributes
    Growth.windowed.delta = Growth.windowed.value.pos - Growth.windowed.value.neg,
    Growth.windowed.SE.delta = sqrt(Growth.windowed.SE.value.pos^2+Growth.windowed.SE.value.neg^2),
    Growth.windowed.ProbPos.delta = 1-pnorm(0,mean=Growth.windowed.delta,sd=Growth.windowed.SE.delta),
    Growth.delta = Growth.value.pos - Growth.value.neg,
    Growth.SE.delta = sqrt(Growth.SE.value.pos^2+Growth.SE.value.neg^2),
    Growth.ProbPos.delta = 1-pnorm(0,mean=Growth.delta,sd=Growth.SE.delta)
  ) %>%
  rename(
    `S- Growth Rate` = Growth.value.neg,
    `S+ Growth Rate` = Growth.value.pos
  )

# thanks stackoverflow: https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
differenceStat <- function(m1,m2,s1,s2,n1,n2)
{
  se <- sqrt( (s1^2/n1) + (s2^2/n2) )
  # welch-satterthwaite df
  df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  t <- (m1-m2)/se 
  # difference of means * (1-p-value)
  out = (m1-m2) * (1-2*pt(-abs(t),df))    
  return(out) 
}

diagonals= tibble(
  intercept = c(0.2,0.1,0,-0.1,-0.2),
  slope = 1,
  label = forcats::as_factor(c("-20%","-10%","same","+10%","+20%"))
)


maxDate = max(tmp2$date)-4
minDate = maxDate-28

topNgrowth = tmp2 %>% filter(date==maxDate) %>% mutate(
  diffStat = differenceStat(
    Growth.windowed.value.pos,
    Growth.windowed.value.neg,
    Growth.windowed.SE.value.pos,
    Growth.windowed.SE.value.neg,
    Growth.windowed.Window.value.pos,
    Growth.windowed.Window.value.neg
  )) %>%
  arrange(desc(diffStat)) %>%
  filter(row_number() <= 10) #%>%
#select(name,code)

# topNgrowth = tmp2 %>% filter(date==maxDate) %>% mutate(
#   #diffStat = Growth.windowed.poisson.pos-Growth.windowed.poisson.neg
#     diffStat = differenceStat(
#       Growth.windowed.poisson.pos,
#       Growth.windowed.poisson.neg,
#       Growth.windowed.SE.poisson.pos,
#       Growth.windowed.SE.poisson.neg,
#       14,14
#     )
#   ) %>% 
#   ungroup() %>%
#   arrange(desc(diffStat)) %>% 
#   filter(row_number() <= 10) %>% 
#   select(name,code)

pathData = tmp2 %>% semi_join(topNgrowth, by="code") %>% 
  filter(date>minDate & date<=maxDate) %>% 
  
  #pointData = tmp2 %>% filter(date==maxDate & abs(Growth.windowed.value.pos) >0.01) %>% mutate(alpha = 1)
  labelData = tmp2 %>% 
  semi_join(topNgrowth, by="code") %>% 
  filter(date==maxDate) %>% 
  mutate(alpha = 1) %>%
  rename(
    `S- Growth Rate` = Growth.value.neg,
    `S+ Growth Rate` = Growth.value.pos
  )

ggplot(pathData, aes(x=`S- Growth Rate`, y=`S+ Growth Rate`, group=name, alpha=alpha))+
  coord_fixed(xlim=c(-0.2,0.2),ylim=c(-0.2,0.2))+
  geom_hline(yintercept = 0,colour="grey40")+
  geom_vline(xintercept = 0,colour="grey40")+
  geom_abline(data = diagonals, aes(intercept = intercept,slope=slope, colour=label), inherit.aes=FALSE)+
  #geom_abline(slope = 1, intercept = 0.1, colour="#FFB0B0")+
  #geom_abline(slope = 1, intercept = -0.1, colour="#B0FFB0")+
  #geom_abline(slope = 1, intercept = 0.2, colour="#FF8080")+
  #geom_abline(slope = 1, intercept = -0.2, colour="#B0FFB0")+
  #geom_point(data=pointData,colour="#B0B0FF")+
  ggrepel::geom_text_repel(data=labelData,aes(label=name),colour="grey40",min.segment.length = 0,box.padding = 2, 
                           max.overlaps = Inf)+
  geom_path(colour="blue")+
  geom_point(data=labelData,colour="blue")+
  
  guides(alpha="none")+
  xlab("S- growth rate")+
  ylab("S+ growth rate")

# ggplot(pathData, aes(x=Growth.windowed.value.neg, y=Growth.windowed.value.pos, group=name, alpha=alpha))+
#   coord_fixed(xlim=c(-0.15,0.15),ylim=c(-0.15,0.15))+
#   geom_hline(yintercept = 0,colour="grey40")+
#   geom_vline(xintercept = 0,colour="grey40")+
#   geom_abline(data = diagonals, aes(intercept = intercept,slope=slope, colour=label))+
#   #geom_abline(slope = 1, intercept = 0.1, colour="#FFB0B0")+
#   #geom_abline(slope = 1, intercept = -0.1, colour="#B0FFB0")+
#   #geom_abline(slope = 1, intercept = 0.2, colour="#FF8080")+
#   #geom_abline(slope = 1, intercept = -0.2, colour="#B0FFB0")+
#   #geom_point(data=pointData,colour="#B0B0FF")+
#   ggrepel::geom_text_repel(data=labelData,aes(label=name),colour="grey40",min.segment.length = 0,box.padding = 2, 
#                            max.overlaps = Inf)+
#   geom_path(colour="blue")+
#   geom_point(data=labelData,colour="blue")+
#   
#   guides(alpha="none")+
#   xlab("S- growth rate")+
#   ylab("S+ growth rate")

topNgrowth = tmp2 %>% filter(date==maxDate) %>%
  arrange(desc(Growth.windowed.ProbPos.delta*Growth.windowed.delta)) %>%
  filter(row_number() <= 10) %>%
  select(name,code)

(ggplot(tmp2 %>% filter(date>"2020-10-01"),
        aes(
          x=date,
          y=Growth.windowed.delta,
          group=name#,
          #colour=Growth.windowed.ProbPos.delta
        )
)+geom_line(alpha=0.2)+geom_hline(yintercept = 0,colour="red"))  %>%
  standardPrintOutput::saveHalfPageFigure(paste0(here::here("output/AllGrowthDelta")))

(ggplot(tmp2 %>% semi_join(topNgrowth) %>% filter(date>"2020-02-01"),
        aes(x=date,
            y=Growth.windowed.delta,
            ymin=Growth.windowed.delta-1.96*Growth.windowed.SE.delta,
            ymax=Growth.windowed.delta+1.96*Growth.windowed.SE.delta,
            group=name)
)+
    geom_ribbon(alpha=0.2,fill="grey50")+
    geom_line(mapping=aes(colour=name)))  %>%
  standardPrintOutput::saveHalfPageFigure(paste0(here::here("output/Top10GrowthDelta")))

rollingPercent = allVoc %>% #combinedSpositives %>% 
  filter(!is.na(type)) %>%
  group_by(LTLA_name, type, date=specimen_date) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  tidyr::complete(LTLA_name,type,date,fill=list(count=0)) %>%
  group_by(LTLA_name, type) %>% 
  arrange(date) %>%
  mutate(Roll.count = stats::filter(count,filter=rep(1,28)/28,sides=1)) %>%
  filter(!is.na(Roll.count)) %>% 
  group_by(LTLA_name,date) %>% 
  mutate(binom::binom.confint(Roll.count, sum(Roll.count), conf.level = 0.95, methods = "wilson")) %>% 
  mutate(weight = sum(Roll.count)) %>% 
  ungroup()

growthVPercent = tmp2 %>% inner_join(rollingPercent, by = c("name"="LTLA_name","date"="date")) %>%
  mutate(growthSigHigher = ifelse(Growth.windowed.ProbPos.delta>0.0975,"sig","ns"))

(ggplot(growthVPercent %>% filter(date==max(date)), aes(x=mean,y=Growth.windowed.delta,size=weight))+geom_point(alpha=0.3)+facet_wrap(vars(type),scales="free")+
    #geom_smooth(method = "lm",fullrange=TRUE)+
    coord_cartesian(xlim = c(0,1))) %>%
  standardPrintOutput::saveHalfPageFigure(paste0(here::here("output/SGeneGrowthDeltaVsPercentSequencing")))

(ggplot(growthVPercent %>% filter(date==max(date)), aes(x=Roll.count,y=Growth.windowed.delta,size=weight))+geom_point(alpha=0.3)+facet_wrap(vars(type))
  #geom_smooth(method = "lm",fullrange=TRUE)+
) %>%
  standardPrintOutput::saveHalfPageFigure(paste0(here::here("output/SGeneGrowthDeltaVsCountSequencing")))

ggplot(growthVPercent %>% filter(type %in% c("B.1.351","B.1.617.1","B.1.617.2","P.1 & P.2") & date==maxDate), aes(x=type,colour=growthSigHigher,y=mean))+geom_boxplot()+scale_y_log10()

```


```{r}
# map = dpc$geog$getMap("LAD19")
# 
# map2 = map %>% inner_join(rollingPercent, by=c("name"="LTLA_name"))
# ggplot(map2 %>% filter(date==maxDate))+geom_sf(aes(fill=mean))+facet_wrap(vars(type))+scale_fill_viridis_c()
#                    
```

```{r}

map3 = map %>% inner_join(tmp2, by="code")
(ggplot(map3 %>% filter(date==maxDate))+geom_sf(aes(fill=Growth.windowed.delta, colour=ifelse(Growth.windowed.ProbPos.delta>0.0975,"sig","ns")),size=0.3)+scale_color_manual(values = c("sig"="white","ns"="black"))+scale_fill_gradient2(high = "yellow", mid="black",low="cyan",limits=c(-0.15,0.15),oob=scales::squish)) %>%
  standardPrintOutput::saveFullPageFigure(paste0(here::here("output/SgeneGrowthDelta-England")))

lapply(names(locations), function(name) {
  limits = locations[[name]]
  xlim = c(limits[2],limits[4])
  ylim = c(limits[1],limits[3])
  (ggplot(map3 %>% filter(date==max(date)))+geom_sf(aes(fill=Growth.windowed.delta))+scale_fill_gradient2(high = "orange",low="blue",limits=c(-0.15,0.15),oob=scales::squish)+coord_sf(xlim=xlim,ylim=ylim)) %>%
    standardPrintOutput::saveHalfPageFigure(paste0(here::here("output/SgeneGrowthDelta-"),name))
})

```

```{r}
topNgrowth %>% inner_join(rollingPercent %>% filter(date==maxDate) %>% mutate(percent =sprintf("%1.1f%%",mean*100)) %>% select(LTLA_name,type,percent)  %>% pivot_wider(values_from = percent, names_from=type), by = c("name"="LTLA_name"))
```


```{r}
ltlaSgeneCounts3 = ltlaSgeneCounts2 %>% dpc$demog$findDemographics()
ltlaSgeneCounts3 = ltlaSgeneCounts3 %>% inner_join(rollingPercent %>% rename(variant=type), by = c("name"="LTLA_name","date"="date"))


tsp$plotGrowthIncidence(ltlaSgeneCounts3 %>% filter(variant == "B.1.617.2") %>% group_by(name, variant,subgroup),plotDates = maxDate,timespan = 28,colour = mean*weight,populationExpr = population,maxSize = 0.6)+facet_wrap(vars(variant,subgroup))+scale_colour_viridis_b()

```


```{r}

sPosll = ll2 %>% select(-specimen_date) %>% inner_join(sgPosUniq %>% filter(sgtf_under30CT == 0) %>% select(FINALID,specimen_date), by="FINALID")

sPosIncid = sPosll %>% dpc$spim$getLineListIncidence(ageBreaks = c(18,35,65,75)) %>% tsp$aggregateGender() %>% dpc$demog$findDemographics()
tsp$plotIncidenceQuantiles(covidTimeseries = sPosIncid %>% filter(name=="England"), dates = "2021-02-01", denominatorExpr = population/1000000, colour=ageCat)+scale_y_continuous(trans="log1p")

tsp$plotRt(covidRtTimeseries = sPosIncid %>% filter(name=="England"), dates = "2021-02-01", colour=ageCat,rtlim = c(0.25,1.5))
```

```{r}
notSNegll = ll %>% semi_join(sgll %>% filter(is.na(sgtf_under30CT) | sgtf_under30CT == 1), by="FINALID")
sNotNegIncid = notSNegll %>% dpc$spim$getLineListIncidence() %>% tsp$aggregateGender() %>% dpc$demog$findDemographics()
tsp$plotIncidenceQuantiles(covidTimeseries = sNotNegIncid, dates = "2021-02-01", denominatorExpr = population/1000000)+facet_wrap(vars(name))+scale_y_continuous(trans="log1p")

tsp$plotRt(covidRtTimeseries = sNotNegIncid, dates = "2021-02-01",rtlim = c(0.25,1.5))+facet_wrap(vars(name))

```

```{r}
sgUniq = sgll %>% group_by(FINALID) %>% 
  filter(CDR_Specimen_Request_SK==min(CDR_Specimen_Request_SK,na.rm=TRUE)) 
llWSgene = ll2 %>% filter(asymptomatic_indicator != "Y") %>% inner_join(sgll %>% filter(!is.na(sgtf_under30CT)) %>% mutate(sGene = ifelse(sgtf_under30CT == 1,"S neg","S pos")) %>% select(FINALID, sGene), by="FINALID")
sGeneIncid = llWSgene %>% dpc$spim$getLineListIncidence(subgroup = sGene) %>% dpc$demog$findDemographics()
tsp$plotIncidenceQuantiles(covidTimeseries = sGeneIncid, dates = "2021-02-01", colour = subgroup, denominatorExpr = population/1000000)+facet_wrap(vars(name))+scale_y_continuous(trans="log1p")

sRT14 = sGeneIncid %>% filter(name != "Unknown (England)") %>% tsp$estimateRt(quick = TRUE,window = 14)
tsp$plotRt(covidRtTimeseries = sRT14, dates = "2021-02-01",rtlim = c(0.25,1.5), colour = subgroup)+facet_wrap(vars(name))

```


# Bolton cluster

```{r}

# sgeneIncidence = sgll %>% select(FINALID,specimen_date,sgtf_under30CT) %>% left_join(ll2 %>% select(-specimen_date), by=c("FINALID"))
# sgeneIncidence2 = sgeneIncidence %>% filter(LTLA_code %in% areasOfConcern$`Bolton etc`)
# #sgeneIncidence2 %>% group_by(FINALID) %>% filter(n() > 1) %>% arrange(specimen_date) %>% filter(lag(specimen_date)<specimen_date-56) %>% View()
# sgeneIncidence3 = sgeneIncidence2 %>% group_by(FINALID) %>% filter(specimen_date == min(specimen_date)) %>% ungroup()
sgSymptIncid = sgllLoc %>% filter(LTLA_code %in% areasOfConcern$`NW Cluster`) %>% dpc$spim$getLineListIncidence(ageBreaks = c(6,18,45,55,75),subgroup = sgtf_under30CT,filterExpr = asymptomatic_indicator=="N",codeTypes = "LAD") 
#sgSymptIncid = sgllLoc %>% filter(LTLA_code %in% areasOfConcern$`Bolton etc`) %>% dpc$spim$getLineListIncidence(ageBreaks = c(6,18,45,55,75),subgroup = sgtf_under30CT,codeTypes = "LAD") 
tmp = sgSymptIncid %>% dpc$demog$findDemographics()
tmp2 = tmp %>% mutate(code="Bolton etc", name="Bolton etc") %>% ukcovidtools::covidStandardDateGrouping() %>% summarise(value=sum(value), population = sum(population))
tmp2 = tmp2 %>% mutate(subgroup=case_when(
  subgroup==0 ~ "S+",
  subgroup==1 ~ "S-",
  TRUE ~ "equivocal"
))
(tsp$plotIncidenceQuantiles(tmp2 %>% filter(date>"2021-03-01" & ageCat !="unknown"),colour = subgroup,denominatorExpr = population/1000000)+facet_wrap(vars(ageCat))+ylab("daily incidence per 1M")) %>%
  standardPrintOutput::saveThirdPageFigure(here::here("output/BoltonIncidencePer1M"))
tmp3 = tmp2 %>% tsp$aggregateAge() %>% tsp$estimateRt(window = 14,quick = FALSE)
tsp$plotRt(tmp3,colour = subgroup,dates =  "2021-03-01", window=14, rtlim = c(0,3)) %>%
  standardPrintOutput::saveThirdPageFigure(here::here("output/BoltonRt"))

ggplot(
  sgllLoc %>% #filter(LTLA_code %in% areasOfConcern$`NW Cluster`) %>% 
    filter(specimen_date>"2021-04-01") %>% mutate(sgene=case_when(
      sgtf_under30CT==0 ~ "S+",
      sgtf_under30CT==1 ~ "S-",
      TRUE ~ "equivocal"
    )) ,aes(x=age, colour=sgene))+geom_histogram(binwidth = 1)+facet_wrap(vars(sgene), scales = "free")

```
```{r}



```

```{r}

ltlaSgeneCountsAll = sgllLoc %>% dpc$spim$getLineListIncidence(ageBreaks = NULL,codeTypes = "LAD", subgroup = sGene)

ltlaToAoc = function(l) l %>% 
  dpc$demog$findDemographics() %>% 
  inner_join(aocDf, by = c("code"="ltlaCode")) %>% 
  select(-code,-name,-name.original) %>%
  rename(name = area) %>%
  mutate(
    name = as.character(name),
    code = name,
    codeType = "AOC"
  ) %>%
  ukcovidtools::covidStandardDateGrouping() %>%
  summarise(value = sum(value), population = sum(population)) %>%
  ungroup()

aocSgeneCounts = bind_rows(
  ltlaToAoc(ltlaSgeneCounts),
  symptomaticPillar2CtrySgeneCounts %>% dpc$demog$findDemographics()
)

aocSgeneCounts %>% filter(subgroup!="equivocal") %>% readr::write_csv(here::here("output/symptomatic-pillar2-cases-by-aoc-and-sgene.csv"))

aocSgeneCountsAll = bind_rows(
  ltlaToAoc(ltlaSgeneCountsAll),
  allPillar2CtrySgeneCounts %>% dpc$demog$findDemographics()
)

aocSgeneCountsAll %>% filter(subgroup!="equivocal") %>% readr::write_csv(here::here("output/all-pillar2-cases-by-aoc-and-sgene.csv"))

aocDf %>% readr::write_csv(here::here("output/ltla-to-aoc.csv"))



tsp$plotIncidenceQuantiles(aocSgeneCounts, colour=subgroup, dates="2021-03-01")+facet_wrap(vars(name))+scale_y_continuous(trans="log1p", breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000))
tsp$plotIncidenceQuantiles(aocSgeneCountsAll, colour=subgroup, dates="2021-03-01")+facet_wrap(vars(name))+scale_y_continuous(trans="log1p", breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000))



clusters = list(
  "CommunityIndia2" = c("NW Cluster","Bedford etc","Sefton & Lpl"),
  "MixedIndia2" = c("Leicester","Nottingham","West London"),
  "MixedSA" = c("East London"),
  "Reference" = c("England")
)

aocSgeneCountsGrowth = aocSgeneCounts %>% tsp$estimateGrowthRate(window = 56,growthRateWindow = 7)

lapply(names(clusters),function(clustName) {
  (
    tsp$plotGrowthRate(
      aocSgeneCountsGrowth %>% filter(name %in% clusters[[clustName]]) %>% filter(subgroup!="equivocal"),
      colour = subgroup,dates="2021-03-01",
      rlim=c(-0.25,0.25),
      growthVar = Growth.poisson,
      growthSEVar = Growth.SE.poisson,
    )+facet_wrap(vars(name))+standardPrintOutput::watermark(lab="provisional")+ylab("r - sliding 8 week window")
  ) %>%   
    standardPrintOutput::saveThirdPageFigure(paste0(here::here("output/sGeneAOCHighGrowth-"),clustName))
})

aocSgeneCountsRt = aocSgeneCounts %>% tsp$estimateRt(window = 14)
lapply(names(clusters),function(clustName) {
  (
    tsp$plotRt(
      aocSgeneCountsRt %>% filter(name %in% clusters[[clustName]]) %>% filter(subgroup!="equivocal"),
      colour = subgroup,
      dates="2021-03-01",
      rtlim=c(0,3.5)
    )+
      facet_wrap(vars(name))+
      standardPrintOutput::watermark(lab="provisional")) %>%   
    standardPrintOutput::saveThirdPageFigure(paste0(here::here("output/sGeneAOCHighRt-"),clustName))
})




# aocSgeneCountsLogGrowth = aocSgeneCounts %>% tsp$logIncidenceStats(smoothingWindow = 21,growthRateWindow = 7, leftSided=TRUE)
# (tsp$plotGrowthRate(aocSgeneCountsLogGrowth %>% 
#               filter(name %in% c("NW Cluster","Bedford etc","Sefton & Lpl")) %>% 
#               filter(subgroup!="equivocal"),colour = subgroup,dates="2021-03-01",rlim=c(-0.25,0.25))+facet_wrap(vars(name))+standardPrintOutput::watermark(lab="provisional")) %>%   
#   standardPrintOutput::saveThirdPageFigure(here::here("output/sGeneAOCHighGrowthLog"))
# 
# 
# aocSgeneCountsRtAll = aocSgeneCountsAll %>% tsp$estimateRt(window = 14)
# (tsp$plotRt(aocSgeneCountsRt %>% filter(subgroup!="equivocal"),colour = subgroup,dates="2021-03-01",rtlim=c(0,3.5))+facet_wrap(vars(name))+standardPrintOutput::watermark(lab="provisional")) %>%   
#   standardPrintOutput::saveHalfPageFigure(here::here("output/sGeneAllAOCRt"))


#ggplot(ltlaSgeneCounts, aes(x=date,y=value))
# ltlasSgeneCounts is asymptomatic == "N"
#ltlaSgeneCounts2 = ltlaSgeneCounts %>% filter(subgroup != "equivocal") %>% tsp$logIncidenceStats(growthRateWindow = 7,nocache=TRUE)
ltlaSgeneCounts2 = ltlaSgeneCounts %>% filter(subgroup != "equivocal") %>% tsp$estimateGrowthRate(window = 14, growthRateWindow = 7)
# ltlaSgeneCounts2 = ltlaSgeneCounts2 %>% rename(
#   Growth.windowed.value = Growth.windowed.poisson,
#   Growth.windowed.SE.value = Growth.windowed.SE.poisson
#   ) #etc.

#ltlaSgeneCounts2 = ltlaSgeneCounts %>% filter(subgroup != "equivocal") %>% tsp$estimateGrowthRate(window = 14, growthRateWindow = 7)
```

```{r}
tmp = ltlaSgeneCounts2 %>% select(name,code,date,subgroup,starts_with("Growth"))

tmp2 = tmp %>% 
  filter(subgroup=="S+") %>% inner_join(
    tmp %>% filter(subgroup=="S-"),
    by = c("code","name","date"),
    suffix = c(".pos",".neg")
  ) %>% 
  mutate(
    alpha = 1-sqrt(as.numeric(maxDate-date)/as.numeric(maxDate-minDate))
  ) %>%
  mutate(
    # assume independence between S- and S+ growth rates and that CIs are normally distributes
    Growth.windowed.delta = Growth.windowed.value.pos - Growth.windowed.value.neg,
    Growth.windowed.SE.delta = sqrt(Growth.windowed.SE.value.pos^2+Growth.windowed.SE.value.neg^2),
    Growth.windowed.ProbPos.delta = 1-pnorm(0,mean=Growth.windowed.delta,sd=Growth.windowed.SE.delta),
    Growth.delta = Growth.value.pos - Growth.value.neg,
    Growth.SE.delta = sqrt(Growth.SE.value.pos^2+Growth.SE.value.neg^2),
    Growth.ProbPos.delta = 1-pnorm(0,mean=Growth.delta,sd=Growth.SE.delta)
  ) %>%
  rename(
    `S- Growth Rate` = Growth.value.neg,
    `S+ Growth Rate` = Growth.value.pos
  )

# thanks stackoverflow: https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
differenceStat <- function(m1,m2,s1,s2,n1,n2)
{
  se <- sqrt( (s1^2/n1) + (s2^2/n2) )
  # welch-satterthwaite df
  df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  t <- (m1-m2)/se 
  # difference of means * (1-p-value)
  out = (m1-m2) * (1-2*pt(-abs(t),df))    
  return(out) 
}

diagonals= tibble(
  intercept = c(0.2,0.1,0,-0.1,-0.2),
  slope = 1,
  label = forcats::as_factor(c("-20%","-10%","same","+10%","+20%"))
)


maxDate = max(tmp2$date)-4
minDate = maxDate-28

topNgrowth = tmp2 %>% filter(date==maxDate) %>% mutate(
  diffStat = differenceStat(
    Growth.windowed.value.pos,
    Growth.windowed.value.neg,
    Growth.windowed.SE.value.pos,
    Growth.windowed.SE.value.neg,
    Growth.windowed.Window.value.pos,
    Growth.windowed.Window.value.neg
  )) %>%
  arrange(desc(diffStat)) %>%
  filter(row_number() <= 10) #%>%
#select(name,code)

# topNgrowth = tmp2 %>% filter(date==maxDate) %>% mutate(
#   #diffStat = Growth.windowed.poisson.pos-Growth.windowed.poisson.neg
#     diffStat = differenceStat(
#       Growth.windowed.poisson.pos,
#       Growth.windowed.poisson.neg,
#       Growth.windowed.SE.poisson.pos,
#       Growth.windowed.SE.poisson.neg,
#       14,14
#     )
#   ) %>% 
#   ungroup() %>%
#   arrange(desc(diffStat)) %>% 
#   filter(row_number() <= 10) %>% 
#   select(name,code)

pathData = tmp2 %>% semi_join(topNgrowth, by="code") %>% 
  filter(date>minDate & date<=maxDate) %>% 
  
  #pointData = tmp2 %>% filter(date==maxDate & abs(Growth.windowed.value.pos) >0.01) %>% mutate(alpha = 1)
  labelData = tmp2 %>% 
  semi_join(topNgrowth, by="code") %>% 
  filter(date==maxDate) %>% 
  mutate(alpha = 1) %>%
  rename(
    `S- Growth Rate` = Growth.value.neg,
    `S+ Growth Rate` = Growth.value.pos
  )

ggplot(pathData, aes(x=`S- Growth Rate`, y=`S+ Growth Rate`, group=name, alpha=alpha))+
  coord_fixed(xlim=c(-0.2,0.2),ylim=c(-0.2,0.2))+
  geom_hline(yintercept = 0,colour="grey40")+
  geom_vline(xintercept = 0,colour="grey40")+
  geom_abline(data = diagonals, aes(intercept = intercept,slope=slope, colour=label), inherit.aes=FALSE)+
  #geom_abline(slope = 1, intercept = 0.1, colour="#FFB0B0")+
  #geom_abline(slope = 1, intercept = -0.1, colour="#B0FFB0")+
  #geom_abline(slope = 1, intercept = 0.2, colour="#FF8080")+
  #geom_abline(slope = 1, intercept = -0.2, colour="#B0FFB0")+
  #geom_point(data=pointData,colour="#B0B0FF")+
  ggrepel::geom_text_repel(data=labelData,aes(label=name),colour="grey40",min.segment.length = 0,box.padding = 2, 
                           max.overlaps = Inf)+
  geom_path(colour="blue")+
  geom_point(data=labelData,colour="blue")+
  
  guides(alpha="none")+
  xlab("S- growth rate")+
  ylab("S+ growth rate")

# ggplot(pathData, aes(x=Growth.windowed.value.neg, y=Growth.windowed.value.pos, group=name, alpha=alpha))+
#   coord_fixed(xlim=c(-0.15,0.15),ylim=c(-0.15,0.15))+
#   geom_hline(yintercept = 0,colour="grey40")+
#   geom_vline(xintercept = 0,colour="grey40")+
#   geom_abline(data = diagonals, aes(intercept = intercept,slope=slope, colour=label))+
#   #geom_abline(slope = 1, intercept = 0.1, colour="#FFB0B0")+
#   #geom_abline(slope = 1, intercept = -0.1, colour="#B0FFB0")+
#   #geom_abline(slope = 1, intercept = 0.2, colour="#FF8080")+
#   #geom_abline(slope = 1, intercept = -0.2, colour="#B0FFB0")+
#   #geom_point(data=pointData,colour="#B0B0FF")+
#   ggrepel::geom_text_repel(data=labelData,aes(label=name),colour="grey40",min.segment.length = 0,box.padding = 2, 
#                            max.overlaps = Inf)+
#   geom_path(colour="blue")+
#   geom_point(data=labelData,colour="blue")+
#   
#   guides(alpha="none")+
#   xlab("S- growth rate")+
#   ylab("S+ growth rate")

topNgrowth = tmp2 %>% filter(date==maxDate) %>%
  arrange(desc(Growth.windowed.ProbPos.delta*Growth.windowed.delta)) %>%
  filter(row_number() <= 10) %>%
  select(name,code)

(ggplot(tmp2 %>% filter(date>"2020-10-01"),
        aes(
          x=date,
          y=Growth.windowed.delta,
          group=name#,
          #colour=Growth.windowed.ProbPos.delta
        )
)+geom_line(alpha=0.2)+geom_hline(yintercept = 0,colour="red"))  %>%
  standardPrintOutput::saveHalfPageFigure(paste0(here::here("output/AllGrowthDelta")))

(ggplot(tmp2 %>% semi_join(topNgrowth) %>% filter(date>"2020-02-01"),
        aes(x=date,
            y=Growth.windowed.delta,
            ymin=Growth.windowed.delta-1.96*Growth.windowed.SE.delta,
            ymax=Growth.windowed.delta+1.96*Growth.windowed.SE.delta,
            group=name)
)+
    geom_ribbon(alpha=0.2,fill="grey50")+
    geom_line(mapping=aes(colour=name)))  %>%
  standardPrintOutput::saveHalfPageFigure(paste0(here::here("output/Top10GrowthDelta")))

rollingPercent = allVoc %>% #combinedSpositives %>% 
  filter(!is.na(type)) %>%
  group_by(LTLA_name, type, date=specimen_date) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  tidyr::complete(LTLA_name,type,date,fill=list(count=0)) %>%
  group_by(LTLA_name, type) %>% 
  arrange(date) %>%
  mutate(Roll.count = stats::filter(count,filter=rep(1,28)/28,sides=1)) %>%
  filter(!is.na(Roll.count)) %>% 
  group_by(LTLA_name,date) %>% 
  mutate(binom::binom.confint(Roll.count, sum(Roll.count), conf.level = 0.95, methods = "wilson")) %>% 
  mutate(weight = sum(Roll.count)) %>% 
  ungroup()

growthVPercent = tmp2 %>% inner_join(rollingPercent, by = c("name"="LTLA_name","date"="date")) %>%
  mutate(growthSigHigher = ifelse(Growth.windowed.ProbPos.delta>0.0975,"sig","ns"))

(ggplot(growthVPercent %>% filter(date==max(date)), aes(x=mean,y=Growth.windowed.delta,size=weight))+geom_point(alpha=0.3)+facet_wrap(vars(type),scales="free")+
    #geom_smooth(method = "lm",fullrange=TRUE)+
    coord_cartesian(xlim = c(0,1))) %>%
  standardPrintOutput::saveHalfPageFigure(paste0(here::here("output/SGeneGrowthDeltaVsPercentSequencing")))

(ggplot(growthVPercent %>% filter(date==max(date)), aes(x=Roll.count,y=Growth.windowed.delta,size=weight))+geom_point(alpha=0.3)+facet_wrap(vars(type))
  #geom_smooth(method = "lm",fullrange=TRUE)+
) %>%
  standardPrintOutput::saveHalfPageFigure(paste0(here::here("output/SGeneGrowthDeltaVsCountSequencing")))

ggplot(growthVPercent %>% filter(type %in% c("B.1.351","B.1.617.1","B.1.617.2","P.1 & P.2") & date==maxDate), aes(x=type,colour=growthSigHigher,y=mean))+geom_boxplot()+scale_y_log10()

```


```{r}
# map = dpc$geog$getMap("LAD19")
# 
# map2 = map %>% inner_join(rollingPercent, by=c("name"="LTLA_name"))
# ggplot(map2 %>% filter(date==maxDate))+geom_sf(aes(fill=mean))+facet_wrap(vars(type))+scale_fill_viridis_c()
#                    
```

```{r}

map3 = map %>% inner_join(tmp2, by="code")
(ggplot(map3 %>% filter(date==maxDate))+geom_sf(aes(fill=Growth.windowed.delta, colour=ifelse(Growth.windowed.ProbPos.delta>0.0975,"sig","ns")),size=0.3)+scale_color_manual(values = c("sig"="white","ns"="black"))+scale_fill_gradient2(high = "yellow", mid="black",low="cyan",limits=c(-0.15,0.15),oob=scales::squish)) %>%
  standardPrintOutput::saveFullPageFigure(paste0(here::here("output/SgeneGrowthDelta-England")))

lapply(names(locations), function(name) {
  limits = locations[[name]]
  xlim = c(limits[2],limits[4])
  ylim = c(limits[1],limits[3])
  (ggplot(map3 %>% filter(date==max(date)))+geom_sf(aes(fill=Growth.windowed.delta))+scale_fill_gradient2(high = "orange",low="blue",limits=c(-0.15,0.15),oob=scales::squish)+coord_sf(xlim=xlim,ylim=ylim)) %>%
    standardPrintOutput::saveHalfPageFigure(paste0(here::here("output/SgeneGrowthDelta-"),name))
})

```

```{r}
topNgrowth %>% inner_join(rollingPercent %>% filter(date==maxDate) %>% mutate(percent =sprintf("%1.1f%%",mean*100)) %>% select(LTLA_name,type,percent)  %>% pivot_wider(values_from = percent, names_from=type), by = c("name"="LTLA_name"))
```


```{r}
ltlaSgeneCounts3 = ltlaSgeneCounts2 %>% dpc$demog$findDemographics()
ltlaSgeneCounts3 = ltlaSgeneCounts3 %>% inner_join(rollingPercent %>% rename(variant=type), by = c("name"="LTLA_name","date"="date"))


tsp$plotGrowthIncidence(ltlaSgeneCounts3 %>% filter(variant == "B.1.617.2") %>% group_by(name, variant,subgroup),plotDates = maxDate,timespan = 28,colour = mean*weight,populationExpr = population,maxSize = 0.6)+facet_wrap(vars(variant,subgroup))+scale_colour_viridis_b()

```