here::i_am("s-pos-clusters.R")
source(here::here("common-setup.R"))

# Define areas ----

tmp = sapply(c("North West","North East and Yorkshire","Midlands","East of England","London","South East"), function(regionName) {
  dpc$codes$getTransitiveClosure() %>% 
    filter(fromCodeType=="LAD",toCodeType=="NHSER") %>% 
    dpc$codes$findNamesByCode(codeVar = toCode) %>% 
    filter(name == regionName) %>% 
    pull(fromCode)
})

# manually defined subregions
areasOfConcern = list(
  
  #`Bolton` = c("E08000001","E06000008"),
  `Bolton and Blackburn` = c("E08000001","E06000008"), 
  #`M65 Corridor` = c("E06000008","E07000120","E07000117"), # Blackburn, Hyndburn, Burnley
  #`Manchester` = c("E08000003", "E08000006", "E08000008", "E08000007", "E08000005", "E08000002", "E07000037"), # Machester, Salford, Tameside, Stockport, Rochdale, Bury, High peak
  
  `Nottingham` = c("E06000018","E07000036","E07000173", "E07000176","E07000172"), # Nottingham, Erewash, Gedling, Rushcliffe, Broxtowe
  `Leicester` = c("E06000016","E07000129","E07000135"), # Leicester, Blaby, Oadby & Wigton
  #`Birmingham` = c("E08000025", "E08000031", "E08000028", "E08000027","E08000029","E08000030"), #Birmingham, Wolverhamptom, Sandwell, Dudley, Solihull, Walsall
  
  
  # `M6 Corridor` = c("E07000123","E08000010", "E07000118"),  # Preston, Wigan, Chorley, Could include St. Helens, Warringtom
  # `Sefton & Lpl` = c("E08000014","E08000012","E07000127", "E08000011","E08000015"), #Sefton, Liverpool, West lancashire, Knowlesley, Wirral
  `Bedford` = c("E06000055","E06000042","E06000056"), # Bedford, Milton Keynes, Contral bedfordshire
  
  `East London` = c("E09000025","E09000002","E09000031","E09000030","E09000012","E09000026"), # Newham,Barking,Waltham Forest,Tower Hamlets, Hackney, Redbridge
  `West London` = c("E09000018","E09000015","E09000005","E09000017","E09000009","E09000018") # Harrow, Brent, Hillingdon, Ealing, Hounslow
  # `Southampton, Portsmouth & IoW` = c("E06000046","E06000045","E07000091","E06000044","E07000086")
) %>% c(tmp)

localClusters = list(
  # "NHS Regions (N)" = c("North West","North East","Midlands"),
  # "NHS Regions (S)" = c("East of England","London","South East"),
  "North West and Midlands" = c("Bolton and Blackburn","Nottingham","Leicester"), #"M6 Corridor","Sefton & Lpl"
  #"Midlands" = c("Nottingham","Leicester","Birmingham"), #"Bedford etc",
  "London and South East" = c("Bedford","East London","West London")
)

nationalClusters = list(
  "Reference" = c("England"),
  "NHS Regions (N)" = c("North West","North East and Yorkshire","Midlands"),
  "NHS Regions (S)" = c("East of England","London","South East")
)

aocDf = lapply(names(areasOfConcern),function(name) tibble(area = rep(name,length(areasOfConcern[[name]])), code = areasOfConcern[[name]])) %>% bind_rows() %>% mutate(area = area %>% ordered(names(areasOfConcern)))
aocDf = aocDf %>% dpc$codes$findNamesByCode(outputNameVar = ltlaName) %>% mutate(ageCat=NA,gender=NA) %>% dpc$demog$findDemographics(codeVar = code)

## Generate map of areas of concern ----

ladMap = dpc$geog$getMap("LAD19")
aocMap = ladMap %>% 
  filter(code %>% stringr::str_starts("E")) %>%
  left_join(
    aocDf %>% select(areaName = area,code), by="code") %>%
  mutate(areaName = ifelse(is.na(areaName),"other",as.character(areaName))) %>%
  group_by(areaName) %>%
  summarise() %>%
  filter(!sf::st_is_empty(.))
