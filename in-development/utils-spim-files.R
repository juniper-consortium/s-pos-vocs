#' 
#' 
#' 
#' augmentLineListWithLSOA = function(ll, ltlaCodeCol = "LTLA_code", imdRankCol="imd_rank") {
#'   ltlaCodeCol = ensym(ltlaCodeCol)
#'   imdRankCol = ensym(imdRankCol)
#'   imd = self$demog$getIMDData() %>% select(
#'     !!ltlaCodeCol := `Local Authority District code (2019)`,
#'     !!imdRankCol := `Index of Multiple Deprivation (IMD) Rank (where 1 is most deprived)`,
#'     LSOA_code = `LSOA code (2011)`,
#'     LSOA_name = `LSOA name (2011)`
#'   )
#'   return(ll %>% left_join(imd, by=c(as_label(ltlaCodeCol),as_label(imdRankCol))))
#' }
#' 
#' getVAMLineList = function(...) {
#'   message("Using: ",path)
#'   tmp2 = readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()))
#'   datecols = colnames(tmp2)[colnames(tmp2) %>% stringr::str_detect("date|_at")]
#'   for(datecol in datecols) {
#'     tmp2[[datecol]] = suppressWarnings(as.Date(tmp2[[datecol]],"%Y-%m-%d"))
#'   }
#'   idcols = colnames(tmp2)[colnames(tmp2) %>% stringr::str_detect("id")]
#'   for(idcol in idcols) {
#'     tmp2[[idcol]] = suppressWarnings(as.integer(tmp2[[idcol]]))
#'   }
#'   tmp2 = tmp2 %>% mutate(
#'     age = suppressWarnings(as.integer(age))
#'   )
#'   return(tmp2)
#' }
#' 
#' getCTASLineList = function(...) {
#'   #/home/terminological/S3/encrypted/2021-03-29/20210329 CTAS SGTF data.zip
#'   path = self$getLatest(self$filter$ctasLineList)
#'   message("Using: ",path)
#'   self$getSaved("CTAS", params = list(path), ...,  orElse = function (...) {
#'     tmp = self$fileProvider$getFile(path)
#'     zipPath = fs::path_file(path) %>% stringr::str_replace("\\.zip",".csv")
#'     tmp2 = readr::read_csv(unz(tmp,filename=zipPath), col_types = readr::cols(.default = readr::col_character()))
#'     # tmp2 = tmp2 %>% mutate(genomic_specimen_date = suppressWarnings(as.Date(genomic_specimen_date,"%Y%m%d")))
#'     datecols = colnames(tmp2)[colnames(tmp2) %>% stringr::str_detect("date|_at")]
#'     for(datecol in datecols) {
#'       tmp2[[datecol]] = suppressWarnings(as.Date(tmp2[[datecol]],"%Y-%m-%d"))
#'     }
#'     idcols = colnames(tmp2)[colnames(tmp2) %>% stringr::str_detect("id")]
#'     for(idcol in idcols) {
#'       tmp2[[idcol]] = suppressWarnings(as.integer(tmp2[[idcol]]))
#'     }
#'     tmp2 = tmp2 %>% mutate(
#'       completed = as.logical(completed),
#'       sex = self$normaliseGender(sex),
#'       sgtf = as.integer(sgtf),
#'       sgtf_under30ct = as.integer(sgtf_under30ct),
#'       p2ch1cq = as.double(p2ch1cq),
#'       p2ch2cq = as.double(p2ch2cq),
#'       p2ch3cq = as.double(p2ch3cq),
#'       p2ch4cq = as.double(p2ch4cq),
#'       age = as.integer(age)
#'     )
#'     if(file.exists(tmp)) unlink(tmp)
#'     return(tmp2)
#'   })
#' }
#' 
#' 
#' #' @description Load line list
#' #' 
#' #' @return raw line list data set
#' getSGeneLineList = function(...) {
#'   path = self$getLatest(self$filter$sgene)
#'   message("Using: ",path)
#'   self$getSaved("SGENE-LINE-LIST", params = list(path), ..., orElse = function (...) {
#'     if (stringr::str_detect(path,"zip")) {
#'       tmpFile = self$fileProvider$getFile(path)
#'       zipPath = fs::path_file(path) %>% stringr::str_replace("\\.zip",".csv")
#'       tmp = readr::read_csv(unz(tmpFile, filename=zipPath))
#'     } else {
#'       tmp = readr::read_csv(self$fileProvider$getFile(path))
#'     }
#'     return(tmp %>% dplyr::ungroup())
#'   })
#' },
#' 
#' #' Interpret S gene status according to various cut off values
#' #' function to help interpret S gene CT values in context of N gene and ORF gene to give S gene status. 
#' #' With the defaults this produces the same result as the sgtf_30 column in the source SGTF line list
#' #' Defaults are S:30,ORF:30,N:30,Control:Inf
#' #'
#' #' @param sGeneLineList - a dataframe includeing 
#' #' @param S_CT - S gene detected when P2CH3CQ <= this value
#' #' @param ORF1ab_CT - ORF1ab gene detected when P2CH1CQ <= this value
#' #' @param N_CT - N gene detected when P2CH2CQ <= this value
#' #' @param Control_CT - control sample is positive when P2CH4CQ <= this value
#' #'
#' #' @return - the same dataframe with additional columns including "sGene" and "result"
#' #'
#' #' @examples coxData = coxData %>% interpretSGene()
#' interpretSGene = function(sGeneLineList, S_CT = 30, ORF1ab_CT = 30, N_CT = 30, Control_CT = Inf) {
#'   sGeneLineList %>% 
#'     mutate(
#'       ORF1ab_CT_threshold = ORF1ab_CT,
#'       N_CT_threshold = N_CT,
#'       S_CT_threshold = S_CT,
#'       S_pos = P2CH3CQ > 0 & P2CH3CQ <= S_CT,
#'       S_undetect = P2CH3CQ == 0,
#'       N_pos = P2CH2CQ > 0 & P2CH2CQ <= N_CT,
#'       ORF1ab_pos = P2CH1CQ > 0 & P2CH1CQ <= ORF1ab_CT,
#'       Control_pos = P2CH4CQ > 0 & P2CH4CQ <= Control_CT,
#'       sGene = case_when(
#'         is.na(P2CH1CQ) ~ "Unknown",
#'         S_pos & N_pos & ORF1ab_pos & Control_pos ~ "Positive",
#'         S_undetect & N_pos & ORF1ab_pos & Control_pos ~ "Negative",
#'         TRUE ~ "Equivocal"
#'       ),
#'       CT_N = ifelse(P2CH2CQ > 0, P2CH2CQ, 40)
#'     ) %>% 
#'     mutate(
#'       result = case_when(
#'         is.na(P2CH1CQ) ~ "Unknown",
#'         !Control_pos ~ "No control",
#'         TRUE ~ paste0(ifelse(S_pos,"S+","S-"),ifelse(N_pos,"N+","N-"),ifelse(ORF1ab_pos,"ORF+","ORF-")))
#'     ) %>%
#'     mutate(
#'       sGene = sGene %>% forcats::fct_relevel("Positive"),
#'       relativeCopyNumber = 2^(median(CT_N,na.rm=TRUE)-CT_N)
#'     )
#' }
#' 
#' 
#' getSGeneEras = function(cutoff = 28, ...) {
#'   path = self$getLatest(self$filter$sgene)
#'   self$getSaved("SGENE-ERAS", params = list(path), ..., orElse = function (...) {
#'     sgll = self$getSGeneLineList()
#'     # group by patient and find time delay between tests (where there are more than one)
#'     tmp = sgll %>% arrange(FINALID,specimen_date) %>% mutate(delay = ifelse(FINALID==lag(FINALID), as.numeric(specimen_date - lag(specimen_date)), NA_real_))
#'     # TODO: there is some interesting properties of the delay
#'     # ggplot(tmp, aes(x=delay))+geom_density()+scale_x_continuous(trans="log1p",breaks=c(0,10,20,50,100,200,500,1000))+facet_wrap(vars(sgtf_under30CT))
#'     # ggplot(tmp, aes(x=delay,y=P2CH1CQ))+geom_density_2d()
#'     # apply some heuristics to determine whether a test is part of the same infection or a new one
#'     tmp2 = tmp %>% mutate(era = case_when(
#'       is.na(delay) ~ "new",
#'       is.na(sgtf_under30CT) & delay < cutoff*2 ~ "same", # prolonged recovery
#'       delay < cutoff ~ "same",
#'       TRUE ~ "new"
#'     ))
#'     # assign an eraIndex - essentially the count of novel infection episodes
#'     tmp3 = tmp2 %>% group_by(FINALID) %>% arrange(specimen_date) %>% mutate(eraIndex = cumsum(ifelse(era=="new",1,0)))
#'     # summarise sGene data into a single value for each era
#'     # TODO: each era may have multiple positive tests there is an opportunity to look at the CT values over time and 
#'     # fit some sort of model here
#'     tmp4 = tmp3 %>% group_by(FINALID,eraIndex) %>% 
#'       summarise(
#'         earliest_specimen_date = min(specimen_date,na.rm=TRUE), 
#'         latest_specimen_date = max(specimen_date,na.rm=TRUE), 
#'         tests=n(), 
#'         minSgtf = min(sgtf_under30CT,na.rm = TRUE), 
#'         maxSgtf = max(sgtf_under30CT,na.rm = TRUE),
#'         min_P2CH1CQ = min(P2CH1CQ,na.rm=TRUE),
#'         min_P2CH2CQ = min(P2CH2CQ,na.rm=TRUE),
#'         min_P2CH3CQ = min(P2CH3CQ,na.rm=TRUE),
#'         min_P2CH4CQ = min(P2CH4CQ,na.rm=TRUE)
#'       ) %>% ungroup() %>% mutate(sgtf_under30CT = case_when(
#'         minSgtf == maxSgtf ~ minSgtf,
#'         TRUE ~ NA_real_
#'       ))
#'     return(tmp4)
#'   })
#' }
#' 
#' 
#' 
#' #' @description Load line list
#' #' 
#' #' @return raw line list data set
#' getLineList = function(...) {
#'   path = self$getLatest(self$filter$lineList)
#'   message("Using: ",path)
#'   self$getSaved("LINE-LIST", params = list(path), ..., orElse = function (...) {
#'     if (stringr::str_detect(path,"zip")) {
#'       tmpFile = self$fileProvider$getFile(path)
#'       zipPath = fs::path_file(path) %>% stringr::str_replace("\\.zip",".csv")
#'       tmp = readr::read_csv(unz(tmpFile, filename=zipPath), col_types = readr::cols(.default = readr::col_character()))
#'       tmp = tmp %>% 
#'         dplyr::mutate(
#'           Onsetdate = maybeDMYorMDY(Onsetdate),
#'           specimen_date = maybeDMYorMDY(specimen_date),
#'           lab_report_date = maybeDMYorMDY(lab_report_date)
#'         ) 
#'       
#'     } else if (stringr::str_detect(path,"csv")) {
#'       tmp = readr::read_csv(self$fileProvider$getFile(path), col_types = readr::cols(.default = readr::col_character()))
#'       tmp = tmp %>% 
#'         dplyr::mutate(
#'           Onsetdate = maybeDMYorMDY(Onsetdate),
#'           specimen_date = maybeDMYorMDY(specimen_date),
#'           lab_report_date = maybeDMYorMDY(lab_report_date)
#'         ) 
#'       
#'     } else {
#'       tmp = readxl::read_excel(path.expand(path), 
#'                                col_types = "text") #c("numeric", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "date", "date", "date"))
#'       tmp = tmp %>% 
#'         dplyr::mutate(
#'           Onsetdate = suppressWarnings(as.Date(as.numeric(Onsetdate),"1899-12-30")),
#'           specimen_date = suppressWarnings(as.Date(as.numeric(specimen_date),"1899-12-30")),
#'           lab_report_date = suppressWarnings(as.Date(as.numeric(lab_report_date),"1899-12-30"))
#'         )
#'     }
#'     
#'     if(any(is.na(tmp$specimen_date))) warning("NA sprecimen dates in cases file")
#'     
#'     return(tmp %>% mutate(
#'       pillar_2_testingkit = tolower(pillar_2_testingkit),
#'       age = suppressWarnings(as.numeric(age)),
#'       FINALID = as.numeric(FINALID),
#'       imd_rank = as.integer(imd_rank),
#'       imd_decile = as.integer(imd_decile),
#'       ethnicity_final = case_when(
#'         ethnicity_final %in% c("African (Black or Black British)","Any other Black background","Caribbean (Black or Black British)") ~ "Afro-caribbean",
#'         ethnicity_final %in% c("Any other Asian background","Bangladeshi (Asian or Asian British)","Indian (Asian or Asian British)","Pakistani (Asian or Asian British)") ~ "Asian",
#'         ethnicity_final %in% c("Any other White background","British (White)","Irish (White)") ~ "White",
#'         ethnicity_final %in% c("Any other Mixed background","Any other ethnic group","White and Black Caribbean (Mixed)","White and Black African (Mixed)","Chinese (other ethnic group)") ~ "Other",
#'         TRUE ~ "Unknown"),
#'       residential_category = case_when(
#'         cat == 'Residential dwelling (including houses, flats, sheltered accommodation)' ~ "Residential",
#'         cat == 'Care/Nursing home' ~ "Care home",
#'         cat == 'Undetermined'~"Other/Unknown",
#'         cat == 'Medical facilities (including hospitals and hospices, and mental health)'~"Other/Unknown",
#'         cat == 'Other property classifications'~"Other/Unknown",
#'         cat == 'House in multiple occupancy (HMO)' ~ "Residential",
#'         cat == 'Prisons, detention centres, secure units'~"Other/Unknown",
#'         cat == 'Residential institution (including residential education)'~"Other/Unknown",
#'         cat == 'No fixed abode'~"Other/Unknown",
#'         cat == 'Overseas address'~"Other/Unknown",
#'         TRUE ~ "Other/Unknown"
#'       )
#'     ) %>% dplyr::ungroup())
#'   })
#'   
#'   # TODO: https://github.com/sarahhbellum/NobBS
#' }
#' 
#' 
#' #' @description Load incidence from line list
#' #' 
#' #' @param ageBreaks - a list of ages which form the cut points for breaking continuous ages into ranges (or NULL for a single age category)
#' 
#' #' @return a covidTimeseriesFormat dataframe
#' 
#' getLineListIncidence = function(ll=NULL, ageBreaks = NULL, gender=FALSE, specimenOrReport="specimen", subgroup="pillar", filterExpr=NULL, codeTypes = c("CTRY","NHSER"), truncate=NULL, ...) {
#'   filterExpr = enexpr(filterExpr)
#'   subgroup = tryCatch(ensym(subgroup), error = function(e) NULL)
#'   if(!identical(ll,NULL)) {
#'     tmp = ll
#'   } else {
#'     tmp = self$getLineList(...) 
#'   }
#'   self$getSaved("LINE-LIST-INCIDENCE", params=list(tmp, ageBreaks, specimenOrReport,as_label(subgroup), as_label(filterExpr), codeTypes, gender), ..., orElse = function (...) covidTimeseriesFormat({
#'     if(!identical(filterExpr,NULL)) 
#'       tmp = tmp %>% filter(!!filterExpr)
#'     tmp = tmp %>% dplyr::mutate(ageCat = age %>% self$cutByAge(ageBreaks)) 
#'     if (gender) {
#'       tmp = tmp %>% dplyr::mutate(gender=self$normaliseGender(sex))
#'     } else {
#'       tmp = tmp %>% dplyr::mutate(gender=NA_character_)
#'     }
#'     if(!identical(subgroup,NULL)) {
#'       tmp = tmp %>% dplyr::mutate(subgroup=!!subgroup)
#'     } else {
#'       tmp = tmp %>% dplyr::mutate(subgroup=NA_character_)
#'     }
#'     if(specimenOrReport == "report")
#'       tmp = tmp %>% dplyr::mutate(date = as.Date(lab_report_date))
#'     else
#'       tmp = tmp %>% dplyr::mutate(date = as.Date(specimen_date))
#'     
#'     selectByRegion = function(df, code, codeType, name) {
#'       code = ensym(code)
#'       name = ensym(name)
#'       # check column exists
#'       if(!(as_label(code) %in% colnames(df))) return(tibble())
#'       df = df %>% dplyr::mutate(code = !!code, codeType=codeType, name=!!name) %>% 
#'         dplyr::mutate(
#'           code = ifelse(is.na(code),"E99999999",code),
#'           name = ifelse(is.na(code),"Unknown (England)",name)
#'         ) %>%
#'         dplyr::group_by( code,codeType,name,date, ageCat, gender,subgroup) %>% 
#'         dplyr::summarise(value = n()) 
#'       return(df)
#'     }
#'     
#'     out = NULL
#'     if ("CTRY" %in% codeTypes) {
#'       england = tmp %>% dplyr::mutate(code = "E92000001", codeType= "CTRY", name="England") %>% 
#'         dplyr::group_by(code,codeType,name,date, ageCat, gender,subgroup) %>% 
#'         dplyr::summarise(value = n())
#'       out = out %>% bind_rows(england)
#'     }
#'     
#'     if ("NHSER" %in% codeTypes) {
#'       nhser = tmp %>% selectByRegion(NHSER_code, "NHSER", NHSER_name)
#'       isNhser = nhser %>% self$codes$allPresentAndCorrect(codeTypes=c("NHSER","PSEUDO"))
#'       
#'       if(!isNhser) {
#'         nhser = tmp %>% selectByRegion(NHSER_code, "NHSER19CDH", NHSER_name) %>% 
#'           dplyr::inner_join(
#'             self$codes$getMappings() %>% dplyr::filter(fromCodeType=="NHSER19CDH" & toCodeType=="NHSER"), 
#'             by=c("code"="fromCode")
#'           ) %>%
#'           dplyr::ungroup() %>%
#'           dplyr::select(-code,-codeType, -fromCodeType,-rel,-weight) %>%
#'           dplyr::rename(code = toCode, codeType=toCodeType)
#'       }
#'       out = out %>% bind_rows(nhser)
#'     }
#'     
#'     if ("PHEC" %in% codeTypes) {out = out %>% bind_rows(tmp %>% selectByRegion(PHEC_code, "PHEC", PHEC_name))}
#'     if ("UA" %in% codeTypes) {out = out %>% bind_rows(tmp %>% selectByRegion(UTLA_code, "UA", UTLA_name))}
#'     if ("LAD" %in% codeTypes) {out = out %>% bind_rows(tmp %>% selectByRegion(LTLA_code, "LAD", LTLA_name))}
#'     if ("LSOA" %in% codeTypes) {out = out %>% bind_rows(tmp %>% selectByRegion(LSOA_code, "LSOA", LSOA_name))}
#'     
#'     out = out %>% dplyr::mutate(source="line list",statistic = "case", type="incidence")
#'     out = out %>% self$fixDatesAndNames(truncate)
#'     out = out %>% self$fillAbsent(completeDates=TRUE)
#'     out = out %>% dplyr::ungroup()
#'     return(out)
#'   }))
#' }