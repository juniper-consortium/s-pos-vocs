#!/usr/bin/Rscript
here::i_am("s-pos-data.R")
source(here::here("common-setup.R"))
source(here::here("s-pos-clusters.R"))

# prevent rerunning this by accident
if (!exists("combinedSpositives")) {

  # Utility function ----
  mapVOCtoLineage = function(df, vocVar = "vam", lineageVar = "type", lineages=c("B.1.1.7","B.1.525","AV.1","B.1.617.1","B.1.617.2","B.1.351","P.1 & P.2"), includeInconclusive = FALSE) {
    vocVar = ensym(vocVar)
    lineageVar = ensym(lineageVar)
    tmp = df %>% mutate(
      !!lineageVar := case_when(
        "B.1.1.7" %in% lineages & !!vocVar %>% stringr::str_detect("20DEC-01") ~ "B.1.1.7",
        "B.1.351" %in% lineages & !!vocVar %>% stringr::str_detect("20DEC-02") ~ "B.1.351",
        "B.1.525" %in% lineages & !!vocVar %>% stringr::str_detect("21FEB-03") ~ "B.1.525",
        "P.1 & P.2" %in% lineages & !!vocVar %>% stringr::str_detect("21JAN-02") ~ "P.1 & P.2",
        "P.1 & P.2" %in% lineages & !!vocVar %>% stringr::str_detect("21JAN-01") ~ "P.1 & P.2",
        "B.1.1.318" %in% lineages & !!vocVar %>% stringr::str_detect("21FEB-04") ~ "B.1.1.318",
        "B.1.617.1" %in% lineages & !!vocVar %>% stringr::str_detect("21APR-01") ~ "B.1.617.1",
        "B.1.617.2" %in% lineages & !!vocVar %>% stringr::str_detect("21APR-02") ~ "B.1.617.2",
        "B.1.617.3" %in% lineages & !!vocVar %>% stringr::str_detect("21APR-03") ~ "B.1.617.3",
        "AV.1" %in% lineages & !!vocVar %>% stringr::str_detect("21MAY-01") ~ "AV.1",
        "C.36.3" %in% lineages & !!vocVar %>% stringr::str_detect("21MAY-02") ~ "C.36.3",
        !!vocVar %>% stringr::str_starts("VOC|VUI") ~ "other VOC/VUI",
        !!vocVar %>% stringr::str_detect("484") ~ "other VOC/VUI",
        !!vocVar %>% stringr::str_detect(stringr::coll("unclassified",ignore_case = TRUE)) ~ "non VOC/VUI",
        #!!vocVar %>% stringr::str_detect(stringr::coll("undetermined",ignore_case = TRUE)) ~ "inconclusive",
        #!!vocVar %>% stringr::str_detect(stringr::coll("inconclusive",ignore_case = TRUE)) ~ "inconclusive",
        TRUE ~ "inconclusive"
      )) 
    if (!includeInconclusive) {
      tmp = tmp %>% filter(!!lineageVar!="inconclusive") %>% mutate(!!lineageVar := !!lineageVar %>% ordered(c("non VOC/VUI",lineages,"other VOC/VUI")))
    } else {
      tmp = tmp %>% mutate(!!lineageVar := !!lineageVar %>% ordered(c("non VOC/VUI",lineages,"other VOC/VUI","inconclusive")))
    }
    return(tmp)
  }
  
  ## Initialise data ----
  
  options("ukcovid.reproduce.at"=as.Date("2021-07-06"))
  # hard date limits defined here:
  maxDate = as.Date("2021-07-01")
  minDate = as.Date("2020-10-01")
  
  loader = dpc$spim
  
  # dpc$spim$getLatestRawFile(dpc$spim$filter$vamLineList,to = "~/Data/s-pos-voc/")
  # dpc$spim$getLatestRawFile(dpc$spim$filter$ctasLineList,to = "~/Data/s-pos-voc/")
  # dpc$spim$getLatestRawFile(dpc$spim$filter$lineList,to = "~/Data/s-pos-voc/")
  # dpc$spim$getLatestRawFile(dpc$spim$filter$sgene,to = "~/Data/s-pos-voc/")
  
  # CTAS line list
  ctasPath = loader$getLatest(dpc$spim$filter$ctasLineList)
  message("Using: ",ctasPath)
  
  ctasGenomics = dpc$getSaved("CTAS-GENOMICS",params = ctasPath, orElse = function(...) {
    ctas = loader$getCTASLineList()
    out = ctas %>% filter(!is.na(genomic_variant)) %>% select(c(starts_with("genomic"))) %>% rename_with(.fn=function(x) stringr::str_remove(x,"genomic_"))  %>% distinct() #,sgtf_finalid, sgtf_under30ct, p2ch3cq, p2ch1cq, p2ch2cq, p2ch4cq))
    rm(ctas)
    return(out)
  })
  
  # VAM line list
  vam = loader$getVAMLineList()
  
  genomics = bind_rows(
    vam %>% select(finalid, variant, exposure_type, specimen_date=specimen_date_sk),
    ctasGenomics %>% anti_join(vam, by="finalid") %>% select(finalid, variant, exposure_type=exp_type,specimen_date)
  ) %>% filter(specimen_date < maxDate & specimen_date >= minDate)
  
  # Cases line list
  ll = loader$getLineList() %>% filter(specimen_date < maxDate & specimen_date >= minDate)
  
  # S gene line list
  sgll = loader$getSGeneLineList() %>% filter(specimen_date < maxDate & specimen_date >= minDate)
  
  dataSources = c(
    ctasPath,
    attr(vam,"paths"),
    attr(ll,"paths"),
    attr(sgll,"paths")
  )
  
  ## Calculate infection episodes ----
  
  # Split line list out into normalised parts.
  llDemog = ll %>% select(FINALID, NHSER_code, NHSER_name, PHEC_code, PHEC_name, UTLA_code, UTLA_name, LTLA_code, LTLA_name, sex, age, ethnicity_final, imd_decile, imd_rank, residential_category, cat) %>% distinct()
  llEpisode = ll %>% select(FINALID, specimen_date, asymptomatic_indicator, pillar, lab_report_date, pillar_2_testingkit, testcentreid, case_category) %>% mutate(episode_type="first positive",episode=1) %>% distinct()
  llTest = ll %>% select(FINALID, specimen_date, case_category,asymptomatic_indicator) %>% mutate(linelist = TRUE)
  
  # Split s gene line list out into normalised parts.
  sgllTest = sgll %>% 
    select(FINALID, specimen_date, CDR_Specimen_Request_SK, sgtf_under30CT) %>% 
    mutate(sglinelist = TRUE)
  
  # combine tests from ll (first only) with tests from S-gene files
  # N.B. this will be a bit of a funny mixture - multiple tests but only for cases that went through taqpath assay.
  tests = llTest %>% 
    full_join(
      sgllTest,
      by = c("FINALID","specimen_date")  
    )  %>%
    mutate(
      linelist = ifelse(is.na(linelist),FALSE,linelist),
      sglinelist = ifelse(is.na(sglinelist),FALSE,sglinelist),
      case_category = ifelse(sglinelist,"TaqPath",case_category)
    ) %>% mutate(
      sGene = case_when(
        is.na(sgtf_under30CT) & sglinelist ~ "equivocal",
        is.na(sgtf_under30CT) & linelist ~ "unknown",
        sgtf_under30CT == 1 ~ "negative",
        sgtf_under30CT == 0 ~ "positive"
      )
    )
  
  # this defines how long between tests before two tests are regarded as a new episode.
  # if the tests are sgtf equivocal double this is allowed.
  cutoff = 28
  
  ## calculate the individual episodes of covid resulting from runs of sequential positive tests <28 days apart.
  allEpisodes = dpc$getSaved(id = "EPISODES",params = list(cutoff, tests), orElse = function(...) {
    
    # use data.table to speed it all up.
    tests2 = dtplyr::lazy_dt(tests)
    
    # look for 
    tests3 = tests2 %>% arrange(FINALID,specimen_date) %>% 
      mutate(
        delay = ifelse(FINALID == lag(FINALID), as.numeric(specimen_date - lag(specimen_date)), NA_integer_)
      ) %>%
      mutate(era = case_when(
      is.na(delay) ~ "new",
      sGene == "equivocal" & delay < cutoff*2 ~ "same", # prolonged recovery
      delay < cutoff ~ "same",
      TRUE ~ "new"
    )) %>%
      # assign an eraIndex - essentially the count of novel infection episodes
      group_by(FINALID) %>% arrange(specimen_date) %>% mutate(eraIndex = cumsum(ifelse(era=="new",1,0)))
      # summarise sGene data into a single value for each era
    
    # TODO: each era may have multiple positive tests there is an opportunity to look at the CT values over time and 
    # fit some sort of model here
    tests4 = tests3 %>% group_by(FINALID,eraIndex) %>% 
      summarise(
        earliest_specimen_date = min(specimen_date,na.rm=TRUE), 
        latest_specimen_date = max(specimen_date,na.rm=TRUE), 
        tests=n(),
        anyPosSGene = any(sGene == "positive"),
        anyNegSGene = any(sGene == "negative"),
        anyEquivSGene = any(sGene == "equivocal"),
        anyUnknSGene = any(sGene == "unknown"),
        asymptomatic_indicator = first(na.omit(asymptomatic_indicator),default="U"),
        lft_only = all(na.omit(case_category=="LFT_Only"))
      ) 
    tests5 = tests4 %>% mutate(
      sGene=case_when(
        anyPosSGene & !anyNegSGene ~ "positive",
        anyNegSGene & !anyPosSGene ~ "negative",
        anyEquivSGene ~ "equivocal",
        TRUE ~ "unknown"
    )) %>% as_tibble()
    
    return(tests5 %>% select(-anyPosSGene, -anyNegSGene, -anyEquivSGene, -anyUnknSGene))
  })
  
  ## Link genomics to episodes ----
  
  # this will lose a few sequenced cases which are not able to be linked to the episodes
  allVoc = genomics %>% 
    left_join(allEpisodes %>% select(FINALID, earliest_specimen_date,latest_specimen_date,eraIndex,sGene), by=c("finalid"="FINALID")) %>%
    mutate(
      # determine smallest time to each era
      offset = pmin(abs(as.numeric(specimen_date-earliest_specimen_date)),abs(as.numeric(specimen_date-latest_specimen_date)))
    ) %>%
    # for every era get sequence with smallest distance to it.
    group_by(finalid,eraIndex) %>% arrange(offset) %>% filter(row_number() == 1) %>% ungroup() %>%
    # make sure the sequence is within 28 days of the era
    filter(is.na(specimen_date) | (specimen_date>=earliest_specimen_date-28 & specimen_date<=latest_specimen_date+28)) %>%
    group_by(finalid) %>%
    filter(!is.na(specimen_date) | eraIndex == max(eraIndex)) %>%
    rename(episode_date = earliest_specimen_date) %>%
    select(-latest_specimen_date,-offset) %>%
    ungroup()
  
  ## Create filtered input broken down by genomics for rest of scripts ----
  
  # filter to remove S gene positives before start of sequencing data
  
  createClassifierInput = function(from = minDate, genomicFilterExpr = NULL, sgeneFilterExpr = NULL, lineages = c("AV.1","B.1.617.2","B.1.351","B.1.1.7"), ...) {
    
    genomicFilterExpr = enexpr(genomicFilterExpr)
    sgeneFilterExpr = enexpr(sgeneFilterExpr)
    if (identical(genomicFilterExpr,NULL)) genomicFilterExpr = TRUE
    if (identical(sgeneFilterExpr,NULL)) sgeneFilterExpr = TRUE
    
    dpc$getSaved(id = "VOC-CLASSIFIER",params = list(from,genomicFilterExpr,sgeneFilterExpr,lineages,allEpisodes,allVoc,mapVOCtoLineage), ..., orElse = function(...) {
      
      #browser()
      
      sgPosUniq = allEpisodes %>% filter(!!sgeneFilterExpr) %>% 
        select(FINALID, eraIndex, sGene, earliest_specimen_date, latest_specimen_date)
    
      tmp = allVoc %>% 
          filter(!!genomicFilterExpr) %>% 
          mapVOCtoLineage(vocVar = variant,lineages = lineages) %>%
          select(finalid = finalid, eraIndex, sGene, type, exposure_type, sequence_date = specimen_date) %>% 
          mutate(
            sequenced=TRUE, 
            exposure_type = case_when(
              exposure_type == "missing_data" ~ NA_character_,
              TRUE ~ stringr::str_remove(exposure_type,", Other")
            )
          )
      
      notTmp = sgPosUniq %>% 
        select(finalid = FINALID, eraIndex, sGene) %>% 
        mutate(type = NA, sequenced = FALSE, exposure_type=NA_character_) %>% 
        anti_join(tmp, by=c("finalid","eraIndex"))
      
      tmp = tmp %>% bind_rows(notTmp) 
      tmp = tmp %>% filter(!is.na(finalid) & !is.na(eraIndex))
      tmp = tmp %>% mutate(
        type = forcats::fct_explicit_na(type,"not sequenced"))
    
      tmp2 = tmp %>% 
        left_join(allEpisodes %>% select(-tests,-sGene), by=c("finalid"="FINALID","eraIndex")) %>%
        filter(earliest_specimen_date >= from)
        
      combinedSpositives = tmp2 %>% as_tibble() %>% 
        inner_join(llDemog, by=c("finalid"="FINALID"))  
      return(combinedSpositives)
    
    }) %>% return()
  
  }
  
  # combinedSpositives = createClassifierInput(from="2021-02-01", genomicFilterExpr = !(variant %in% c("VOC-20DEC-01","VUI-21FEB-03")), sgeneFilterExpr = sGene == "positive", lineages = c("B.1.617.1","B.1.617.2","B.1.351"))
  # combinedSnotNegatives = createClassifierInput(from="2021-02-01", genomicFilterExpr = !(variant %in% c("VOC-20DEC-01","VUI-21FEB-03")), sgeneFilterExpr = sGene != "negative", lineages = c("B.1.617.1","B.1.617.2","B.1.351"))
  # combinedCases = createClassifierInput(from="2021-02-01",lineages = c("B.1.617.1","B.1.617.2","B.1.351","B.1.1.7"))
  # combinedSpositives = createClassifierInput(from="2020-12-31", filterExpr = sGene == "positive", lineages = c("B.1.617.1","B.1.617.2","B.1.351"))
  # combinedCases = createClassifierInput(from="2021-02-01",lineages = c("B.1.617.1","B.1.617.2","B.1.351","B.1.1.7"))
  combinedSpositives = createClassifierInput(from="2021-02-01", genomicFilterExpr = !(variant %in% c("VOC-20DEC-01","VUI-21FEB-03")), sgeneFilterExpr = sGene == "positive", lineages = c("B.1.351","B.1.617.1","B.1.617.2","P.1 & P.2"))
  # combinedSnotNegatives = createClassifierInput(from="2021-02-01", genomicFilterExpr = !(variant %in% c("VOC-20DEC-01","VUI-21FEB-03")), sgeneFilterExpr = sGene != "negative", lineages = c("AV.1","B.1.617.2","B.1.351"))
  combinedCases = createClassifierInput(from="2021-02-01",lineages = c("B.1.351","B.1.617.1","B.1.617.2","P.1 & P.2"))
  
  ## Set important dates for the figures - TODO: need to rename this  ----
  lastDate = max(combinedSpositives$earliest_specimen_date,na.rm = TRUE)
  earliestDate = lastDate-28
  
  
  ## Free up some memory
  dpc$unloadCaches()
  
  ## Demographics data for Chris ----
  rawDemographics = allEpisodes %>% 
    left_join(llDemog, by="FINALID") %>%
    left_join(aocDf %>% filter(area %in% as.vector(unlist(localClusters))), by=c("LTLA_code"="code")) %>% 
    select(FINALID,eraIndex,date=earliest_specimen_date,sGene,area,LTLA_code,LTLA_name,age,sex,ethnicity_final,imd_decile,residential_category) %>%
    mutate(area = ifelse(is.na(area),"other",as.character(area))) %>%
    filter(date>=minDate)
  
  
  
  rawDemographics %>% readr::write_csv(paste0(here::here("output/raw-demographics-line-list-by-area-of-concern-"),maxDate,".csv"))
  zip::zip(
    zipfile = paste0(here::here("output/raw-demographics-line-list-by-area-of-concern-"),maxDate,".zip"),
    files = paste0(here::here("output/raw-demographics-line-list-by-area-of-concern-"),maxDate,".csv"),
    mode="cherry-pick"
  )
  unlink(paste0(here::here("output/raw-demographics-line-list-by-area-of-concern-"),maxDate,".csv"))
  
  ## Growth rate data ----

  
  # Data set for S-gene growth rates
  
  # * How do we construct the case count of S-positive and S-negative cases by area of concern & other & England?
  # * problem in that sgene line list is test level and contains >1 test per person
  # * scale of this problem
  # * sgll %>% filter(specimen_date>"2021-03-01") %>% group_by(FINALID,sgtf_under30CT) %>% summarise(count=n()) %>% group_by(sgtf_under30CT,count) %>% summarise(freq = n())
  #   
  # 1) line list left join s-gene line list by final id and specimen_date. filter on S-gene status defined (pos,neg,equivocal).
  # * this is most naive approach and will give you the subset of Pillar 2 first infections that have a s-gene status.
  # * will lose people who have a first infection that does not have a taqpath result
  # * can make this a bit flexible in terms of date matching (s-gene date within 4 days line list date...) 
  # * either way still need to deduplicate by finalid, in case there was >1 test for that finalid on that day.  in which case must also de-duplicate. 
  # * cannot tell you about re-infections.
  # * Will lose people that have long infection episode with a later test being done at TaqPath
  # * can tell you about Pillar 1 & 2 tests where taqpath results not available.
  # * FIRST_INFECTION
  # 
  # 2) s-gene line list left-join demographics from line list, - order by date and get first occurrence for each finalid. Filter by S-positive (sgtf_under30CT == 0).
  # * probably gives you similar to 1) but removes need for S-gene to be first ever test for an individual. 
  # * Instead gives first ever taqpath test for patient, which might not be the first test in an infection episode.
  # * Will ensure uniqueness - only one test per patient.
  # * Nothing for re-infections
  # * FIRST_TAQPATH
  # 
  # 3) s-gene line list left-join demographics from line list, - order by reverse date and get last occurrence for each finalid. Filter by S-positive.
  # * Last ever taqpath test for patient.
  # * Dates for patients with >1 test shifted backwards.
  # * if test represents 2 infection episodes will give you most recent episode
  # * LAST_TAQPATH
  # 
  # 4) combined ll & sgll line-list with infection episodes defined as per s-pos-data allEpisodes
  # * Date as per earliest_specimen_date
  # * Episodes are unique finalid+eraIndex but not finalid alone
  # * Potential re-infection detected.
  # * May potentially shift people much earlier than other methods if they have had a few S-gene unknown P1 & P2 tests followed by TaqPath test.
  # * INFECTION_EPISODES

  # 1) FIRST_INFECTION ----
  
  
  firstInfection = ll %>% left_join(sgll, by="FINALID", suffix=c("",".sgene")) %>% 
    filter(is.na(specimen_date.sgene) | (specimen_date <= specimen_date.sgene & specimen_date+4 > specimen_date.sgene)) %>% 
    arrange(FINALID) %>% filter(FINALID != lag(FINALID)) %>% mutate(
      sGene = case_when(
        sgtf_under30CT == 1 ~ "negative",
        sgtf_under30CT == 0 ~ "positive",
        is.na(sgtf) ~ "unknown",
        is.na(sgtf_under30CT) ~ "equivocal",
        TRUE ~ NA_character_
      )
    ) %>% select(all_of(colnames(ll)),sGene)
  
  if(any(duplicated(firstInfection$FINALID))) stop("duplicates detected")

  # 2) FIRST_TAQPATH ----
  
  firstTaqpath = sgll %>% left_join(llDemog, by="FINALID", suffix=c(".sgll","")) %>%
    left_join(llTest, by=c("FINALID","specimen_date"), suffix=c("",".test")) %>% 
    filter(!is.na(specimen_date)) %>% 
    arrange(FINALID,specimen_date) %>% 
    filter(FINALID != lag(FINALID)) %>% 
    mutate(
      pillar="Pillar 2",
      sGene = case_when(
        sgtf_under30CT == 1 ~ "negative",
        sgtf_under30CT == 0 ~ "positive",
        is.na(sgtf) ~ "unknown",
        is.na(sgtf_under30CT) ~ "equivocal",
        TRUE ~ NA_character_
      ),
      asymptomatic_indicator = ifelse(is.na(asymptomatic_indicator),"U",asymptomatic_indicator)
    )
  
  if(any(duplicated(firstTaqpath$FINALID))) stop("duplicates detected")

  # 3) LAST_TAQPATH ----
  
  lastTaqpath = sgll %>% left_join(llDemog, by="FINALID", suffix=c(".sgll","")) %>% 
    left_join(llTest, by=c("FINALID","specimen_date"), suffix=c("",".test")) %>% 
    filter(!is.na(specimen_date)) %>% 
    arrange(FINALID,desc(specimen_date)) %>% 
    filter(FINALID != lag(FINALID)) %>% 
    mutate(
      pillar="Pillar 2",
      sGene = case_when(
        sgtf_under30CT == 1 ~ "negative",
        sgtf_under30CT == 0 ~ "positive",
        is.na(sgtf) ~ "unknown",
        is.na(sgtf_under30CT) ~ "equivocal",
        TRUE ~ NA_character_
      ),
      asymptomatic_indicator = ifelse(is.na(asymptomatic_indicator),"U",asymptomatic_indicator)
    )
  
  if(any(duplicated(lastTaqpath$FINALID))) stop("duplicates detected")
  
  

  # 4) INFECTION_EPISODES
  
  infectionEpisodes = allEpisodes %>% inner_join(llDemog,by="FINALID") %>% rename(specimen_date=earliest_specimen_date)
  

  # Convert line lists into timeseries of counts ----
  
  buildTimeseries = function(filterExpr = TRUE) {
    
    filterExpr = enexpr(filterExpr)
    
    # lastTrustworthyDate = min(c(
    #   max(sgll$specimen_date,na.rm = TRUE),
    #   max(ll$specimen_date,na.rm = TRUE),
    #   max(allEpisodes$earliest_specimen_date)
    # ))-4
    
    # Since we are doing this analysis retrospectively now this is the date defined in s-pos-data.R
    lastTrustworthyDate = maxDate-1
    
    # Do combination defined sub-regional areas of England
    tmp = bind_rows(
      firstInfection %>% dpc$spim$getLineListIncidence(ageBreaks = NULL,codeTypes = "LAD", subgroup = sGene, filterExpr = !!filterExpr) %>% mutate(source="first infection"),
      firstTaqpath %>% dpc$spim$getLineListIncidence(ageBreaks = NULL,codeTypes = "LAD", subgroup = sGene, filterExpr = !!filterExpr) %>% mutate(source="first taqpath"),
      lastTaqpath %>% dpc$spim$getLineListIncidence(ageBreaks = NULL,codeTypes = "LAD", subgroup = sGene, filterExpr = !!filterExpr) %>% mutate(source="last taqpath"),
      infectionEpisodes %>% dpc$spim$getLineListIncidence(ageBreaks = NULL,codeTypes = "LAD", subgroup = sGene, filterExpr = !!filterExpr) %>% mutate(source="infection episodes")
    )
    
    tmp = tmp %>% dpc$demog$findDemographics()
    
    tmp2 = tmp %>% 
      left_join(aocDf %>% select(area,code),by = "code") %>% 
      mutate(area = forcats::fct_explicit_na(area,na_level = "other")) %>% 
      select(-code,-name) %>% 
      mutate(code=area,name=area) %>% 
      covidStandardDateGrouping() %>% 
      summarise(
        value = sum(value),
        population = sum(population,na.rm = TRUE),
        Implicit = all(Implicit)
      ) 
    
    # Do combination defined England national area
    tmp3 = bind_rows(
      firstInfection %>% dpc$spim$getLineListIncidence(ageBreaks = NULL,codeTypes = "CTRY", subgroup = sGene, filterExpr = !!filterExpr) %>% mutate(source="first infection"),
      firstTaqpath %>% dpc$spim$getLineListIncidence(ageBreaks = NULL,codeTypes = "CTRY", subgroup = sGene, filterExpr = !!filterExpr) %>% mutate(source="first taqpath"),
      lastTaqpath %>% dpc$spim$getLineListIncidence(ageBreaks = NULL,codeTypes = "CTRY", subgroup = sGene, filterExpr = !!filterExpr) %>% mutate(source="last taqpath"),
      infectionEpisodes %>% dpc$spim$getLineListIncidence(ageBreaks = NULL,codeTypes = "CTRY", subgroup = sGene, filterExpr = !!filterExpr) %>% mutate(source="infection episodes")
    )
    
    tmp3 = tmp3 %>% dpc$demog$findDemographics()
    
    tmp4 =   bind_rows(tmp2,tmp3) %>% filter(date <= lastTrustworthyDate)
    return(tmp4)
  }
  
  allTimeseries = buildTimeseries()
  symptomaticTimeseries = buildTimeseries(asymptomatic_indicator != "Y")

  allTimeseries %>% readr::write_csv(paste0(here::here("output/all-case-counts-by-area-of-concern-and-method-"),maxDate,".csv"))
  zip::zip(
    zipfile = paste0(here::here("output/all-case-counts-by-area-of-concern-and-method-"),maxDate,".zip"),
    files = paste0(here::here("output/all-case-counts-by-area-of-concern-and-method-"),maxDate,".csv"),
    mode="cherry-pick"
  )
  unlink(paste0(here::here("output/all-case-counts-by-area-of-concern-and-method-"),maxDate,".csv"))

  symptomaticTimeseries %>% readr::write_csv(paste0(here::here("output/symptomatic-case-counts-by-area-of-concern-and-method-"),maxDate,".csv"))
  zip::zip(
    zipfile = paste0(here::here("output/symptomatic-case-counts-by-area-of-concern-and-method-"),maxDate,".zip"),
    files = paste0(here::here("output/symptomatic-case-counts-by-area-of-concern-and-method-"),maxDate,".csv"),
    mode="cherry-pick"
  )
  unlink(paste0(here::here("output/symptomatic-case-counts-by-area-of-concern-and-method-"),maxDate,".csv"))
  
  
}


