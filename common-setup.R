library(data.table)
library(dtplyr)
library(tidyverse, warn.conflicts = FALSE)
library(patchwork)
library(rgdal)
library(ggplot2)
library(ggspatial)
library(rgeos)
library(maptools)
library(patchwork)
library(sp)
library(sf)

here::i_am("common-setup.R")

if (!"ukcovidtools" %in% c(devtools::dev_packages(),rownames(installed.packages()))) {
  if (fs::dir_exists("~/Git/uk-covid-datatools/")) {
    devtools::load_all("~/Git/uk-covid-datatools/")
  } else {
    devtools::install_github("terminological/uk-covid-datatools")
  }
}
ukcovidtools::setup(here::here("config.yml"))


if (!"standardPrintOutput" %in% c(devtools::dev_packages(),rownames(installed.packages()))) {
  if (fs::dir_exists("~/Git/standard-print-output/")) {
    devtools::load_all("~/Git/standard-print-output/")
  } else {
    devtools::install_github("terminological/standard-print-output")
  }
}
standardPrintOutput::setDefaults()

vocLevels = c("B.1.1.7","B.1.351","B.1.617.1","B.1.617.2","P.1 & P.2","other VOC/VUI","non VOC/VUI")
palVoc = scales::brewer_pal(palette="Set2")(7)
names(palVoc) = vocLevels
vocLevels2 = c("B.1.1.7","B.1.617.2","B.1.351","B.1.617.1","P.1 & P.2","other VOC/VUI","non VOC/VUI")

#sgenePal = c("#E41A1C","#377EB8","#808080","#B0B0B0")
sgenePal = c(
  negative="#E41A1C",
  positive="#377EB8",
  equivocal="#C020A0",
  unknown="#B0B0B0"
)