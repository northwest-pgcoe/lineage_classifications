#Links to reference:
#CDC: https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-classifications.html and https://covid.cdc.gov/covid-data-tracker/#variant-proportions
#WHO: https://www.who.int/activities/tracking-SARS-CoV-2-variants

################################
#   Load packages
################################

if (!suppressPackageStartupMessages(require(pacman))) {
  install.packages("pacman",
                   repos = "http://cran.us.r-project.org")
}

pacman::p_load(odbc,
               tidyverse,
               readxl,
               lubridate,
               fs,
               dplyr)

################################
#   Import most recent lineage file created from DIQA
################################

lineages<- read_csv("Y:/Confidential/DCHS/CDE/01_Linelists_Cross Coverage/Novel CoV/01 - Epi/Sequence Data and Reporting/Data_Objects/Lineages/Lineages.csv")

names(lineages)
length(unique(lineages$lineage_extracted))


#####################
# Check to see if there are duplicate lineages in file 
# This step done by DIQA
####################

lineages %>%
  group_by(lineage_extracted) %>%
  dplyr::summarise(count = n()) %>%
  filter(count >1)


#####################
#Filter out Withdrawn lineages
####################
active_lineages <- lineages %>%
  filter(status =="Active")


################################
#   New variables:
# Variable indicating VOC: cdc_class
# variable indicating WHO name: who_name
# variable indicating grouping in DOH S & V Report: doh_variant_name
# variable indicating hex color for doh_variant_name group: hex_code
# variable indicating reporting group of lineage: lineage_reporting_group
#                                   1 : Currently monitoring
#                                   2 : Formerly monitoring
#                                   3 : Formerly circulating, not monitored
# variable indicating variable name in numerical/pango form for tables: report_table_name
################################



#############################
# Creating new variables
############################

lineage_data_1 <- active_lineages %>%
# 'cdc_class' variable code
  mutate(vbm_class = ifelse(grepl(c("B.1.617.2|^AY.|^B\\.1\\.1\\.7$|^Q.|
                                    B.1.351|B.1.351.|^P.1|^P.1.|^B.1.427|
                                    ^B.1.429|B.1.525$|B.1.526$|B.1.617.1$|
                                    B.1.617.3$|B.1.621$|B.1.621.1$|P.2") , lineage_extracted), "VBM", "non VBM"),
        voc_class = ifelse(
                            grepl(c("B.1.1.529|XBB|XBC|XB"), lineage_extracted) |
                            grepl(c("B.1.1.529|XBB|XBC|XB") , description), 
                            "VOC", "non VOC"),
        #If adding in recombinant omicron
        cdc_class = case_when( vbm_class == "VBM" ~ "VBM",
                               voc_class == "VOC" ~ "VOC", 
                               TRUE  ~ "non VOC/VBM")) %>%
  dplyr::select(-c("vbm_class", "voc_class")) %>%
  # 'who_name' variable code
  mutate(who_name = case_when( 
      # Variants being monitored (VBM)
    #Alpha
    lineage_extracted == "B.1.1.7" | grepl("^Q.", lineage_extracted) ~ "Alpha",       #exact match to "B.1.1.7" or starts with "Q."
    #Beta
    lineage_extracted == "B.1.351" | grepl("B.1.351", lineage_extracted) ~ "Beta",   #exact match to "B.1.351" or starts with "B.1.351."
    #Gamma
    lineage_extracted == "P.1" | grepl("^P.1", lineage_extracted) ~ "Gamma",          #exact match to "P.1" or starts with "P.1."
    #Epsilon
    grepl("^B.1.427|^B.1.429", lineage_extracted) ~ "Epsilon",     
    #Eta
    lineage_extracted == "B.1.525" ~ "Eta",         #exact match to "B.1.525"
    #Iota
    lineage_extracted == "B.1.526" ~ "Iota",         #exact match to "B.1.526"
    #Kappa
    lineage_extracted == "B.1.617.1" ~ "Kappa",
    #Mu
    lineage_extracted == "B.1.621" | lineage_extracted == "B.1.621.1"~ "Mu",
    #Zeta
    lineage_extracted == "P.2" ~ "Zeta",
    #Delta
    lineage_extracted == "B.1.617.2" | grepl("^AY.", lineage_extracted) ~ "Delta",
      # Variants of concern (VOC)
    #Omicron
    grepl("B\\.1\\.1\\.529", description) | lineage_extracted =="B.1.1.529" ~ "Omicron",
    #Recombinant
    grepl("^X", lineage_extracted) ~ "Recombinant",
    grepl("XBB", description) ~ "Recombinant",
    #all non VOC/VBMs as "Other",
    TRUE ~ "Other"),
# 'doh_variant_name' variable code
    doh_variant_name = case_when( 
    #Recombinants
      #XBB.1.5
      grepl("^XBB.1.5.1", lineage_extracted) ~ "XBB.1.5.1",
      grepl("^XBB.1.5", lineage_extracted) | grepl("XBB\\.1\\.5", description) ~ "XBB.1.5",
      grepl("XBW", description) ~ "Other",
      #XBB Recombinant
      grepl("^XBB", lineage_extracted) ~ "XBB",
      grepl("^XB", lineage_extracted) ~ "Other",
      grepl("^X", lineage_extracted) ~ "Other",
    #Variants being monitored (VBM)
      #Alpha
      lineage_extracted == "B.1.1.7" | grepl("^Q.", lineage_extracted) ~ "Alpha",       #exact match to "B.1.1.7" or starts with "Q."
      #Beta
      lineage_extracted == "B.1.351" | grepl("B.1.351", lineage_extracted) ~ "Beta",   #exact match to "B.1.351" or starts with "B.1.351."
      #Gamma
      lineage_extracted == "P.1" | grepl("^P.1", lineage_extracted) ~ "Gamma",          #exact match to "P.1" or starts with "P.1."
      #Epsilon
      grepl("^B.1.427|^B.1.429", lineage_extracted) ~ "Epsilon",
      #Eta
      lineage_extracted == "B.1.525" ~ "Eta",         #exact match to "B.1.525"
      #Iota
      lineage_extracted == "B.1.526" ~ "Iota",         #exact match to "B.1.526"
      #Kappa
      lineage_extracted == "B.1.617.1" ~ "Kappa",
      #Mu
      lineage_extracted == "B.1.621" | lineage_extracted == "B.1.621.1"~ "Mu",
      #Zeta
      lineage_extracted == "P.2" ~ "Zeta",
      #Delta
      lineage_extracted == "B.1.617.2" | grepl("^AY.", lineage_extracted) ~ "Delta",
    #Variants of concern (VOC)
      #BA.1.1
      lineage_extracted == "BA.1.1" | grepl("B\\.1\\.1\\.529\\.1\\.1\\.", description) ~ "BA.1.1",
      #BA.2.12.1
      grepl("B.1.1.529.2.12.1", description) ~ "BA.2.12.1",
    #CH.1.1
      grepl("B.1.1.529.2.75.3.4.1.1.1.1", description) ~ "CH.1.1",
      #BN.1
      grepl("B.1.1.529.2.75.5.1", description) ~ "BN.1",
      #BA.2.75.2
      grepl("B.1.1.529.2.75.2", description) ~ "BA.2.75.2" ,
      #BA.2.75
      grepl("B.1.1.529.2.75", description) ~ "BA.2.75" ,
      #BA.2
      grepl("B.1.1.529.2", description) ~ "BA.2" ,
      #BA.4.6
      grepl("B.1.1.529.4.6", description) ~ "BA.4.6" ,
      #BA.4
      grepl("B.1.1.529.4", description) ~ "BA.4" ,
      #BF.7
      grepl("B.1.1.529.5.2.1.7", description) ~ "BF.7",
      #BF.11
      grepl("B.1.1.529.5.2.1.11", description) ~ "BF.11",
      #BQ.1.1
      lineage_extracted == "BQ.1.1" | grepl("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1\\.1\\.1\\.", description) ~ "BQ.1.1",
      #BQ.1
      grepl("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1\\.1", description) ~ "BQ.1",
      #BA.5.2.6
      grepl("B\\.1\\.1\\.529\\.5\\.2\\.6", description) ~ "BA.5.2.6",
      #BA.5
      grepl("B.1.1.529.5", description) | grepl("BA.5", description) ~ "BA.5" ,
      #Other Omicron
      lineage_extracted == "B.1.1.529" | grepl("B.1.1.529", description) | grepl("BA", description) ~ "Other Omicron",
      #all non VOC/VBMs as "Other",
      TRUE ~ "Other") )


###########################
# Adding in hex_code column
###########################   

lineage_data_2 <- lineage_data_1 %>%
  mutate(hex_code = case_when(doh_variant_name == "Alpha" ~ "#8dd3c7",
                  doh_variant_name == "Beta" ~ "#ffffb3",
                  doh_variant_name == "Delta" ~ "#b39ddb",
                  doh_variant_name == "Epsilon" ~ "#bebada",
                  doh_variant_name == "Eta" ~ "#fb8072",
                  doh_variant_name == "Gamma" ~ "#fdb462",
                  doh_variant_name == "Iota" ~ "#b3de69",
                  doh_variant_name == "Kappa" ~ "#fccde5",
                  doh_variant_name == "Mu" ~ "#bc80bd",
                  doh_variant_name == "Zeta" ~ "#ffed6f",
                  doh_variant_name == "Other Omicron" ~ "#e26028",
                  doh_variant_name == "BA.1.1" ~ "#ff824c",
                  doh_variant_name == "BA.2" ~ "#9ccc65",
                  doh_variant_name == "BA.2.12.1" ~ "#7cb342",
                  doh_variant_name == "BA.2.75" ~ "#d4e157",
                  doh_variant_name == "BA.2.75.2" ~ "#c0ca33",
                  doh_variant_name == "BA.4" ~ "#ffd54f",
                  doh_variant_name == "BA.4.6" ~ "#ffb300",
                  doh_variant_name == "BA.5" ~ "#80cbc4",
                  doh_variant_name == "BF.7" ~ "#81d4fa",
                  doh_variant_name == "BF.11" ~ "#29b6f6",
                  doh_variant_name == "BN.1" ~ "#9e9d24",
                  doh_variant_name == "BA.5.2.6" ~ "#009688",
                  doh_variant_name == "BQ.1" ~ "#006064",
                  doh_variant_name == "BQ.1.1" ~ "#00838f",
                  doh_variant_name == "XBB" ~ "#9fa8da",
                  doh_variant_name == "XBB.1.5" ~ "#5363bb",
                  doh_variant_name == "XBB.1.5.1" ~ "#1a237e",
                  doh_variant_name == "CH.1.1" ~ "#827717",
                  TRUE ~ "#797979"))


###########################
# Adding in 'lineage_reporting_group' column
##########################

#These two lists should be mutually exclusive
currently_monitoring_list <- c( "Other Omicron",
                          "BA.1.1",
                          "BA.2",
                          "BA.2.12.1",
                          "BA.2.75",
                          "BA.2.75.2",
                          "BA.4",
                          "BA.4.6",
                          "BA.5",
                          "BF.7",
                          "BF.11",
                          "BN.1",
                          "BA.5.2.6",
                          "BQ.1",
                          "BQ.1.1",
                          "XBB",
                          "XBB.1.5",
                          "XBB.1.5.1",
                          "CH.1.1",
                          "Other")

formerly_monitoring_list <- c("Alpha",
                       "Beta",
                       "Delta",
                       "Epsilon",
                       "Eta",
                       "Gamma",
                       "Iota",
                       "Kappa",
                       "Mu",
                       "Zeta")


lineage_data_3 <- lineage_data_2 %>%
  mutate(lineage_reporting_group = case_when(
    doh_variant_name %in% currently_monitoring_list ~ 1, 
    doh_variant_name %in% formerly_monitoring_list ~ 2, 
    TRUE ~ 3))

###########################
# Adding in report_table_name' variable 
##########################
lineage_data_final <- lineage_data_3 %>%
  mutate(doh_variant_name_tables = case_when(doh_variant_name == "Delta" ~ "B.1.617.2",
                                      doh_variant_name == "Alpha" ~ "B.1.1.7",
                                      doh_variant_name == "Beta" ~ "B.1.351",
                                      doh_variant_name == "Epsilon" ~ "B.1.427 / B.1.429",
                                      doh_variant_name == "Eta" ~ "B.1.525",
                                      doh_variant_name == "Iota" ~ "B.1.526",
                                      doh_variant_name == "Kappa" ~ "B.1.617.1",
                                      doh_variant_name == "Gamma" ~ "P.1",
                                      doh_variant_name == "Mu" ~ "B.1.621",
                                      doh_variant_name == "Zeta" ~ "P.2",
                                      doh_variant_name == "Other Omicron" ~ "B.1.1.529",
                                      TRUE ~ doh_variant_name))


###########################
# What are the new sublineages from last time
###########################   

#read in last days file
previous_lineage_data <- read_csv("Y:/Confidential/DCHS/CDE/01_Linelists_Cross Coverage/Novel CoV/01 - Epi/Sequence Data and Reporting/Data_Objects/Lineages/lineage_classifications.csv")

lineage_data_final$lineage_extracted <- as.character(lineage_data_final$lineage_extracted)
previous_lineage_data$lineage_extracted   <- as.character(previous_lineage_data$lineage_extracted)

length(previous_lineage_data)
length(lineage_data_final)

nrow(previous_lineage_data)
nrow(lineage_data_final)

#new_lineage_data <-anti_join(previous_lineage_data, lineage_data_final)

#list of ones not in previous list
new_lineage_data <- lineage_data_final %>%
  filter(!lineage_extracted %in% previous_lineage_data$lineage_extracted)

new_lineage_data



###########################
# Part 3: write out as csv
###########################   

#save(lineage_data_final, file="Y:/Confidential/DCHS/CDE/01_Linelists_Cross Coverage/Novel CoV/01 - Epi/Sequence Data and Reporting/Data_Objects/Lineages/lineage_classifications.rData")
#load("Y:/Confidential/DCHS/CDE/01_Linelists_Cross Coverage/Novel CoV/01 - Epi/Sequence Data and Reporting/Data_Objects/Lineages/lineage_classifications.rData")


write.csv(lineage_data_final, file="Y:/Confidential/DCHS/CDE/01_Linelists_Cross Coverage/Novel CoV/01 - Epi/Sequence Data and Reporting/Data_Objects/Lineages/lineage_classifications.csv", row.names=FALSE)

