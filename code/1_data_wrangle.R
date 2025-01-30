# Read in data, extract and calculate variables needed for analysis

library(dplyr)
library(lubridate)
library(summarytools)
library(mice)
library(tidylog, warn.conflicts = FALSE)
library(table1)

#read in data
diabetes_dat<-read.csv("data/Bence_dashboarddataset_wtihDMtypes_20220917.csv")
dim(diabetes_dat) #[1]904766 rows   241 columns

# view(dfSummary(diabetes_dat), file = "data/20220917_summary.html")

countries_of_origin=read.csv("data/countries.csv") %>% 
  filter(include.as.English.speaking.non.indigenous==1) #only those with 1
medsurg_codes <- read.csv("data/MedSurg_codes.csv", strip.white = TRUE)
emerg_codes <- read.csv("data/Emerg_codes.csv")
reclassify <- read.csv("data/diabetes_type_queries_for_Gordana_310522.csv") %>% 
  mutate(reclassify_diabetes = ifelse(reclassify_diabetes == 1, "Type1",
                                      ifelse(reclassify_diabetes == 2, "Type2", "No diabetes")))




# mutate first then remove
diabetes_full <- diabetes_dat %>%                                               #C1
  left_join(reclassify) %>% 
  mutate(ethnicity = ifelse(patient_country_of_origin=="Australia","aust",            #C2
                            ifelse(patient_country_of_origin %in% 
                                     countries_of_origin$patient_country_of_origin,
                                   "other_english_sp","non_english_speaking")),
         type_1_diabetes_flag = ifelse(is.na(type_1_diabetes_flag),0,type_1_diabetes_flag),
         type_2_diabetes_flag = ifelse(is.na(type_2_diabetes_flag),0,type_2_diabetes_flag),
         diabetes_in_pregnancy_flag = ifelse(is.na(diabetes_in_pregnancy_flag),0,diabetes_in_pregnancy_flag),
         insulin_first24hrs_flag = ifelse(is.na(insulin_first24hrs_flag),0,insulin_first24hrs_flag),
         dmtypes = ifelse(type_1_diabetes_flag==1, "Type1",   #overlap check    #C10  
                          ifelse(type_2_diabetes_flag==1, "Type2","No diabetes")),  
         dmtypes = ifelse(!is.na(reclassify_diabetes),reclassify_diabetes,dmtypes),
         admit_dttm=date(ymd_hms(admit_dttm)),
         ID = paste(patient_mrn,substr(osp_full_name,1,2))) %>%                 #C3
  left_join(emerg_codes) %>%   
  mutate(first_specialty_unit_description = ifelse(                             #F2
    first_specialty_unit_description %in% c("Emergency Services","Emergency"), 
    final_specialty_unit_description, first_specialty_unit_description)) %>% 
  left_join(medsurg_codes) 


table1(~osp_full_name + factor(MedSurg) + absolute_length_of_stay_days +
         referred_to_on_separation_display + factor(diabetes_in_pregnancy_flag)+
         factor(type_1_diabetes_flag) + factor(type_2_diabetes_flag), data = diabetes_full)

diabetes_full <- diabetes_full %>% 
  filter(osp_full_name %in% 
           c("Prince of Wales Hospital","St George Hospital")) %>%              #F1
  #recode first emergency to last  
  filter(MedSurg %in% 1:2) %>%   #Emergency is coded as 1 medical
  filter(absolute_length_of_stay_days>1) %>%  #exclude <= 1 day stay            #F3
  filter(absolute_length_of_stay_days<=60) %>%    #exclude >60 days stay)         #F4
  filter(!referred_to_on_separation_display %in% 
           c("Died", "Organ Procurement")) %>%                                  #F5
  filter(diabetes_in_pregnancy_flag != 1) %>% 
  filter(dmtypes %in% c("Type1", "Type2")) %>%                                #F6
  filter(admit_dttm < dmy("15/03/20")) %>%  
  filter(!(osp_full_name == "Prince of Wales Hospital" &                        #F7
             admit_dttm < dmy("01/01/17"))) %>% 
  filter(!(osp_full_name == "St George Hospital" & 
             admit_dttm < dmy("01/06/18"))) %>%
  transmute(ID,#,  ID includes hospital
            patient_mrn, 
            admit_dttm,
            osp_full_name,
            readmit = ifelse(is.na(unplanned_readmission_length_of_stay),0,1),  #C4   
            represent = ifelse(is.na(days_until_representation),0,1),           #email 14/05/2022  
            age=patient_age_at_admission,                                       #C5
            sex=ifelse(patient_sex != "Male","Not_Male",patient_sex),#no NA     #C6
            ethnicity,                                                          #C2
            dmtypes,
            lcharlson=log(charlson_comorbidity_index_excl_diabetes + 1),        #C7
            insulin=ifelse(insulin_first24hrs_flag == 1,"Insulin","No_insulin"),#C8 
            last90days=ifelse(!is.na(days_since_previous_stay),"Yes","No"),     #C9
            med_surg=MedSurg,                                                   #C11
            emergency=Emerg_Elective,                                           #C12
            MaritalStatus=ifelse(patient_marital_status=="Married or De facto", #C13
                                 "Married or De facto","Other"),
            bgl_max = ifelse(bgl_bedside_max_stay<10,"<10",
                             ifelse(bgl_bedside_max_stay <= 15,"10-15",">15")), #C14
            bgl_min=ifelse(bgl_bedside_min_stay <4,"<4",">=4"),                 #C15
            disadvantage = dis_score,                                           #C16
            length_of_stay_days) %>%                             #Meeting Sep 8 2023
  na.omit() 

# table1(~osp_full_name + factor(MedSurg) + absolute_length_of_stay_days +
#          referred_to_on_separation_display + factor(diabetes_in_pregnancy_flag)+
#          factor(type_1_diabetes_flag) + factor(type_2_diabetes_flag), data = diabetes_full)


write.csv(diabetes_full,"data/clean_alladmit.csv",row.names = FALSE)


diabetes_full <-  diabetes_full %>% 
  arrange(desc(admit_dttm)) %>% #arrange in order of date
  distinct(ID,.keep_all = TRUE)   #take last visit                             
#last admission per Peter convo
# table1(~osp_full_name + factor(MedSurg) + absolute_length_of_stay_days +
#          referred_to_on_separation_display + factor(diabetes_in_pregnancy_flag)+
#          factor(type_1_diabetes_flag) + factor(type_2_diabetes_flag), data = diabetes_full)

dim(diabetes_full) #[1] 16417    20



diabetes_full %>% 
  group_by(osp_full_name, readmit) %>% 
  count()

write.csv(diabetes_full,"data/clean_modsel.csv",row.names = FALSE)

# print summary of dataset
summarytools::view(dfSummary(diabetes_full), file="data/data_summary_clean_modsel.html")
#print summary of missing data


## separate training (POW) and test (St George) sample 

diabetes_clean_POW=diabetes_full %>% 
  filter(osp_full_name == "Prince of Wales Hospital")
diabetes_clean_StGeorge=diabetes_full %>% 
  filter(osp_full_name == "St George Hospital")


write.csv(diabetes_clean_POW,"data/clean_POW.csv", row.names = FALSE)
write.csv(diabetes_clean_StGeorge,"data/clean_StGeorge.csv", row.names = FALSE)
write.csv(diabetes_clean_StGeorge,"data/diabetes_clean_validation.csv", row.names = FALSE)





set.seed(05102022)
train_test = sample(c("train","test"), size = nrow(diabetes_clean_POW), replace = TRUE, prob = c(0.8,0.2) )
diabetes_clean_train = diabetes_clean_POW[train_test=="train",]
diabetes_clean_test = diabetes_clean_POW[train_test=="test",] 

write.csv(diabetes_clean_train,"data/diabetes_clean_train.csv", row.names = FALSE)
write.csv(diabetes_clean_test,"data/diabetes_clean_test.csv", row.names = FALSE)







