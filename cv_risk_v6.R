## Data Analysis: Daniel García-Ramírez & Omar Yaxmehen Bello-Chavolla
## Latest version of Analysis: March 2025
## For any question regarding analysis contact:
## Omar Yaxmehen Bello-Chavolla at oyaxbell@yahoo.com.mx

#### Packages ####
library(tidyverse);library(survey);library(ggpubr); library(ggfittext); library(ggsci);
library(ggsci);library(ggplot2);library(naniar); library(officer); library(flextable)
library(visdat);library(broom);library(sf); library(nephro); library(ggbreak);
library(rmapshaper);library(spdep);library(gtsummary); library(RiskScorescvd)

#setwd("~/Downloads/ENSANUT")
#setwd("~/Google Drive/Mi unidad/Proyectos/Diabetes Cardiovascular Risk/ENSANUT")
#setwd("G:/.shortcut-targets-by-id/1UwNpmg3AtIDHHkE3eHcm3kciiITS_T0h/Diabetes Cardiovascular Risk/ENSANUT_fin")

#### Custom function####
wt_prop <- function(x, design) {
  f <- as.formula(paste0("~", x))
  prop <- svymean(f, design, na.rm = TRUE)
  ci <- confint(prop)
  table <- cbind(prop, ci)
  table <- data.frame(table)
  table
}

list_extr <- function(years, variable, suffix = NULL) {
  n <- length(years)
  datalist <- vector("list", length = n)
  
  for (i in years) {
    data <- eval(parse(text = paste0("diab_", i, suffix, "[[", "\"", variable, "\"", "]]")))
    data$year <- factor(i)
    data$var <- eval(parse(text = paste0("levels(ensanut", i, "$", variable, ")")))
    datalist[[i]] <- data
  }
  
  data_final <- do.call(rbind, datalist)
  data_final <- data_final %>% rename(lower = X2.5..,
                                      upper = X97.5..)
}

wt_prop_by <- function(variable, by, design, suffix, year) {
  var <- as.formula(paste0("~", variable))
  by <- as.formula(paste0("~", by))
  table <- svyby(var, by, design, svymean, na.rm = TRUE) %>% 
    pivot_longer(-1, names_to = c(".value", "type"), names_sep = "_") %>% 
    rename(prop = suffix,
           se = paste0("se.", suffix)) %>% 
    mutate(lower = prop - (se * 1.96),
           upper = prop + (se * 1.96),
           year = year)
  table
} #Use with care

score_diab <- function(x, y) {
  library(RiskScorescvd)
  x <- as.data.frame(x)
  score <- numeric(nrow(x)) # Pre-allocate the score vector for better performance
  
  for (i in 1:nrow(x)) { # Correct loop
    score[i] <- SCORE2_Diabetes(
      Risk.region = x$Risk.region[i],
      Age = x$Age[i],
      Gender = x$Gender[i],
      smoker = x$smoker[i],
      systolic.bp = x$systolic.bp[i],
      total.chol = x$total.chol[i],
      total.hdl = x$total.hdl[i],
      diabetes = x$diabetes[i],
      diabetes.age = x$diabetes.age[i],
      HbA1c = x$HbA1c[i],
      eGFR = x$eGFR[i],
      classify = y
    )
  }
  
  return(score)
}


#### Datasets####
setwd("~/Mi unidad (obello@facmed.unam.mx)/Datasets/ENSANUT")

ensanut2016_fin <- read_csv("ensanut2016_fin.csv")
ensanut2018_fin <- read_csv("ensanut2018_fin.csv")
ensanut2021_fin <- read_csv("ensanut2021_fin.csv")
ensanut2022_fin <- read_csv("ensanut2022_fin.csv")
ensanut2023_fin <- read_csv("ensanut2023_fin.csv")

base <- read_delim("CS_RESIDENTES.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(id_edu = paste0(as.numeric(UPM), "_", as.numeric(VIV_SEL), "_", as.numeric(HOGAR), "_", as.numeric(NUMREN))) %>% 
  select(id_edu, NIVEL)

setwd("~/Mi unidad (obello@facmed.unam.mx)/CV Risk Diabetes")

#### Variables ENSANUT 2016####
ensanut2016 <- ensanut2016_fin %>% mutate(
  male = case_when(sexo.x == 1 ~ 1,
                   sexo.x == 2 ~ 0),
  imss = case_when(h211a == 1 ~ 1,
                   h211a != 1 ~ 0),
  issste = case_when(h211a == 2 |  h211a == 3 | h211b == 2 | h211b == 3 ~ 1,
                     h211a != 2 |  h211a != 3 | h211b != 2 | h211b != 3 ~ 0),
  pemex = case_when(h211a == 4 | h211b == 4 ~ 1,
                    h211a != 4 | h211b != 4 ~ 0),
  defensa = case_when(h211a == 5 | h211b == 5 ~ 1,
                      h211a != 5 | h211b != 5 ~ 0),
  seg_popular = case_when(h211a == 6 | h211b == 6 ~ 1,
                          h211a != 6 | h211b != 6 ~ 0),
  seguro_privado = case_when(h211a == 7 | h211b == 7 ~ 1,
                             h211a != 7 | h211b != 7 ~ 0),
  sin_seguro = case_when(h211a == 9 | h211b == 9 ~ 1,
                         h211a != 9 | h211b != 9 ~ 0),
  education = case_when(h218a == 0 | h218a == 1 ~ "No education",
                        h218a == 2 ~ "Elementary school",
                        h218a == 3 | h218a == 4 ~ "Middle/High school",
                        h218a == 9 | h218a == 10 | h218a == 11 | h218a == 12 ~ "University",
                        h218a == 5 |h218a == 6 | h218a == 7 | h218a == 8 ~ "Other") %>% 
    factor(levels = c("No education", "Elementary school", "Middle/High school", "University", "Other")),
  #Weight control
  wt_loss = if_else(a102 == 2, 1, 0), #1 = Perdió peso, 0 = No perdió peso
  intent_wt_loss = case_when(a104 == 1 ~ 1, #1 = Pérdida de peso intencional, 0 = Pérdida de peso no intencional
                             a104 == 2 ~ 0),
  wt_loss_qt = a103, #¿Cuántos kg fueron los que perdió?
  #Diabetes diagnosis and management
  HX_T2D_AGE = a3025,
  plan_alim = case_when(a309a == 1 ~ 1,
                        is.na(a309a) ~ 0) %>% factor(labels = c("No meal plan", "Meal plan")),
  plan_ejerc = case_when(a309b == 1 ~ 1,
                         is.na(a309b) ~ 0) %>% factor(labels = c("No excercise plan", "Excercise plan")),
  aspirin = a312b,
  smoke_quit = case_when(Smoking %in% c(0, 1) ~ "Non-smoker",
                         Smoking == 2 ~ "Smoker") %>% factor(),
  #Hypertension diagnosis and management
  HX_HBP_AGE = as.numeric(edad.x) - a402b.x,
  salt = case_when(a407e.x == 1 ~ 1,
                   is.na(a407e.x) ~ 0) %>% factor(labels = c("No salt reduction", "Salt reduction")),
  #Cholesterol management
  chol_diag = case_when(a607.x == 1 ~ 1, #1 = Diagnóstico médico de hipercolesterolemia, 0 = Sin diagnóstico médico de hipercolesterolemia
                        a607.x == 2 ~ 0),
  statin = case_when(a312i == 1 ~ 1,
                     is.na(a312i) ~ 0) %>% factor(labels = c("No statin", "Statin")),
  fibrates = case_when(a610a == 1 ~ 1,
                       is.na(a610a) ~ 0) %>% factor(labels = c("No fibrates", "Fibrates")),
  #Previous cardiovascular disease
  prev_cv = case_when((a502a == 1 | a502b == 1| a502c == 1 | a611 == 1) ~ 1,
                      (a502a != 1 | a502b != 1| a502c != 1 | a611 != 1) ~ 0),
  #s- ldl
  s_ldl = (valor.COLEST/0.948) -(valor.COL_HDL/9.971)-((valor.TRIG/8.56)+((valor.TRIG*(valor.COLEST-valor.COL_HDL))/2140)- ((valor.TRIG*valor.TRIG)/16100))-9.4

)

#Age
ensanut2016$edad.x <- if_else(ensanut2016$edad.x > 105, NA, ensanut2016$edad.x)
ensanut2016$age_group <- cut(as.numeric(ensanut2016$edad.x), breaks = c(20, 40, 60, Inf), right = FALSE, labels = c("<40 years", "40-59 years", "≥60 years"))

#Sex
ensanut2016$male1 <- factor(ensanut2016$male, labels = c("Women", "Men"))

#Diabetes
ensanut2016$diab[ensanut2016$HX_T2D == 1 | ensanut2016$HBA1C >= 6.5 | ensanut2016$TX_T2D %in% c(1, 2, 3)] <- 1

#Diagnosed and undiagnosed diabetes
ensanut2016$diag <- 0
ensanut2016$diag[ensanut2016$HX_T2D == 1 | ensanut2016$TX_T2D %in% c(1, 2, 3)] <- 1
ensanut2016$diag[ensanut2016$diag == 0 & ensanut2016$HBA1C >= 6.5] <- 2

#Diabetes treatment
ensanut2016$TX_T2D <- factor(ensanut2016$TX_T2D, labels = c("Insulin", "Pills", "Both", "None"))
ensanut2016$TX_T2D <- factor(ensanut2016$TX_T2D, levels = c("None", "Pills", "Insulin", "Both"))

#Glycemic control
ensanut2016$gly_control[(ensanut2016$edad.x < 65 & ensanut2016$HBA1C < 7) | (ensanut2016$edad.x >= 65 & ensanut2016$HBA1C < 7.5)] <- 1
ensanut2016$gly_control[(ensanut2016$edad.x < 65 & ensanut2016$HBA1C >= 7) | (ensanut2016$edad.x >= 65 & ensanut2016$HBA1C >= 7.5)] <- 0
ensanut2016$gly_control <- factor(ensanut2016$gly_control, labels = c("Uncontrolled", "Controlled"))

#Hypertension treatment (in individuals with diagnosed hypertension)
ensanut2016$TX_HBP <- factor(ensanut2016$TX_HBP, labels = c("Untreated", "Treated"))

#Uncontrolled hypertension
ensanut2016$high_bp[ensanut2016$SBP >= 130 | ensanut2016$DBP >= 80] <- 1
ensanut2016$high_bp[ensanut2016$SBP < 130 & ensanut2016$DBP < 80] <- 0
ensanut2016$high_bp <- factor(ensanut2016$high_bp, labels = c("Controlled", "Uncontrolled"))

#Uncontrolled LDL
ensanut2016$high_ldl2[ensanut2016$LDL >= 70] <- 1
ensanut2016$high_ldl2[ensanut2016$LDL < 70] <- 0
ensanut2016$high_ldl2 <- factor(ensanut2016$high_ldl2, labels = c("<70 mg/dL", "≥70 mg/dL"))

#Uncontrolled S-LDL
ensanut2016$high_sldl[ensanut2016$s_ldl >= 70] <- 1
ensanut2016$high_sldl[ensanut2016$s_ldl < 70] <- 0
ensanut2016$high_sldl <- factor(ensanut2016$high_sldl, labels = c("<70 mg/dL", "≥70 mg/dL"))

#Uncontrolled S-LDL (>100 mg/dL)
ensanut2016$high_sldl100[ensanut2016$s_ldl >= 100] <- 1
ensanut2016$high_sldl100[ensanut2016$s_ldl < 100] <- 0
ensanut2016$high_sldl100 <- factor(ensanut2016$high_sldl100, labels = c("<100 mg/dL", "≥100 mg/dL"))

#Uncontrolled HDL
ensanut2016$high_hdl[(ensanut2016$HDL <40 & ensanut2016$male == 1) | (ensanut2016$HDL <50 & ensanut2016$male == 0)] <- 1
ensanut2016$high_hdl[(ensanut2016$HDL >=40 & ensanut2016$male == 1) | (ensanut2016$HDL >=50 & ensanut2016$male == 0)] <- 0
ensanut2016$high_hdl <- factor(ensanut2016$high_hdl, labels = c("Controlled", "Uncontrolled"))

#Smoking status
ensanut2016$Smoking <- factor(ensanut2016$Smoking, labels = c("Never smoker", "Former smoker", "Current smoker"))

#Glycemic, blood pressure and cholesterol control (LDL-C <70 mg/dL)
ensanut2016$comb_abc[ensanut2016$gly_control == "Controlled" & ensanut2016$high_bp == "Controlled" & ensanut2016$high_sldl == "<70 mg/dL"] <- 1
ensanut2016$comb_abc[ensanut2016$gly_control == "Uncontrolled" | ensanut2016$high_bp == "Uncontrolled" | ensanut2016$high_sldl == "≥70 mg/dL"] <- 0
ensanut2016$comb_abc <- factor(ensanut2016$comb_abc, labels = c("Uncontrolled", "Control"))

#Blood pressure and cholesterol control (LDL-C <70 mg/dL)
ensanut2016$comb_bc[ensanut2016$high_bp == "Controlled" & ensanut2016$high_sldl == "<70 mg/dL"] <- 1
ensanut2016$comb_bc[ensanut2016$high_bp == "Uncontrolled" | ensanut2016$high_sldl == "≥70 mg/dL"] <- 0
ensanut2016$comb_bc <- factor(ensanut2016$comb_bc, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure, cholesterol control and no smoker (LDL-C <70 mg/dL)
ensanut2016$comb_abcn[ensanut2016$gly_control == "Uncontrolled" | ensanut2016$high_bp == "Uncontrolled" | ensanut2016$high_sldl == "≥70 mg/dL" | ensanut2016$smoke_quit == "Smoker"] <- 0
ensanut2016$comb_abcn[ensanut2016$gly_control == "Controlled" & ensanut2016$high_bp == "Controlled" & ensanut2016$high_sldl == "<70 mg/dL" & ensanut2016$smoke_quit == "Non-smoker"] <- 1
ensanut2016$comb_abcn <- factor(ensanut2016$comb_abcn, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure and cholesterol control (LDL-C <100 mg/dL)
ensanut2016$comb_abc100[ensanut2016$gly_control == "Controlled" & ensanut2016$high_bp == "Controlled" & ensanut2016$high_sldl100 == "<100 mg/dL"] <- 1
ensanut2016$comb_abc100[ensanut2016$gly_control == "Uncontrolled" | ensanut2016$high_bp == "Uncontrolled" | ensanut2016$high_sldl100 == "≥100 mg/dL"] <- 0
ensanut2016$comb_abc100 <- factor(ensanut2016$comb_abc100, labels = c("Uncontrolled", "Control"))

#Blood pressure and cholesterol control (LDL-C <100 mg/dL)
ensanut2016$comb_bc100[ensanut2016$high_bp == "Controlled" & ensanut2016$high_sldl100 == "<100 mg/dL"] <- 1
ensanut2016$comb_bc100[ensanut2016$high_bp == "Uncontrolled" | ensanut2016$high_sldl100 == "≥100 mg/dL"] <- 0
ensanut2016$comb_bc100 <- factor(ensanut2016$comb_bc100, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure, cholesterol control and no smoker (LDL-C <100 mg/dL)
ensanut2016$comb_abcn100[ensanut2016$gly_control == "Uncontrolled" | ensanut2016$high_bp == "Uncontrolled" | ensanut2016$high_sldl100 == "≥100 mg/dL" | ensanut2016$smoke_quit == "Smoker"] <- 0
ensanut2016$comb_abcn100[ensanut2016$gly_control == "Controlled" & ensanut2016$high_bp == "Controlled" & ensanut2016$high_sldl100 == "<100 mg/dL" & ensanut2016$smoke_quit == "Non-smoker"] <- 1
ensanut2016$comb_abcn100 <- factor(ensanut2016$comb_abcn100, labels = c("Uncontrolled", "Control"))

#Type of settlement (rural or urban)
ensanut2016$urban <- factor(ensanut2016$Area, labels = c("Rural", "Urban"))

#Indigenous
ensanut2016$indigenous_fct <- factor(ensanut2016$Indigenous, labels = c("Non-indigenous", "Indigenous"))

#Social security
ensanut2016$social_sec[ensanut2016$imss == 1 | ensanut2016$issste == 1 | ensanut2016$pemex == 1 | ensanut2016$defensa == 1 | ensanut2016$seg_popular == 1 | ensanut2016$seguro_privado == 1] <- 1
ensanut2016$social_sec[ensanut2016$sin_seguro == 1] <- 0
ensanut2016$social_sec <- factor(ensanut2016$social_sec, labels = c("Without social security", "With social security"))

#Diagnostic time
ensanut2016$edad_diag <- if_else(ensanut2016$a3025 == 99, NA, ensanut2016$a3025)
ensanut2016$diag_time <- as.numeric(ensanut2016$edad.x) - as.numeric(ensanut2016$edad_diag)
ensanut2016$diag_time <- if_else(ensanut2016$diag_time < 0, NA, ensanut2016$diag_time)
ensanut2016$cronic_diag[ensanut2016$diag_time >= 10] <- 1
ensanut2016$cronic_diag[ensanut2016$diag_time < 10 ] <- 0
ensanut2016$cronic_diag <- factor(ensanut2016$cronic_diag, labels = c("<10 years of diagnosis", "≥10 years of diagnosis"))

#High blood pressure diagnosis
ensanut2016$bp_diag[ensanut2016$HX_HBP == 0 & (ensanut2016$SBP < 140 & ensanut2016$DBP < 90)] <- 0 #Sin hipertensión
ensanut2016$bp_diag[ensanut2016$HX_HBP == 1] <- 1 #Diagnosticado
ensanut2016$bp_diag[ensanut2016$HX_HBP == 0 & (ensanut2016$SBP >= 140 | ensanut2016$DBP >= 90)] <- 2 #No diagnosticado
ensanut2016$bp_diag <- factor(ensanut2016$bp_diag, labels = c("Without hypertension", "Diagnosed hypertension", "Undiagnosed hypertension"))

#Primary and secondary prevention
ensanut2016$prevent[ensanut2016$edad.x <40] <- 0
ensanut2016$prevent[ensanut2016$edad.x >= 40 & ensanut2016$prev_cv == 0] <- 1
ensanut2016$prevent[ensanut2016$edad.x >= 40 & ensanut2016$prev_cv == 1] <- 2
ensanut2016$prevent <- factor(ensanut2016$prevent, labels = c("<40 years", "Elegible for primary prevention", "Elegible for secondary prevention"))

#Hypertension treatment (regardless of hypertension diagnosis)
ensanut2016$tx_bp_all <- if_else((ensanut2016$bp_diag == "Without hypertension" & is.na(ensanut2016$TX_HBP)), "Untreated", ensanut2016$TX_HBP)
ensanut2016$tx_bp_all <- if_else((ensanut2016$bp_diag == "Undiagnosed hypertension" & is.na(ensanut2016$TX_HBP)), "Untreated", ensanut2016$tx_bp_all)
ensanut2016$tx_bp_all <- factor(ensanut2016$tx_bp_all)

#Diabetes complications
ensanut2016$diabetes_retinopatia<-if_else((ensanut2016$a313c==1 | ensanut2016$a313d==1 | ensanut2016$a313e==1)==T,"Retinopathy", "No retinopathy")
ensanut2016$diabetes_nefropatia<-if_else(ensanut2016$a313f==1,"Nephropathy", "No nephropathy")
ensanut2016$diabetes_neuropatia<-if_else(ensanut2016$a313i==1,"Neuropathy", "No neuropathy")
ensanut2016$diabetes_anycomplication<-ifelse((ensanut2016$diabetes_retinopatia=="Retinopathy"|
                                                ensanut2016$diabetes_nefropatia=="Nephropathy"|
                                                ensanut2016$diabetes_neuropatia=="Neuropathy"), "Yes", "No")

ensanut2016$diabetes_allcomplications<-ifelse((ensanut2016$diabetes_retinopatia=="Retinopathy" &
                                                ensanut2016$diabetes_nefropatia=="Nephropathy" &
                                                ensanut2016$diabetes_neuropatia=="Neuropathy"), "Yes", "No")

#Cardiovascular outcomes and risk
ensanut2016$ascvd<-ifelse((ensanut2016$a502a==1 | ensanut2016$a502b==1)==T, "ASCVD", "No ASCVD")
ensanut2016$ethnicity <- rep(0, nrow(ensanut2016))
ensanut2016$eGFR<- with(ensanut2016, CKDEpi.creat(
  creatinine=Creatinine, sex=Sex, age=Age, ethnicity=ethnicity))
ensanut2016$severe_ckd<- ifelse(ensanut2016$eGFR<30, "eGFR <30mL/min/1.73m2", "eGFR ≥30mL/min/1.73m2")
ensanut2016$Risk.region<-"Moderate"
ensanut2016$total.chol<-ensanut2016$valor.COLEST/38.67
ensanut2016$total.hdl<-ensanut2016$valor.COL_HDL/38.67
ensanut2016$HbA1c<-10.929*(ensanut2016$valor.HB1AC-2.15)
ensanut2016$systolic.bp<-ensanut2016$SBP
ensanut2016$diabetes.age<-ensanut2016$edad_diag
ensanut2016$diabetes<-ensanut2016$diag
ensanut2016$Gender<-as.character(factor(ensanut2016$male1, labels = c("female", "male")))
ensanut2016$smoker<-ifelse(ensanut2016$Smoking=="Current smoker",1,0)
ensanut2016$classify<-TRUE
ensanut2016$SCORE2_diabetes<-score_diab(ensanut2016, FALSE)
ensanut2016$SCORE2_diabetes_cat<-score_diab(ensanut2016, TRUE)
ensanut2016$SCORE2_diabetes_cat2<-ensanut2016$SCORE2_diabetes_cat
ensanut2016$SCORE2_diabetes_cat2[ensanut2016$ascvd=="ASCVD"]<-"Very high risk"
ensanut2016$SCORE2_diabetes_cat2[ensanut2016$diabetes_allcomplications=="Yes"]<-"Very high risk"
ensanut2016$high_risk<-ifelse(ensanut2016$SCORE2_diabetes_cat2 %in% c("High risk", "Very high risk"), 1, 0)

## LDL guidelines
ensanut2016$high_ldl<-ifelse(((ensanut2016$SCORE2_diabetes_cat2=="Very high risk" & ensanut2016$s_ldl<55)|
                                      (ensanut2016$SCORE2_diabetes_cat2=="High risk" & ensanut2016$s_ldl<70) |
                                      (ensanut2016$SCORE2_diabetes_cat2 %in% c("Low risk", "Moderate risk") & ensanut2016$s_ldl<100)), 1, 0)
ensanut2016$high_ldl<-factor(ensanut2016$high_ldl, labels = c("Uncontrolled", "coNTROLLED"))

#### Variables ENSANUT 2018####
ensanut2018 <- ensanut2018_fin %>% mutate(
  EDAD.x = as.numeric(EDAD.x),
  male = case_when(SEXO.x == 1 ~ 1,
                   SEXO.x == 2 ~ 0),
  imss = case_when(P3_10_01 == 1 ~ 1,  
                   P3_10_01 != 1 ~ 0),
  imss_prospera = case_when(P3_10_08 == 1 ~ 1,  
                            P3_10_08 != 1 ~ 0),
  issste = case_when(P3_10_02 == 1 |  P3_10_03 == 1  ~ 1,
                     P3_10_02 != 1 |  P3_10_03 != 1 ~ 0),
  pemex = case_when(P3_10_04 == 1 ~ 1,  
                    P3_10_04 != 1 ~ 0),
  defensa = case_when(P3_10_05 == 1 ~ 1,  
                      P3_10_05 != 1 ~ 0),
  marina = case_when(P3_10_06 == 1 ~ 1,  
                     P3_10_06 != 1 ~ 0),
  seg_popular = case_when(P3_10_07 == 1 ~ 1,  
                          P3_10_07 != 1 ~ 0),
  seguro_privado = case_when(P3_10_09 == 1 ~ 1,  
                             P3_10_09 != 1 ~ 0),
  otro = case_when(P3_10_10 == 1 ~ 1,  
                   P3_10_10 != 1 ~ 0),
  sin_seguro = case_when(P3_10_11 == 1 ~ 1,  
                         P3_10_11 != 1 ~ 0),
  #Weight control
  wt_loss = if_else(P1_7 == 2, 1, 0),
  intent_wt_loss = case_when(P1_9 == 1 ~ 1, #1 = Sí, 0 = No
                             P1_9 == 2 ~ 0),
  wt_loss_qt = P1_8, #¿Cuántos kg fueron los que perdió?
  #Diabetes diagnosis and management
  HX_T2D_AGE = P3_2,
  plan_alim = case_when(P3_13_1 == 1 ~ 1, 
                        P3_13_1 != 1 ~ 0) %>% factor(labels = c("No meal plan", "Meal plan")), 
  plan_ejerc = case_when(P3_13_2 == 1  ~ 1,
                         P3_13_1 != 1 ~ 0) %>% factor(labels = c("No excercise plan", "Excercise plan")),
  aspirin = case_when(P3_16_1 == 2 | P3_16_2 == 2 ~ 1, 
                      P3_16_1 != 2 | P3_16_2 != 2  ~ 0),
  smoke_quit = case_when(Smoking %in% c(0, 1) ~ "Non-smoker",
                         Smoking == 2 ~ "Smoker") %>% factor(),
  #Hypertension diagnosis and management
  HX_HBP_AGE = as.numeric(EDAD.x) - as.numeric(P4_2A),
  salt = case_when(P4_8_3 == 1 ~ 1,
                   P4_8_3 != 1  ~ 0) %>% factor(labels = c("No salt reduction", "Salt reduction")),
  #Cholesterol management
  chol_diag = if_else(P6_4 == 1, 1, 0), #1 = Diagnóstico médico de hipercolesterolemia, 0 = Sin diagnóstico médico de hipercolesterolemia
  statin = case_when(P6_8_3 == 1 ~ 1,
                     P6_8_3 != 1  ~ 0) %>% factor(labels = c("No statin", "Statin")),
  fibrates = case_when(P6_8_4 == 1 ~ 1, 
                      P6_8_4 != 1 ~ 0) %>% factor(labels = c("No fibrates", "Fibrates")),
  #Previous cardiovascular disease 
  prev_cv = case_when(P5_2_1 == 1 | P5_2_2 == 1 | P5_2_3 == 1  ~ 1, #1 = Con CVD previa, 0 = Sin CVD
                    P5_2_1 == 2 | P5_2_2 == 2 | P5_2_3 == 2  ~ 0),
  #s_ldl
  s_ldl = (VALOR_COLEST/0.948) -(VALOR_COL_HDL/9.971)-((VALOR_TRIG/8.56)+((VALOR_TRIG*(VALOR_COLEST-VALOR_COL_HDL))/2140)- ((VALOR_TRIG*VALOR_TRIG)/16100))-9.4
)

#Age
ensanut2018$age_group <- cut(as.numeric(ensanut2018$EDAD.x), breaks = c(20, 40, 60, Inf), right = FALSE, labels = c("<40 years", "40-59 years", "≥60 years"))

#Sex
ensanut2018$male1 <- factor(ensanut2018$male, labels = c("Women", "Men"))

#Diabetes
ensanut2018$diab[ensanut2018$HX_T2D == 1 | ensanut2018$HBA1C >= 6.5 | ensanut2018$TX_T2D %in% c(1, 2, 3)] <- 1

#Diagnosed and undiagnosed diabetes
ensanut2018$diag <- 0 
ensanut2018$diag[ensanut2018$HX_T2D == 1 | ensanut2018$TX_T2D %in% c(1, 2, 3)] <- 1
ensanut2018$diag[ensanut2018$diag == 0 & ensanut2018$HBA1C >= 6.5] <- 2

#Diabetes treatment
ensanut2018$TX_T2D <- factor(ensanut2018$TX_T2D, labels = c("Insulin", "Pills", "Both", "None"))
ensanut2018$TX_T2D <- factor(ensanut2018$TX_T2D, levels = c("None", "Pills", "Insulin", "Both"))

#Glycemic control
ensanut2018$gly_control[(ensanut2018$EDAD.x < 65 & ensanut2018$HBA1C < 7) | (ensanut2018$EDAD.x >= 65 & ensanut2018$HBA1C < 7.5)] <- 1
ensanut2018$gly_control[(ensanut2018$EDAD.x < 65 & ensanut2018$HBA1C >= 7) | (ensanut2018$EDAD.x >= 65 & ensanut2018$HBA1C >= 7.5)] <- 0
ensanut2018$gly_control <- factor(ensanut2018$gly_control, labels = c("Uncontrolled", "Controlled"))

#Hypertension treatment
ensanut2018$TX_HBP <- factor(ensanut2018$TX_HBP, labels = c("Untreated", "Treated"))

#Uncontrolled hypertension
ensanut2018$high_bp[ensanut2018$SBP >= 130 | ensanut2018$DBP >= 80] <- 1
ensanut2018$high_bp[ensanut2018$SBP < 130 & ensanut2018$DBP < 80] <- 0
ensanut2018$high_bp <- factor(ensanut2018$high_bp, labels = c("Controlled", "Uncontrolled"))

#Uncontrolled LDL
ensanut2018$high_ldl2[ensanut2018$LDL >= 70] <- 1
ensanut2018$high_ldl2[ensanut2018$LDL < 70] <- 0
ensanut2018$high_ldl2 <- factor(ensanut2018$high_ldl2, labels = c("<70 mg/dL", "≥70 mg/dL"))

#Uncontrolled S-LDL
ensanut2018$high_sldl[ensanut2018$s_ldl >= 70] <- 1
ensanut2018$high_sldl[ensanut2018$s_ldl < 70] <- 0
ensanut2018$high_sldl <- factor(ensanut2018$high_sldl, labels = c("<70 mg/dL", "≥70 mg/dL"))

#Uncontrolled S-LDL (>100 mg/dL)
ensanut2018$high_sldl100[ensanut2018$s_ldl >= 100] <- 1
ensanut2018$high_sldl100[ensanut2018$s_ldl < 100] <- 0
ensanut2018$high_sldl100 <- factor(ensanut2018$high_sldl100, labels = c("<100 mg/dL", "≥100 mg/dL"))

#Uncontrolled HDL
ensanut2018$high_hdl[(ensanut2018$HDL <40 & ensanut2018$male == 1) | (ensanut2018$HDL <50 & ensanut2018$male == 0)] <- 1
ensanut2018$high_hdl[(ensanut2018$HDL >=40 & ensanut2018$male == 1) | (ensanut2018$HDL >=50 & ensanut2018$male == 0)] <- 0
ensanut2018$high_hdl <- factor(ensanut2018$high_hdl, labels = c("Controlled", "Uncontrolled"))

#Smoking status
ensanut2018$Smoking <- factor(ensanut2018$Smoking, labels = c("Never smoker", "Former smoker", "Current smoker"))

#Glycemic, blood pressure and cholesterol control
ensanut2018$comb_abc[ensanut2018$gly_control == "Controlled" & ensanut2018$high_bp == "Controlled" & ensanut2018$high_sldl == "<70 mg/dL"] <- 1
ensanut2018$comb_abc[ensanut2018$gly_control == "Uncontrolled" | ensanut2018$high_bp == "Uncontrolled" | ensanut2018$high_sldl == "≥70 mg/dL"] <- 0
ensanut2018$comb_abc <- factor(ensanut2018$comb_abc, labels = c("Uncontrolled", "Control"))

#Blood pressure and cholesterol control
ensanut2018$comb_bc[ensanut2018$high_bp == "Controlled" & ensanut2018$high_sldl == "<70 mg/dL"] <- 1
ensanut2018$comb_bc[ensanut2018$high_bp == "Uncontrolled" | ensanut2018$high_sldl == "≥70 mg/dL"] <- 0
ensanut2018$comb_bc <- factor(ensanut2018$comb_bc, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure, cholesterol control and no smoker
ensanut2018$comb_abcn[ensanut2018$gly_control == "Uncontrolled" | ensanut2018$high_bp == "Uncontrolled" | ensanut2018$high_sldl == "≥70 mg/dL" | ensanut2018$smoke_quit == "Smoker"] <- 0
ensanut2018$comb_abcn[ensanut2018$gly_control == "Controlled" & ensanut2018$high_bp == "Controlled" & ensanut2018$high_sldl == "<70 mg/dL" & ensanut2018$smoke_quit == "Non-smoker"] <- 1
ensanut2018$comb_abcn <- factor(ensanut2018$comb_abcn, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure and cholesterol control (LDL-C <100 mg/dL)
ensanut2018$comb_abc100[ensanut2018$gly_control == "Controlled" & ensanut2018$high_bp == "Controlled" & ensanut2018$high_sldl100 == "<100 mg/dL"] <- 1
ensanut2018$comb_abc100[ensanut2018$gly_control == "Uncontrolled" | ensanut2018$high_bp == "Uncontrolled" | ensanut2018$high_sldl100 == "≥100 mg/dL"] <- 0
ensanut2018$comb_abc100 <- factor(ensanut2018$comb_abc100, labels = c("Uncontrolled", "Control"))

#Blood pressure and cholesterol control (LDL-C <100 mg/dL)
ensanut2018$comb_bc100[ensanut2018$high_bp == "Controlled" & ensanut2018$high_sldl100 == "<100 mg/dL"] <- 1
ensanut2018$comb_bc100[ensanut2018$high_bp == "Uncontrolled" | ensanut2018$high_sldl100 == "≥100 mg/dL"] <- 0
ensanut2018$comb_bc100 <- factor(ensanut2018$comb_bc100, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure, cholesterol control and no smoker (LDL-C <100 mg/dL)
ensanut2018$comb_abcn100[ensanut2018$gly_control == "Uncontrolled" | ensanut2018$high_bp == "Uncontrolled" | ensanut2018$high_sldl100 == "≥100 mg/dL" | ensanut2018$smoke_quit == "Smoker"] <- 0
ensanut2018$comb_abcn100[ensanut2018$gly_control == "Controlled" & ensanut2018$high_bp == "Controlled" & ensanut2018$high_sldl100 == "<100 mg/dL" & ensanut2018$smoke_quit == "Non-smoker"] <- 1
ensanut2018$comb_abcn100 <- factor(ensanut2018$comb_abcn100, labels = c("Uncontrolled", "Control"))

#Type of settlement (rural or urban)
ensanut2018$urban <- factor(ensanut2018$Area, labels = c("Rural", "Urban"))

#Indigenous
ensanut2018$indigenous_fct <- factor(ensanut2018$Indigenous, labels = c("Non-indigenous", "Indigenous"))

#Education
ensanut2018 <- ensanut2018 %>%
  mutate(id_edu = paste0(as.numeric(UPM.x), "_", as.numeric(VIV_SEL.x), "_", as.numeric(HOGAR.x), "_", as.numeric(NUMREN.x)))

ensanut2018 <- left_join(ensanut2018, base, by = "id_edu") %>% 
  mutate(NIVEL = as.numeric(NIVEL.y),
         education = case_when(NIVEL == 0 | NIVEL == 1 ~ "No education",
                               NIVEL == 2 ~ "Elementary school",
                               NIVEL == 3 | NIVEL == 4 ~ "Middle/High school",
                               NIVEL == 9 | NIVEL == 10 | NIVEL == 11 | NIVEL == 12 ~ "University",
                               NIVEL == 5 | NIVEL == 6 | NIVEL == 7 | NIVEL == 8 ~ "Other") %>%
           factor(levels = c("No education", "Elementary school", "Middle/High school", "University", "Other")))

#Social security
ensanut2018$social_sec[ensanut2018$imss == 1 | ensanut2018$imss_prospera == 1 | ensanut2018$issste == 1 | ensanut2018$pemex == 1 | ensanut2018$defensa == 1 | ensanut2018$marina == 1 | ensanut2018$seg_popular == 1 | ensanut2018$seguro_privado == 1] <- 1
ensanut2018$social_sec[ensanut2018$sin_seguro == 1] <- 0
ensanut2018$social_sec <- factor(ensanut2018$social_sec, labels = c("Without social security", "With social security"))

#Diagnostic time
ensanut2018$edad_diag <- if_else(ensanut2018$P3_2 == 99, NA, ensanut2018$P3_2)
ensanut2018$diag_time <- as.numeric(ensanut2018$EDAD.x) - as.numeric(ensanut2018$edad_diag)
ensanut2018$diag_time <- if_else(ensanut2018$diag_time < 0, NA, ensanut2018$diag_time)
ensanut2018$cronic_diag[ensanut2018$diag_time >= 10] <- 1
ensanut2018$cronic_diag[ensanut2018$diag_time < 10 ] <- 0
ensanut2018$cronic_diag <- factor(ensanut2018$cronic_diag, labels = c("<10 years of diagnosis", "≥10 years of diagnosis"))

#High blood pressure diagnosis
ensanut2018$bp_diag[ensanut2018$HX_HBP == 0 & (ensanut2018$SBP < 140 & ensanut2018$DBP < 90)] <- 0 #Sin hipertensión
ensanut2018$bp_diag[ensanut2018$HX_HBP == 1] <- 1 #Diagnosticado
ensanut2018$bp_diag[ensanut2018$HX_HBP == 0 & (ensanut2018$SBP >= 140 | ensanut2018$DBP >= 90)] <- 2 #No diagnosticado
ensanut2018$bp_diag <- factor(ensanut2018$bp_diag, labels = c("Without hypertension", "Diagnosed hypertension", "Undiagnosed hypertension"))

#Primary and secondary prevention
ensanut2018$prevent[ensanut2018$EDAD.x <40] <- 0
ensanut2018$prevent[ensanut2018$EDAD.x >= 40 & ensanut2018$prev_cv == 0] <- 1
ensanut2018$prevent[ensanut2018$EDAD.x >= 40 & ensanut2018$prev_cv == 1] <- 2
ensanut2018$prevent <- factor(ensanut2018$prevent, labels = c("<40 years", "Elegible for primary prevention", "Elegible for secondary prevention"))

#Hypertension treatment (regardless of hypertension diagnosis)
ensanut2018$tx_bp_all <- if_else((ensanut2018$bp_diag == "Without hypertension" & is.na(ensanut2018$TX_HBP)), "Untreated", ensanut2018$TX_HBP)
ensanut2018$tx_bp_all <- if_else((ensanut2018$bp_diag == "Undiagnosed hypertension" & is.na(ensanut2018$TX_HBP)), "Untreated", ensanut2018$tx_bp_all)
ensanut2018$tx_bp_all <- factor(ensanut2018$tx_bp_all)

#Geographical region
ensanut2018$ENT.x <- as.numeric(ensanut2018$ENT.x)
ensanut2018$SUBREGION[ensanut2018$ENT.x %in% c(2,3,18,25,26)]<-1
ensanut2018$SUBREGION[ensanut2018$ENT.x %in% c(5,8,19,28)]<-2
ensanut2018$SUBREGION[ensanut2018$ENT.x %in% c(6,14,16)]<-3
ensanut2018$SUBREGION[ensanut2018$ENT.x %in% c(1,10,11,22,24,32)]<-4
ensanut2018$SUBREGION[ensanut2018$ENT.x %in% c(13,29,30)]<-5
ensanut2018$SUBREGION[ensanut2018$ENT.x %in% c(9)]<-6
ensanut2018$SUBREGION[ensanut2018$ENT.x %in% c(15)]<-7
ensanut2018$SUBREGION[ensanut2018$ENT.x %in% c(12,17,20,21)]<-8
ensanut2018$SUBREGION[ensanut2018$ENT.x %in% c(4,7,27,23,31)]<-9

#Diabetes complications
ensanut2018$diabetes_retinopatia<-if_else((ensanut2018$P3_18_3 ==1 | ensanut2018$P3_18_4==1)==T,"Retinopathy", "No retinopathy")
ensanut2018$diabetes_nefropatia<-if_else(ensanut2018$P3_18_5==1,"Nephropathy", "No nephropathy")
ensanut2018$diabetes_neuropatia<-if_else(ensanut2018$P3_18_1==1,"Neuropathy", "No neuropathy")
ensanut2018$diabetes_anycomplication<-ifelse((ensanut2018$diabetes_retinopatia=="Retinopathy"|
                                                ensanut2018$diabetes_nefropatia=="Nephropathy"|
                                                ensanut2018$diabetes_neuropatia=="Neuropathy"), "Yes", "No")

ensanut2018$diabetes_allcomplications<-ifelse((ensanut2018$diabetes_retinopatia=="Retinopathy" &
                                                 ensanut2018$diabetes_nefropatia=="Nephropathy" &
                                                 ensanut2018$diabetes_neuropatia=="Neuropathy"), "Yes", "No")

#Cardiovascular outcomes and risk
ensanut2018$ascvd<-ifelse((ensanut2018$P5_2_1==1 | ensanut2018$P5_2_2==1)==T, "ASCVD", "No ASCVD")
ensanut2018$ethnicity <- rep(0, nrow(ensanut2018))
ensanut2018$eGFR<- with(ensanut2018, CKDEpi.creat(
  creatinine=Creatinine, sex=Sex, age=Age, ethnicity=ethnicity))
ensanut2018$severe_ckd<- ifelse(ensanut2018$eGFR<30, "eGFR <30mL/min/1.73m2", "eGFR ≥30mL/min/1.73m2")
ensanut2018$Risk.region<-"Moderate"
ensanut2018$total.chol<-ensanut2018$VALOR_COLEST/38.67
ensanut2018$total.hdl<-ensanut2018$VALOR_COL_HDL/38.67
ensanut2018$HbA1c<-10.929*(ensanut2018$VALOR_HB1AC-2.15)
ensanut2018$systolic.bp<-ensanut2018$SBP
ensanut2018$diabetes.age<-as.numeric(ensanut2018$edad_diag)
ensanut2018$diabetes<-ensanut2018$diag
ensanut2018$Gender<-as.character(factor(ensanut2018$male1, labels = c("female", "male")))
ensanut2018$smoker<-ifelse(ensanut2018$Smoking=="Current smoker",1,0)
ensanut2018$classify<-TRUE
ensanut2018$SCORE2_diabetes<-score_diab(ensanut2018, FALSE)
ensanut2018$SCORE2_diabetes_cat<-score_diab(ensanut2018, TRUE)
ensanut2018$SCORE2_diabetes_cat2<-ensanut2018$SCORE2_diabetes_cat
ensanut2018$SCORE2_diabetes_cat2[ensanut2018$ascvd=="ASCVD"]<-"Very high risk"
ensanut2018$SCORE2_diabetes_cat2[ensanut2018$diabetes_allcomplications=="Yes"]<-"Very high risk"
ensanut2018$high_risk<-ifelse(ensanut2018$SCORE2_diabetes_cat2 %in% c("High risk", "Very high risk"), 1, 0)
## LDL guidelines
ensanut2018$high_ldl<-ifelse(((ensanut2018$SCORE2_diabetes_cat2=="Very high risk" & ensanut2018$s_ldl<55)|
                                      (ensanut2018$SCORE2_diabetes_cat2=="High risk" & ensanut2018$s_ldl<70) |
                                      (ensanut2018$SCORE2_diabetes_cat2 %in% c("Low risk", "Moderate risk") & ensanut2018$s_ldl<100)), 1, 0)
ensanut2018$high_ldl<-factor(ensanut2018$high_ldl, labels = c("Uncontrolled", "coNTROLLED"))

#### Variables ENSANUT 2021####
ensanut2021 <- ensanut2021_fin %>% mutate(
  male = case_when(sexo == 1 ~ 1,
                   sexo == 2 ~ 0),
  imss = case_when(H0310A == 1 ~ 1,
                   H0310A != 1 ~ 0),
  issste = case_when(H0310A == 2 |  H0310A == 3 | H0310B == 2 | H0310B == 3 ~ 1,
                     H0310A != 2 |  H0310A != 3 | H0310B != 2 | H0310B != 3 ~ 0),
  pemex = case_when(H0310A == 4 | H0310B == 4 ~ 1,
                    H0310A != 4 | H0310B != 4 ~ 0),
  defensa = case_when(H0310A == 5 | H0310B == 5 ~ 1,
                      H0310A != 5 | H0310B != 5 ~ 0),
  marina = case_when(H0310A == 6 | H0310B == 6 ~ 1,
                     H0310A != 6 | H0310B != 6 ~ 0),
  imss_bienestar = case_when(H0310A == 7 | H0310B == 7 ~ 1,
                             H0310A != 7 | H0310B != 7 ~ 0),
  seguro_privado = case_when(H0310A == 8 | H0310B == 8 ~ 1,
                             H0310A != 8 | H0310B != 8 ~ 0),
  sin_seguro = case_when(H0310A == 10 | H0310B == 10 ~ 1,
                         H0310A != 10 | H0310B != 10 ~ 0),
  education = case_when(h0317a == 0 | h0317a == 1 ~ "No education",
                        h0317a == 2 ~ "Elementary school",
                        h0317a == 3 | h0317a == 4 ~ "Middle/High school",
                        h0317a == 9 | h0317a == 10 | h0317a == 11 | h0317a == 12 ~ "University",
                        h0317a == 5 |h0317a == 6 | h0317a == 7 | h0317a == 8 ~ "Other") %>% 
    factor(levels = c("No education", "Elementary school", "Middle/High school", "University", "Other")),
  #Weight control
  wt_loss = if_else(a0107 == 2, 1, 0), #1 = Perdió peso, 0 = No perdió peso
  intent_wt_loss = case_when(a0109 == 1 ~ 1, #1 = Sí, 0 = No
                             a0109 == 2 ~ 0),
  wt_loss_qt = a0108,
  #Diabetes diagnosis and management
  HX_T2D_AGE = a0302,
  plan_alim = case_when(A0312A == 1 ~ 1, 
                        A0312A %in% c(2:5) ~ 0) %>% factor(labels = c("No meal plan", "Meal plan")),
  plan_ejerc = case_when(A0312A == 2 | A0312B == 2 | A0312C == 2 | A0312D == 2 ~ 1,
                         A0312A != 2 | A0312B != 2 | A0312C != 2 | A0312D != 2 ~ 0) %>% factor(labels = c("No excercise plan", "Excercise plan")),
  aspirin = case_when(A0315A == 2 | A0315B == 2 ~ 1, 
                      A0315A != 2 | A0315B != 2  ~ 0),
  smoke_quit = case_when(Smoking %in% c(0, 1) ~ "Non-smoker",
                         Smoking == 2 ~ "Smoker") %>% factor(),
  #Hypertension diagnosis and management
  HX_HBP_AGE = edad - a0402a,
  salt = case_when(A0408A == 3 | A0408B == 3 | A0408C == 3 ~ 1, 
                   A0408A != 3 | A0408B != 3 | A0408C != 3  ~ 0) %>% factor(labels = c("No salt reduction", "Salt reduction")),
  #Cholesterol management
  chol_diag = if_else(a0604 == 1, 1, 0), #1 = Diagnóstico médico de hipercolesterolemia, 0 = Sin diagnóstico médico de hipercolesterolemia
  statin = case_when(A0608A == 1 ~ 1, 
                     A0608A != 1 ~ 0) %>% factor(labels = c("No statin", "Statin")),
  fibrates = case_when(A0608A == 2 | A0608B == 2 ~ 1,
                       A0608A != 2 | A0608B != 2 ~ 0) %>% factor(labels = c("No fibrates", "Fibrates")),
  #Previous cardiovascular disease
  prev_cv = case_when(a0502a == 1 | a0502b == 1 | a0502c == 1 | a0506 == 1 ~ 1, #1 = Con CVD previa, 0 = Sin CVD
                      a0502a == 2 | a0502b == 2 | a0502c == 2 | a0506 == 2 ~ 0),
  #s_ldl
  s_ldl = (valor_COLEST/0.948)-(valor_COL_HDL/9.971)-((valor_TRIG/8.56)+((valor_TRIG*(valor_COLEST-valor_COL_HDL))/2140)-((valor_TRIG*valor_TRIG)/16100))-9.4
)

#Age
ensanut2021$age_group <- cut(as.numeric(ensanut2021$edad), breaks = c(20, 40, 60, Inf), right = FALSE, labels = c("<40 years", "40-59 years", "≥60 years"))

#Sex
ensanut2021$male1 <- factor(ensanut2021$male, labels = c("Women", "Men"))

#Diabetes
ensanut2021$diab[ensanut2021$HX_T2D == 1 | ensanut2021$HBA1C >= 6.5 | ensanut2021$TX_T2D %in% c(1, 2, 3)] <- 1

#Diagnosed and undiagnosed diabetes
ensanut2021$diag <- 0
ensanut2021$diag[ensanut2021$HX_T2D == 1 | ensanut2021$TX_T2D %in% c(1, 2, 3)] <- 1
ensanut2021$diag[ensanut2021$diag == 0 & ensanut2021$HBA1C >= 6.5] <- 2

#Diabetes treatment
ensanut2021$TX_T2D <- factor(ensanut2021$TX_T2D, labels = c("Insulin", "Pills", "Both", "None"))
ensanut2021$TX_T2D <- factor(ensanut2021$TX_T2D, levels = c("None", "Pills", "Insulin", "Both"))

#Glycemic control
ensanut2021$gly_control[(ensanut2021$edad < 65 & ensanut2021$HBA1C < 7) | (ensanut2021$edad >= 65 & ensanut2021$HBA1C < 7.5)] <- 1
ensanut2021$gly_control[(ensanut2021$edad < 65 & ensanut2021$HBA1C >= 7) | (ensanut2021$edad >= 65 & ensanut2021$HBA1C >= 7.5)] <- 0
ensanut2021$gly_control <- factor(ensanut2021$gly_control, labels = c("Uncontrolled", "Controlled"))

#Hypertension treatment
ensanut2021$TX_HBP <- factor(ensanut2021$TX_HBP, labels = c("Untreated", "Treated"))

#Uncontrolled hypertension
ensanut2021$high_bp[ensanut2021$SBP >= 130 | ensanut2021$DBP >= 80] <- 1
ensanut2021$high_bp[ensanut2021$SBP < 130 & ensanut2021$DBP < 80] <- 0
ensanut2021$high_bp <- factor(ensanut2021$high_bp, labels = c("Controlled", "Uncontrolled"))

#Uncontrolled LDL
ensanut2021$high_ldl2[ensanut2021$LDL >= 70] <- 1
ensanut2021$high_ldl2[ensanut2021$LDL < 70] <- 0
ensanut2021$high_ldl2 <- factor(ensanut2021$high_ldl2, labels = c("<70 mg/dL", "≥70 mg/dL"))

#Uncontrolled S-LDL
ensanut2021$high_sldl[ensanut2021$s_ldl >= 70] <- 1
ensanut2021$high_sldl[ensanut2021$s_ldl < 70] <- 0
ensanut2021$high_sldl <- factor(ensanut2021$high_sldl, labels = c("<70 mg/dL", "≥70 mg/dL"))

#Uncontrolled S-LDL (>100 mg/dL)
ensanut2021$high_sldl100[ensanut2021$s_ldl >= 100] <- 1
ensanut2021$high_sldl100[ensanut2021$s_ldl < 100] <- 0
ensanut2021$high_sldl100 <- factor(ensanut2021$high_sldl100, labels = c("<100 mg/dL", "≥100 mg/dL"))

#Uncontrolled HDL
ensanut2021$high_hdl[(ensanut2021$HDL <40 & ensanut2021$male == 1) | (ensanut2021$HDL <50 & ensanut2021$male == 0)] <- 1
ensanut2021$high_hdl[(ensanut2021$HDL >=40 & ensanut2021$male == 1) | (ensanut2021$HDL >=50 & ensanut2021$male == 0)] <- 0
ensanut2021$high_hdl <- factor(ensanut2021$high_hdl, labels = c("Controlled", "Uncontrolled"))

#Smoking status
ensanut2021$Smoking <- factor(ensanut2021$Smoking, labels = c("Never smoker", "Former smoker", "Current smoker"))

#Glycemic, blood pressure and cholesterol control
ensanut2021$comb_abc[ensanut2021$gly_control == "Controlled" & ensanut2021$high_bp == "Controlled" & ensanut2021$high_sldl == "<70 mg/dL"] <- 1
ensanut2021$comb_abc[ensanut2021$gly_control == "Uncontrolled" | ensanut2021$high_bp == "Uncontrolled" | ensanut2021$high_sldl == "≥70 mg/dL"] <- 0
ensanut2021$comb_abc <- factor(ensanut2021$comb_abc, labels = c("Uncontrolled", "Control"))

#Blood pressure and cholesterol control
ensanut2021$comb_bc[ensanut2021$high_bp == "Controlled" & ensanut2021$high_sldl == "<70 mg/dL"] <- 1
ensanut2021$comb_bc[ensanut2021$high_bp == "Uncontrolled" | ensanut2021$high_sldl == "≥70 mg/dL"] <- 0
ensanut2021$comb_bc <- factor(ensanut2021$comb_bc, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure, cholesterol control and no smoker
ensanut2021$comb_abcn[ensanut2021$gly_control == "Uncontrolled" | ensanut2021$high_bp == "Uncontrolled" | ensanut2021$high_sldl == "≥70 mg/dL" | ensanut2021$smoke_quit == "Smoker"] <- 0
ensanut2021$comb_abcn[ensanut2021$gly_control == "Controlled" & ensanut2021$high_bp == "Controlled" & ensanut2021$high_sldl == "<70 mg/dL" & ensanut2021$smoke_quit == "Non-smoker"] <- 1
ensanut2021$comb_abcn <- factor(ensanut2021$comb_abcn, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure and cholesterol control (LDL-C <100 mg/dL)
ensanut2021$comb_abc100[ensanut2021$gly_control == "Controlled" & ensanut2021$high_bp == "Controlled" & ensanut2021$high_sldl100 == "<100 mg/dL"] <- 1
ensanut2021$comb_abc100[ensanut2021$gly_control == "Uncontrolled" | ensanut2021$high_bp == "Uncontrolled" | ensanut2021$high_sldl100 == "≥100 mg/dL"] <- 0
ensanut2021$comb_abc100 <- factor(ensanut2021$comb_abc100, labels = c("Uncontrolled", "Control"))

#Blood pressure and cholesterol control (LDL-C <100 mg/dL)
ensanut2021$comb_bc100[ensanut2021$high_bp == "Controlled" & ensanut2021$high_sldl100 == "<100 mg/dL"] <- 1
ensanut2021$comb_bc100[ensanut2021$high_bp == "Uncontrolled" | ensanut2021$high_sldl100 == "≥100 mg/dL"] <- 0
ensanut2021$comb_bc100 <- factor(ensanut2021$comb_bc100, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure, cholesterol control and no smoker (LDL-C <100 mg/dL)
ensanut2021$comb_abcn100[ensanut2021$gly_control == "Uncontrolled" | ensanut2021$high_bp == "Uncontrolled" | ensanut2021$high_sldl100 == "≥100 mg/dL" | ensanut2021$smoke_quit == "Smoker"] <- 0
ensanut2021$comb_abcn100[ensanut2021$gly_control == "Controlled" & ensanut2021$high_bp == "Controlled" & ensanut2021$high_sldl100 == "<100 mg/dL" & ensanut2021$smoke_quit == "Non-smoker"] <- 1
ensanut2021$comb_abcn100 <- factor(ensanut2021$comb_abcn100, labels = c("Uncontrolled", "Control"))

#Type of settlement (rural or urban)
ensanut2021$urban <- factor(ensanut2021$Area, labels = c("Rural", "Urban"))

#Indigenous
ensanut2021$indigenous_fct <- factor(ensanut2021$Indigenous, labels = c("Non-indigenous", "Indigenous"))

#Social security
ensanut2021$social_sec[ensanut2021$imss == 1 | ensanut2021$issste == 1 | ensanut2021$pemex == 1 | ensanut2021$defensa == 1 | ensanut2021$marina == 1 | ensanut2021$imss_bienestar == 1 | ensanut2021$seguro_privado == 1] <- 1
ensanut2021$social_sec[ensanut2021$sin_seguro == 1] <- 0
ensanut2021$social_sec <- factor(ensanut2021$social_sec, labels = c("Without social security", "With social security"))

#Geographic distribution
ensanut2021$entidad.x <- as.numeric(ensanut2021$entidad.x)
ensanut2021$SUBREGION[ensanut2021$entidad.x %in% c(2,3,18,25,26)]<-1
ensanut2021$SUBREGION[ensanut2021$entidad.x %in% c(5,8,19,28)]<-2
ensanut2021$SUBREGION[ensanut2021$entidad.x %in% c(6,14,16)]<-3
ensanut2021$SUBREGION[ensanut2021$entidad.x %in% c(1,10,11,22,24,32)]<-4
ensanut2021$SUBREGION[ensanut2021$entidad.x %in% c(13,29,30)]<-5
ensanut2021$SUBREGION[ensanut2021$entidad.x %in% c(9)]<-6
ensanut2021$SUBREGION[ensanut2021$entidad.x %in% c(15)]<-7
ensanut2021$SUBREGION[ensanut2021$entidad.x %in% c(12,17,20,21)]<-8
ensanut2021$SUBREGION[ensanut2021$entidad.x %in% c(4,7,27,23,31)]<-9

#Diagnostic time
ensanut2021$edad_diag <- if_else(ensanut2021$a0302 == 99, NA, ensanut2021$a0302)
ensanut2021$diag_time <- as.numeric(ensanut2021$edad) - as.numeric(ensanut2021$edad_diag)
ensanut2021$diag_time <- if_else(ensanut2021$diag_time < 0, NA, ensanut2021$diag_time)
ensanut2021$cronic_diag[ensanut2021$diag_time >= 10] <- 1
ensanut2021$cronic_diag[ensanut2021$diag_time < 10 ] <- 0
ensanut2021$cronic_diag <- factor(ensanut2021$cronic_diag, labels = c("<10 years of diagnosis", "≥10 years of diagnosis"))

#High blood pressure diagnosis
ensanut2021$bp_diag[ensanut2021$HX_HBP == 0 & (ensanut2021$SBP < 140 & ensanut2021$DBP < 90)] <- 0 #Sin hipertensión
ensanut2021$bp_diag[ensanut2021$HX_HBP == 1] <- 1 #Diagnosticado
ensanut2021$bp_diag[ensanut2021$HX_HBP == 0 & (ensanut2021$SBP >= 140 | ensanut2021$DBP >= 90)] <- 2 #No diagnosticado
ensanut2021$bp_diag <- factor(ensanut2021$bp_diag, labels = c("Without hypertension", "Diagnosed hypertension", "Undiagnosed hypertension"))

#Primary and secondary prevention
ensanut2021$prevent[ensanut2021$edad <40] <- 0
ensanut2021$prevent[ensanut2021$edad >= 40 & ensanut2021$prev_cv == 0] <- 1
ensanut2021$prevent[ensanut2021$edad >= 40 & ensanut2021$prev_cv == 1] <- 2
ensanut2021$prevent <- factor(ensanut2021$prevent, labels = c("<40 years", "Elegible for primary prevention", "Elegible for secondary prevention"))

#Hypertension treatment (regardless of hypertension diagnosis)
ensanut2021$tx_bp_all <- if_else((ensanut2021$bp_diag == "Without hypertension" & is.na(ensanut2021$TX_HBP)), "Untreated", ensanut2021$TX_HBP)
ensanut2021$tx_bp_all <- if_else((ensanut2021$bp_diag == "Undiagnosed hypertension" & is.na(ensanut2021$TX_HBP)), "Untreated", ensanut2021$tx_bp_all)
ensanut2021$tx_bp_all <- factor(ensanut2021$tx_bp_all)

#Diabetes complications
ensanut2021$diabetes_retinopatia<-if_else((ensanut2021$a0316c ==1 | ensanut2021$a0316d==1)==T,"Retinopathy", "No retinopathy")
ensanut2021$diabetes_nefropatia<-if_else(ensanut2021$a0316e==1,"Nephropathy", "No nephropathy")
ensanut2021$diabetes_neuropatia<-if_else(ensanut2021$a0316a==1,"Neuropathy", "No neuropathy")
ensanut2021$diabetes_anycomplication<-ifelse((ensanut2021$diabetes_retinopatia=="Retinopathy"|
                                                ensanut2021$diabetes_nefropatia=="Nephropathy"|
                                                ensanut2021$diabetes_neuropatia=="Neuropathy"), "Yes", "No")

ensanut2021$diabetes_allcomplications<-ifelse((ensanut2021$diabetes_retinopatia=="Retinopathy" &
                                                 ensanut2021$diabetes_nefropatia=="Nephropathy" &
                                                 ensanut2021$diabetes_neuropatia=="Neuropathy"), "Yes", "No")

#Cardiovascular outcomes and risk
ensanut2021$ascvd<-ifelse((ensanut2021$a0502a==1 | ensanut2021$a0502b==1)==T, "ASCVD", "No ASCVD")
ensanut2021$ethnicity <- rep(0, nrow(ensanut2021))
ensanut2021$eGFR<- with(ensanut2021, CKDEpi.creat(
  creatinine=Creatinine, sex=Sex, age=Age, ethnicity=ethnicity))
ensanut2021$severe_ckd<- ifelse(ensanut2021$eGFR<30, "eGFR <30mL/min/1.73m2", "eGFR ≥30mL/min/1.73m2")
ensanut2021$Risk.region<-"Moderate"
ensanut2021$total.chol<-ensanut2021$valor_COLEST/38.67
ensanut2021$total.hdl<-ensanut2021$valor_COL_HDL/38.67
ensanut2021$HbA1c<-10.929*(ensanut2021$valor_HB1AC-2.15)
ensanut2021$systolic.bp<-ensanut2021$SBP
ensanut2021$diabetes.age<-as.numeric(ensanut2021$edad_diag)
ensanut2021$diabetes<-ensanut2021$diag
ensanut2021$Gender<-as.character(factor(ensanut2021$male1, labels = c("female", "male")))
ensanut2021$smoker<-ifelse(ensanut2021$Smoking=="Current smoker",1,0)
ensanut2021$classify<-TRUE
ensanut2021$SCORE2_diabetes<-score_diab(ensanut2021, FALSE)
ensanut2021$SCORE2_diabetes_cat<-score_diab(ensanut2021, TRUE)
ensanut2021$SCORE2_diabetes_cat2<-ensanut2021$SCORE2_diabetes_cat
ensanut2021$SCORE2_diabetes_cat2[ensanut2021$ascvd=="ASCVD"]<-"Very high risk"
ensanut2021$SCORE2_diabetes_cat2[ensanut2021$diabetes_allcomplications=="Yes"]<-"Very high risk"
ensanut2021$high_risk<-ifelse(ensanut2021$SCORE2_diabetes_cat2 %in% c("High risk", "Very high risk"), 1, 0)

## LDL guidelines
ensanut2021$high_ldl<-ifelse(((ensanut2021$SCORE2_diabetes_cat2=="Very high risk" & ensanut2021$s_ldl<55)|
                                      (ensanut2021$SCORE2_diabetes_cat2=="High risk" & ensanut2021$s_ldl<70) |
                                      (ensanut2021$SCORE2_diabetes_cat2 %in% c("Low risk", "Moderate risk") & ensanut2021$s_ldl<100)), 1, 0)
ensanut2021$high_ldl<-factor(ensanut2021$high_ldl, labels = c("Uncontrolled", "coNTROLLED"))

#### Variables ENSANUT 2022####
ensanut2022 <- ensanut2022_fin %>% mutate(
  male = case_when(h0302.x == 1 ~ 1,
                   h0302.x == 2 ~ 0),
  H0310A = as.numeric(H0310A),
  H0310B = as.numeric(H0310B),
  imss = case_when(H0310A == 1 ~ 1,
                   H0310A != 1 ~ 0),
  issste = case_when(H0310A == 2 |  H0310A == 3 | H0310B == 2 | H0310B == 3 ~ 1,
                     H0310A != 2 |  H0310A != 3 | H0310B != 2 | H0310B != 3 ~ 0),
  pemex = case_when(H0310A == 4 | H0310B == 4 ~ 1,
                    H0310A != 4 | H0310B != 4 ~ 0),
  defensa = case_when(H0310A == 5 | H0310B == 5 ~ 1,
                      H0310A != 5 | H0310B != 5 ~ 0),
  marina = case_when(H0310A == 6 | H0310B == 6 ~ 1,
                     H0310A != 6 | H0310B != 6 ~ 0),
  imss_bienestar = case_when(H0310A == 7 | H0310B == 7 ~ 1,
                             H0310A != 7 | H0310B != 7 ~ 0),
  seguro_privado = case_when(H0310A == 8 | H0310B == 8 ~ 1,
                             H0310A != 8 | H0310B != 8 ~ 0),
  sin_seguro = case_when(H0310A == 10 | H0310B == 10 ~ 1,
                         H0310A != 10 | H0310B != 10 ~ 0),
  education = case_when(h0317a == 0 | h0317a == 1 ~ "No education",
                        h0317a == 2 ~ "Elementary school",
                        h0317a == 3 | h0317a == 4 ~ "Middle/High school",
                        h0317a == 9 | h0317a == 10 | h0317a == 11 | h0317a == 12 ~ "University",
                        h0317a == 5 |h0317a == 6 | h0317a == 7 | h0317a == 8 ~ "Other") %>% 
    factor(levels = c("No education", "Elementary school", "Middle/High school", "University", "Other")),
  #Weight control
  wt_loss = if_else(a0107 == 2, 1, 0), #1 = Perdió peso, 0 = No perdió peso
  intent_wt_loss = case_when(a0109 == 1 ~ 1, #1 = Sí, 0 = No
                             a0109 == 2 ~ 0),
  wt_loss_qt = a0108,
  #Diabetes diagnosis and management
  HX_T2D_AGE = a0302,
  plan_alim = case_when(A0312A == 1 ~ 1, 
                        A0312A %in% c(2:5) ~ 0) %>% factor(labels = c("No meal plan", "Meal plan")),
  plan_ejerc = case_when(A0312A == 2 | A0312B == 2 | A0312C == 2 | A0312D == 2 ~ 1,
                         A0312A != 2 | A0312B != 2 | A0312C != 2 | A0312D != 2 ~ 0) %>% factor(labels = c("No excercise plan", "Excercise plan")),
  aspirin = case_when(A0315A == 2 | A0315B == 2 ~ 1, 
                      A0315A != 2 | A0315B != 2  ~ 0),
  smoke_quit = case_when(Smoking %in% c(0, 1) ~ "Non-smoker",
                         Smoking == 2 ~ "Smoker") %>% factor(),
  #Hypertension diagnosis and management
  HX_HBP_AGE = edad - a0402a,
  salt = case_when(A0408A == 3 | A0408B == 3 | A0408C == 3 ~ 1, 
                   A0408A != 3 | A0408B != 3 | A0408C != 3  ~ 0) %>% factor(labels = c("No salt reduction", "Salt reduction")),
  #Cholesterol management
  chol_diag = if_else(a0604 == 1, 1, 0), #1 = Diagnóstico médico de hipercolesterolemia, 0 = Sin diagnóstico médico de hipercolesterolemia
  statin = case_when(A0605A == 1 ~ 1, 
                     A0605A != 1 ~ 0) %>% factor(labels = c("No statin", "Statin")),
  fibrates = case_when(A0607A1 == 2  ~ 1,
                       A0607A1 != 2 ~ 0) %>% factor(labels = c("No fibrates", "Fibrates")),
  #Previous cardiovascular disease
  prev_cv = case_when(a0502a == 1 | a0502b == 1 | a0502c == 1  ~ 1, #1 = Con CVD previa, 0 = Sin CVD
                      a0502a == 2 | a0502b == 2  | a0502c == 2 ~ 0), 
 #s_ldl
  s_ldl= (valor_COLEST/0.948) -(valor_COL_HDL/9.971)-((valor_TRIG/8.56)+((valor_TRIG*(valor_COLEST-valor_COL_HDL))/2140)- ((valor_TRIG*valor_TRIG)/16100))-9.4 
)

#Age
ensanut2022$age_group <- cut(as.numeric(ensanut2022$edad), breaks = c(20, 40, 60, Inf), right = FALSE, labels = c("<40 years", "40-59 years", "≥60 years"))

#Sex
ensanut2022$male1 <- factor(ensanut2022$male, labels = c("Women", "Men"))

#Diabetes
ensanut2022$diab[ensanut2022$HX_T2D == 1 | ensanut2022$HBA1C >= 6.5 | ensanut2022$TX_T2D %in% c(1, 2, 3)] <- 1

#Diagnosed and undiagnosed diabetes
ensanut2022$diag <- 0
ensanut2022$diag[ensanut2022$HX_T2D == 1 | ensanut2022$TX_T2D %in% c(1, 2, 3)] <- 1
ensanut2022$diag[ensanut2022$diag == 0 & ensanut2022$HBA1C >= 6.5] <- 2

#Diabetes treatment
ensanut2022$TX_T2D <- factor(ensanut2022$TX_T2D, labels = c("Insulin", "Pills", "Both", "None"))
ensanut2022$TX_T2D <- factor(ensanut2022$TX_T2D, levels = c("None", "Pills", "Insulin", "Both"))

#Glycemic control
ensanut2022$gly_control[(ensanut2022$h0303.x < 65 & ensanut2022$HBA1C < 7) | (ensanut2022$h0303.x >= 65 & ensanut2022$HBA1C < 7.5)] <- 1
ensanut2022$gly_control[(ensanut2022$h0303.x < 65 & ensanut2022$HBA1C >= 7) | (ensanut2022$h0303.x >= 65 & ensanut2022$HBA1C >= 7.5)] <- 0
ensanut2022$gly_control <- factor(ensanut2022$gly_control, labels = c("Uncontrolled", "Controlled"))

#Hypertension treatment
ensanut2022$TX_HBP <- factor(ensanut2022$TX_HBP, labels = c("Untreated", "Treated"))

#Uncontrolled hypertension
ensanut2022$high_bp[ensanut2022$SBP >= 130 | ensanut2022$DBP >= 80] <- 1
ensanut2022$high_bp[ensanut2022$SBP < 130 & ensanut2022$DBP < 80] <- 0
ensanut2022$high_bp <- factor(ensanut2022$high_bp, labels = c("Controlled", "Uncontrolled"))

#Uncontrolled LDL
ensanut2022$high_ldl2[ensanut2022$LDL >= 70] <- 1
ensanut2022$high_ldl2[ensanut2022$LDL < 70] <- 0
ensanut2022$high_ldl2 <- factor(ensanut2022$high_ldl2, labels = c("<70 mg/dL", "≥70 mg/dL"))

#Uncontrolled S-LDL
ensanut2022$high_sldl[ensanut2022$s_ldl >= 70] <- 1
ensanut2022$high_sldl[ensanut2022$s_ldl < 70] <- 0
ensanut2022$high_sldl <- factor(ensanut2022$high_sldl, labels = c("<70 mg/dL", "≥70 mg/dL"))

#Uncontrolled S-LDL (>100 mg/dL)
ensanut2022$high_sldl100[ensanut2022$s_ldl >= 100] <- 1
ensanut2022$high_sldl100[ensanut2022$s_ldl < 100] <- 0
ensanut2022$high_sldl100 <- factor(ensanut2022$high_sldl100, labels = c("<100 mg/dL", "≥100 mg/dL"))

#Uncontrolled HDL
ensanut2022$high_hdl[(ensanut2022$HDL <40 & ensanut2022$male == 1) | (ensanut2022$HDL <50 & ensanut2022$male == 0)] <- 1
ensanut2022$high_hdl[(ensanut2022$HDL >=40 & ensanut2022$male == 1) | (ensanut2022$HDL >=50 & ensanut2022$male == 0)] <- 0
ensanut2022$high_hdl <- factor(ensanut2022$high_hdl, labels = c("Controlled", "Uncontrolled"))

#Smoking status
ensanut2022$Smoking <- factor(ensanut2022$Smoking, labels = c("Never smoker", "Former smoker", "Current smoker"))

#Glycemic, blood pressure and cholesterol control
ensanut2022$comb_abc[ensanut2022$gly_control == "Controlled" & ensanut2022$high_bp == "Controlled" & ensanut2022$high_sldl == "<70 mg/dL"] <- 1
ensanut2022$comb_abc[ensanut2022$gly_control == "Uncontrolled" | ensanut2022$high_bp == "Uncontrolled" | ensanut2022$high_sldl == "≥70 mg/dL"] <- 0
ensanut2022$comb_abc <- factor(ensanut2022$comb_abc, labels = c("Uncontrolled", "Control"))

#Blood pressure and cholesterol control
ensanut2022$comb_bc[ensanut2022$high_bp == "Controlled" & ensanut2022$high_sldl == "<70 mg/dL"] <- 1
ensanut2022$comb_bc[ensanut2022$high_bp == "Uncontrolled" | ensanut2022$high_sldl == "≥70 mg/dL"] <- 0
ensanut2022$comb_bc <- factor(ensanut2022$comb_bc, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure, cholesterol control and no smoker
ensanut2022$comb_abcn[ensanut2022$gly_control == "Uncontrolled" | ensanut2022$high_bp == "Uncontrolled" | ensanut2022$high_sldl == "≥70 mg/dL" | ensanut2022$smoke_quit == "Smoker"] <- 0
ensanut2022$comb_abcn[ensanut2022$gly_control == "Controlled" & ensanut2022$high_bp == "Controlled" & ensanut2022$high_sldl == "<70 mg/dL" & ensanut2022$smoke_quit == "Non-smoker"] <- 1
ensanut2022$comb_abcn <- factor(ensanut2022$comb_abcn, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure and cholesterol control (LDL-C <100 mg/dL)
ensanut2022$comb_abc100[ensanut2022$gly_control == "Controlled" & ensanut2022$high_bp == "Controlled" & ensanut2022$high_sldl100 == "<100 mg/dL"] <- 1
ensanut2022$comb_abc100[ensanut2022$gly_control == "Uncontrolled" | ensanut2022$high_bp == "Uncontrolled" | ensanut2022$high_sldl100 == "≥100 mg/dL"] <- 0
ensanut2022$comb_abc100 <- factor(ensanut2022$comb_abc100, labels = c("Uncontrolled", "Control"))

#Blood pressure and cholesterol control (LDL-C <100 mg/dL)
ensanut2022$comb_bc100[ensanut2022$high_bp == "Controlled" & ensanut2022$high_sldl100 == "<100 mg/dL"] <- 1
ensanut2022$comb_bc100[ensanut2022$high_bp == "Uncontrolled" | ensanut2022$high_sldl100 == "≥100 mg/dL"] <- 0
ensanut2022$comb_bc100 <- factor(ensanut2022$comb_bc100, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure, cholesterol control and no smoker (LDL-C <100 mg/dL)
ensanut2022$comb_abcn100[ensanut2022$gly_control == "Uncontrolled" | ensanut2022$high_bp == "Uncontrolled" | ensanut2022$high_sldl100 == "≥100 mg/dL" | ensanut2022$smoke_quit == "Smoker"] <- 0
ensanut2022$comb_abcn100[ensanut2022$gly_control == "Controlled" & ensanut2022$high_bp == "Controlled" & ensanut2022$high_sldl100 == "<100 mg/dL" & ensanut2022$smoke_quit == "Non-smoker"] <- 1
ensanut2022$comb_abcn100 <- factor(ensanut2022$comb_abcn100, labels = c("Uncontrolled", "Control"))

#Type of settlement (rural or urban)
ensanut2022$urban <- factor(ensanut2022$Area, labels = c("Rural", "Urban"))

#Indigenous
ensanut2022$indigenous_fct <- factor(ensanut2022$Indigenous, labels = c("Non-indigenous", "Indigenous"))

#Social security
ensanut2022$social_sec[ensanut2022$imss == 1 | ensanut2022$issste == 1 | ensanut2022$pemex == 1 | ensanut2022$defensa == 1 | ensanut2022$marina == 1 | ensanut2022$imss_bienestar == 1 | ensanut2022$seguro_privado == 1] <- 1
ensanut2022$social_sec[ensanut2022$sin_seguro == 1] <- 0
ensanut2022$social_sec <- factor(ensanut2022$social_sec, labels = c("Without social security", "With social security"))

#Geographic regions
ensanut2022$entidad.x <- as.numeric(ensanut2022$entidad.x)
ensanut2022$SUBREGION[ensanut2022$entidad.x %in% c(2,3,18,25,26)]<-1
ensanut2022$SUBREGION[ensanut2022$entidad.x %in% c(5,8,19,28)]<-2
ensanut2022$SUBREGION[ensanut2022$entidad.x %in% c(6,14,16)]<-3
ensanut2022$SUBREGION[ensanut2022$entidad.x %in% c(1,10,11,22,24,32)]<-4
ensanut2022$SUBREGION[ensanut2022$entidad.x %in% c(13,29,30)]<-5
ensanut2022$SUBREGION[ensanut2022$entidad.x %in% c(9)]<-6
ensanut2022$SUBREGION[ensanut2022$entidad.x %in% c(15)]<-7
ensanut2022$SUBREGION[ensanut2022$entidad.x %in% c(12,17,20,21)]<-8
ensanut2022$SUBREGION[ensanut2022$entidad.x %in% c(4,7,27,23,31)]<-9

#Diagnostic time
ensanut2022$edad_diag <- if_else(ensanut2022$a0302 == 99, NA, ensanut2022$a0302)
ensanut2022$diag_time <- as.numeric(ensanut2022$edad) - as.numeric(ensanut2022$edad_diag)
ensanut2022$diag_time <- if_else(ensanut2022$diag_time < 0, NA, ensanut2022$diag_time)
ensanut2022$cronic_diag[ensanut2022$diag_time >= 10] <- 1
ensanut2022$cronic_diag[ensanut2022$diag_time < 10 ] <- 0
ensanut2022$cronic_diag <- factor(ensanut2022$cronic_diag, labels = c("<10 years of diagnosis", "≥10 years of diagnosis"))

#High blood pressure diagnosis
ensanut2022$bp_diag[ensanut2022$HX_HBP == 0 & (ensanut2022$SBP < 140 & ensanut2022$DBP < 90)] <- 0 #Sin hipertensión
ensanut2022$bp_diag[ensanut2022$HX_HBP == 1] <- 1 #Diagnosticado
ensanut2022$bp_diag[ensanut2022$HX_HBP == 0 & (ensanut2022$SBP >= 140 | ensanut2022$DBP >= 90)] <- 2 #No diagnosticado
ensanut2022$bp_diag <- factor(ensanut2022$bp_diag, labels = c("Without hypertension", "Diagnosed hypertension", "Undiagnosed hypertension"))

#Primary and secondary prevention
ensanut2022$prevent[ensanut2022$edad <40] <- 0
ensanut2022$prevent[ensanut2022$edad >= 40 & ensanut2022$prev_cv == 0] <- 1
ensanut2022$prevent[ensanut2022$edad >= 40 & ensanut2022$prev_cv == 1] <- 2
ensanut2022$prevent <- factor(ensanut2022$prevent, labels = c("<40 years", "Elegible for primary prevention", "Elegible for secondary prevention"))

#Hypertension treatment (regardless of hypertension diagnosis)
ensanut2022$tx_bp_all <- if_else((ensanut2022$bp_diag == "Without hypertension" & is.na(ensanut2022$TX_HBP)), "Untreated", ensanut2022$TX_HBP)
ensanut2022$tx_bp_all <- if_else((ensanut2022$bp_diag == "Undiagnosed hypertension" & is.na(ensanut2022$TX_HBP)), "Untreated", ensanut2022$tx_bp_all)
ensanut2022$tx_bp_all <- factor(ensanut2022$tx_bp_all)

#Diabetes complications
ensanut2022$diabetes_retinopatia<-if_else((ensanut2022$a0316c ==1 | ensanut2022$a0316d==1)==T,"Retinopathy", "No retinopathy")
ensanut2022$diabetes_nefropatia<-if_else(ensanut2022$a0316e==1,"Nephropathy", "No nephropathy")
ensanut2022$diabetes_neuropatia<-if_else(ensanut2022$a0316a==1,"Neuropathy", "No neuropathy")
ensanut2022$diabetes_anycomplication<-ifelse((ensanut2022$diabetes_retinopatia=="Retinopathy"|
                                                ensanut2022$diabetes_nefropatia=="Nephropathy"|
                                                ensanut2022$diabetes_neuropatia=="Neuropathy"), "Yes", "No")

ensanut2022$diabetes_allcomplications<-ifelse((ensanut2022$diabetes_retinopatia=="Retinopathy" &
                                                 ensanut2022$diabetes_nefropatia=="Nephropathy" &
                                                 ensanut2022$diabetes_neuropatia=="Neuropathy"), "Yes", "No")

#Cardiovascular outcomes and risk
ensanut2022$ascvd<-ifelse((ensanut2022$a0502a==1 | ensanut2022$a0502b==1)==T, "ASCVD", "No ASCVD")
ensanut2022$ethnicity <- rep(0, nrow(ensanut2022))
ensanut2022$eGFR<- with(ensanut2022, CKDEpi.creat(
  creatinine=Creatinine, sex=Sex, age=Age, ethnicity=ethnicity))
ensanut2022$severe_ckd<- ifelse(ensanut2022$eGFR<30, "eGFR <30mL/min/1.73m2", "eGFR ≥30mL/min/1.73m2")
ensanut2022$Risk.region<-"Moderate"
ensanut2022$total.chol<-ensanut2022$valor_COLEST/38.67
ensanut2022$total.hdl<-ensanut2022$valor_COL_HDL/38.67
ensanut2022$HbA1c<-10.929*(ensanut2022$valor_HB1AC-2.15)
ensanut2022$systolic.bp<-ensanut2022$SBP
ensanut2022$diabetes.age<-as.numeric(ensanut2022$edad_diag)
ensanut2022$diabetes<-ensanut2022$diag
ensanut2022$Gender<-as.character(factor(ensanut2022$male1, labels = c("female", "male")))
ensanut2022$smoker<-ifelse(ensanut2022$Smoking=="Current smoker",1,0)
ensanut2022$classify<-TRUE
ensanut2022$SCORE2_diabetes<-score_diab(ensanut2022, FALSE)
ensanut2022$SCORE2_diabetes_cat<-score_diab(ensanut2022, TRUE)
ensanut2022$SCORE2_diabetes_cat2<-ensanut2022$SCORE2_diabetes_cat
ensanut2022$SCORE2_diabetes_cat2[ensanut2022$ascvd=="ASCVD"]<-"Very high risk"
ensanut2022$SCORE2_diabetes_cat2[ensanut2022$diabetes_allcomplications=="Yes"]<-"Very high risk"
ensanut2022$high_risk<-ifelse(ensanut2022$SCORE2_diabetes_cat2 %in% c("High risk", "Very high risk"), 1, 0)

## LDL guidelines
ensanut2022$high_ldl<-ifelse(((ensanut2022$SCORE2_diabetes_cat2=="Very high risk" & ensanut2022$s_ldl<55)|
                                      (ensanut2022$SCORE2_diabetes_cat2=="High risk" & ensanut2022$s_ldl<70) |
                                      (ensanut2022$SCORE2_diabetes_cat2 %in% c("Low risk", "Moderate risk") & ensanut2022$s_ldl<100)), 1, 0)
ensanut2022$high_ldl<-factor(ensanut2022$high_ldl, labels = c("Uncontrolled", "coNTROLLED"))

#### Variables ENSANUT 2023####
ensanut2023 <- ensanut2023_fin %>% mutate(
  Year="2023",
  male = case_when(h0302.x == 1 ~ 1,
                   h0302.x == 2 ~ 0),
  H0310A = as.numeric(H0310A),
  H0310B = as.numeric(H0310B),
  imss = case_when(H0310A == 1 ~ 1,
                   H0310A != 1 ~ 0),
  issste = case_when(H0310A == 2 |  H0310A == 3 | H0310B == 2 | H0310B == 3 ~ 1,
                     H0310A != 2 |  H0310A != 3 | H0310B != 2 | H0310B != 3 ~ 0),
  pemex = case_when(H0310A == 4 | H0310B == 4 ~ 1,
                    H0310A != 4 | H0310B != 4 ~ 0),
  defensa = case_when(H0310A == 5 | H0310B == 5 ~ 1,
                      H0310A != 5 | H0310B != 5 ~ 0),
  marina = case_when(H0310A == 6 | H0310B == 6 ~ 1,
                     H0310A != 6 | H0310B != 6 ~ 0),
  imss_bienestar = case_when(H0310A == 7 | H0310B == 7 ~ 1,
                             H0310A != 7 | H0310B != 7 ~ 0),
  seguro_privado = case_when(H0310A == 8 | H0310B == 8 ~ 1,
                             H0310A != 8 | H0310B != 8 ~ 0),
  sin_seguro = case_when(H0310A == 10 | H0310B == 10 ~ 1,
                         H0310A != 10 | H0310B != 10 ~ 0),
  education = case_when(h0317a == 0 | h0317a == 1 ~ "No education",
                        h0317a == 2 ~ "Elementary school",
                        h0317a == 3 | h0317a == 4 ~ "Middle/High school",
                        h0317a == 9 | h0317a == 10 | h0317a == 11 | h0317a == 12 ~ "University",
                        h0317a == 5 |h0317a == 6 | h0317a == 7 | h0317a == 8 ~ "Other") %>% 
    factor(levels = c("No education", "Elementary school", "Middle/High school", "University", "Other")),
  #Weight control
  wt_loss = if_else(a0107 == 2, 1, 0), #1 = Perdió peso, 0 = No perdió peso
  intent_wt_loss = case_when(a0109 == 1 ~ 1, #1 = Sí, 0 = No
                             a0109 == 2 ~ 0),
  wt_loss_qt = a0108,
  #Diabetes diagnosis and management
  HX_T2D_AGE = a0302,
  plan_alim = case_when(A0312A == 1 ~ 1, 
                        A0312A %in% c(2:5) ~ 0) %>% factor(labels = c("No meal plan", "Meal plan")),
  plan_ejerc = case_when(A0312A == 2 | A0312B == 2 | A0312C == 2 | A0312D == 2 ~ 1,
                         A0312A != 2 | A0312B != 2 | A0312C != 2 | A0312D != 2 ~ 0) %>% factor(labels = c("No excercise plan", "Excercise plan")),
  aspirin = case_when(A0315A == 2 | A0315B == 2 ~ 1, 
                      A0315A != 2 | A0315B != 2  ~ 0),
  smoke_quit = case_when(Smoking %in% c(0, 1) ~ "Non-smoker",
                         Smoking == 2 ~ "Smoker") %>% factor(),
  #Hypertension diagnosis and management
  HX_HBP_AGE = Age - a0402a,
  salt = case_when(A0408A == 3 | A0408B == 3 | A0408C == 3 ~ 1, 
                   A0408A != 3 | A0408B != 3 | A0408C != 3  ~ 0) %>% factor(labels = c("No salt reduction", "Salt reduction")),
  #Cholesterol management
  chol_diag = if_else(a0604 == 1, 1, 0), #1 = Diagnóstico médico de hipercolesterolemia, 0 = Sin diagnóstico médico de hipercolesterolemia
  statin = case_when(A0605A == 1 ~ 1, 
                     A0605A != 1 ~ 0) %>% factor(labels = c("No statin", "Statin")),
  fibrates = case_when(A0607A1 == 2  ~ 1,
                       A0607A1 != 2 ~ 0) %>% factor(labels = c("No fibrates", "Fibrates")),
  #Previous cardiovascular disease
  prev_cv = case_when(a0502a == 1 | a0502b == 1 | a0502c == 1  ~ 1, #1 = Con CVD previa, 0 = Sin CVD
                      a0502a == 2 | a0502b == 2  | a0502c == 2 ~ 0), 
  #s_ldl
  s_ldl= (COLEST/0.948) -(COL_HDL/9.971)-((TRIG/8.56)+((TRIG*(COLEST-COL_HDL))/2140)- ((TRIG*TRIG)/16100))-9.4 
)

#Age
ensanut2023$age_group <- cut(as.numeric(ensanut2023$Age), breaks = c(20, 40, 60, Inf), right = FALSE, labels = c("<40 years", "40-59 years", "≥60 years"))

#Sex
ensanut2023$male1 <- factor(ensanut2023$male, labels = c("Women", "Men"))

#Diabetes
ensanut2023$diab[ensanut2023$HX_T2D == 1 | ensanut2023$HBA1C >= 6.5 | ensanut2023$TX_T2D %in% c(1, 2, 3)] <- 1

#Diagnosed and undiagnosed diabetes
ensanut2023$diag <- 0
ensanut2023$diag[ensanut2023$a0301 == 1 | ensanut2023$TX_T2D %in% c(1, 2, 3)] <- 1
ensanut2023$diag[ensanut2023$diag == 0 & ensanut2023$HBA1C >= 6.5] <- 2

#Diabetes treatment
ensanut2023$TX_T2D <- factor(ensanut2023$TX_T2D, labels = c("Insulin", "Pills", "Both", "None"))
ensanut2023$TX_T2D <- factor(ensanut2023$TX_T2D, levels = c("None", "Pills", "Insulin", "Both"))

#Glycemic control
ensanut2023$gly_control[(ensanut2023$h0303.x < 65 & ensanut2023$HBA1C < 7) | (ensanut2023$h0303.x >= 65 & ensanut2023$HBA1C < 7.5)] <- 1
ensanut2023$gly_control[(ensanut2023$h0303.x < 65 & ensanut2023$HBA1C >= 7) | (ensanut2023$h0303.x >= 65 & ensanut2023$HBA1C >= 7.5)] <- 0
ensanut2023$gly_control <- factor(ensanut2023$gly_control, labels = c("Uncontrolled", "Controlled"))

#Hypertension treatment
ensanut2023$TX_HBP <- factor(ensanut2023$TX_HBP, labels = c("Untreated", "Treated"))

#Uncontrolled hypertension
ensanut2023$high_bp[ensanut2023$SBP >= 130 | ensanut2023$DBP >= 80] <- 1
ensanut2023$high_bp[ensanut2023$SBP < 130 & ensanut2023$DBP < 80] <- 0
ensanut2023$high_bp <- factor(ensanut2023$high_bp, labels = c("Controlled", "Uncontrolled"))

#Uncontrolled LDL
ensanut2023$high_ldl2[ensanut2023$LDL >= 70] <- 1
ensanut2023$high_ldl2[ensanut2023$LDL < 70] <- 0
ensanut2023$high_ldl2 <- factor(ensanut2023$high_ldl2, labels = c("<70 mg/dL", "≥70 mg/dL"))

#Uncontrolled S-LDL
ensanut2023$high_sldl[ensanut2023$s_ldl >= 70] <- 1
ensanut2023$high_sldl[ensanut2023$s_ldl < 70] <- 0
ensanut2023$high_sldl <- factor(ensanut2023$high_sldl, labels = c("<70 mg/dL", "≥70 mg/dL"))

#Uncontrolled S-LDL (>100 mg/dL)
ensanut2023$high_sldl100[ensanut2023$s_ldl >= 100] <- 1
ensanut2023$high_sldl100[ensanut2023$s_ldl < 100] <- 0
ensanut2023$high_sldl100 <- factor(ensanut2023$high_sldl100, labels = c("<100 mg/dL", "≥100 mg/dL"))

#Uncontrolled HDL
ensanut2023$high_hdl[(ensanut2023$HDL <40 & ensanut2023$male == 1) | (ensanut2023$HDL <50 & ensanut2023$male == 0)] <- 1
ensanut2023$high_hdl[(ensanut2023$HDL >=40 & ensanut2023$male == 1) | (ensanut2023$HDL >=50 & ensanut2023$male == 0)] <- 0
ensanut2023$high_hdl <- factor(ensanut2023$high_hdl, labels = c("Controlled", "Uncontrolled"))

#Smoking status
ensanut2023$Smoking <- factor(ensanut2023$Smoking, labels = c("Never smoker", "Former smoker", "Current smoker"))

#Glycemic, blood pressure and cholesterol control
ensanut2023$comb_abc[ensanut2023$gly_control == "Controlled" & ensanut2023$high_bp == "Controlled" & ensanut2023$high_sldl == "<70 mg/dL"] <- 1
ensanut2023$comb_abc[ensanut2023$gly_control == "Uncontrolled" | ensanut2023$high_bp == "Uncontrolled" | ensanut2023$high_sldl == "≥70 mg/dL"] <- 0
ensanut2023$comb_abc <- factor(ensanut2023$comb_abc, labels = c("Uncontrolled", "Control"))

#Blood pressure and cholesterol control
ensanut2023$comb_bc[ensanut2023$high_bp == "Controlled" & ensanut2023$high_sldl == "<70 mg/dL"] <- 1
ensanut2023$comb_bc[ensanut2023$high_bp == "Uncontrolled" | ensanut2023$high_sldl == "≥70 mg/dL"] <- 0
ensanut2023$comb_bc <- factor(ensanut2023$comb_bc, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure, cholesterol control and no smoker
ensanut2023$comb_abcn[ensanut2023$gly_control == "Uncontrolled" | ensanut2023$high_bp == "Uncontrolled" | ensanut2023$high_sldl == "≥70 mg/dL" | ensanut2023$smoke_quit == "Smoker"] <- 0
ensanut2023$comb_abcn[ensanut2023$gly_control == "Controlled" & ensanut2023$high_bp == "Controlled" & ensanut2023$high_sldl == "<70 mg/dL" & ensanut2023$smoke_quit == "Non-smoker"] <- 1
ensanut2023$comb_abcn <- factor(ensanut2023$comb_abcn, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure and cholesterol control (LDL-C <100 mg/dL)
ensanut2023$comb_abc100[ensanut2023$gly_control == "Controlled" & ensanut2023$high_bp == "Controlled" & ensanut2023$high_sldl100 == "<100 mg/dL"] <- 1
ensanut2023$comb_abc100[ensanut2023$gly_control == "Uncontrolled" | ensanut2023$high_bp == "Uncontrolled" | ensanut2023$high_sldl100 == "Uncontrolled"] <- 0
ensanut2023$comb_abc100 <- factor(ensanut2023$comb_abc100, labels = c("Uncontrolled", "Control"))

#Blood pressure and cholesterol control (LDL-C <100 mg/dL)
ensanut2023$comb_bc100[ensanut2023$high_bp == "Controlled" & ensanut2023$high_sldl100 == "<100 mg/dL"] <- 1
ensanut2023$comb_bc100[ensanut2023$high_bp == "Uncontrolled" | ensanut2023$high_sldl100 == "Uncontrolled"] <- 0
ensanut2023$comb_bc100 <- factor(ensanut2023$comb_bc100, labels = c("Uncontrolled", "Control"))

#Glycemic, blood pressure, cholesterol control and no smoker (LDL-C <100 mg/dL)
ensanut2023$comb_abcn100[ensanut2023$gly_control == "Uncontrolled" | ensanut2023$high_bp == "Uncontrolled" | ensanut2023$high_sldl100 == "≥100 mg/dL" | ensanut2023$smoke_quit == "Smoker"] <- 0
ensanut2023$comb_abcn100[ensanut2023$gly_control == "Controlled" & ensanut2023$high_bp == "Controlled" & ensanut2023$high_sldl100 == "<100 mg/dL" & ensanut2023$smoke_quit == "Non-smoker"] <- 1
ensanut2023$comb_abcn100 <- factor(ensanut2023$comb_abcn100, labels = c("Uncontrolled", "Control"))

#Type of settlement (rural or urban)
ensanut2023$urban <- factor(ensanut2023$Area, labels = c("Rural", "Urban"))

#Indigenous
ensanut2023$indigenous_fct <- factor(ensanut2023$Indigenous, labels = c("Non-indigenous", "Indigenous"))

#Social security
ensanut2023$social_sec[ensanut2023$imss == 1 | ensanut2023$issste == 1 | ensanut2023$pemex == 1 | ensanut2023$defensa == 1 | ensanut2023$marina == 1 | ensanut2023$imss_bienestar == 1 | ensanut2023$seguro_privado == 1] <- 1
ensanut2023$social_sec[ensanut2023$sin_seguro == 1] <- 0
ensanut2023$social_sec <- factor(ensanut2023$social_sec, labels = c("Without social security", "With social security"))

#Geographic regions
ensanut2023$entidad.x <- as.numeric(ensanut2023$entidad.x)
ensanut2023$SUBREGION[ensanut2023$entidad.x %in% c(2,3,18,25,26)]<-1
ensanut2023$SUBREGION[ensanut2023$entidad.x %in% c(5,8,19,28)]<-2
ensanut2023$SUBREGION[ensanut2023$entidad.x %in% c(6,14,16)]<-3
ensanut2023$SUBREGION[ensanut2023$entidad.x %in% c(1,10,11,22,24,32)]<-4
ensanut2023$SUBREGION[ensanut2023$entidad.x %in% c(13,29,30)]<-5
ensanut2023$SUBREGION[ensanut2023$entidad.x %in% c(9)]<-6
ensanut2023$SUBREGION[ensanut2023$entidad.x %in% c(15)]<-7
ensanut2023$SUBREGION[ensanut2023$entidad.x %in% c(12,17,20,21)]<-8
ensanut2023$SUBREGION[ensanut2023$entidad.x %in% c(4,7,27,23,31)]<-9

#Diagnostic time
ensanut2023$edad_diag <- if_else(ensanut2023$a0302 == 99, NA, ensanut2023$a0302)
ensanut2023$diag_time <- as.numeric(ensanut2023$Age) - as.numeric(ensanut2023$edad_diag)
ensanut2023$diag_time <- if_else(ensanut2023$diag_time < 0, NA, ensanut2023$diag_time)
ensanut2023$cronic_diag[ensanut2023$diag_time >= 10] <- 1
ensanut2023$cronic_diag[ensanut2023$diag_time < 10 ] <- 0
ensanut2023$cronic_diag <- factor(ensanut2023$cronic_diag, labels = c("<10 years of diagnosis", "≥10 years of diagnosis"))

#High blood pressure diagnosis
ensanut2023$bp_diag[ensanut2023$HX_HBP == 0 & (ensanut2023$SBP < 140 & ensanut2023$DBP < 90)] <- 0 #Sin hipertensión
ensanut2023$bp_diag[ensanut2023$HX_HBP == 1] <- 1 #Diagnosticado
ensanut2023$bp_diag[ensanut2023$HX_HBP == 0 & (ensanut2023$SBP >= 140 | ensanut2023$DBP >= 90)] <- 2 #No diagnosticado
ensanut2023$bp_diag <- factor(ensanut2023$bp_diag, labels = c("Without hypertension", "Diagnosed hypertension", "Undiagnosed hypertension"))

#Primary and secondary prevention
ensanut2023$prevent[ensanut2023$Age <40] <- 0
ensanut2023$prevent[ensanut2023$Age >= 40 & ensanut2023$prev_cv == 0] <- 1
ensanut2023$prevent[ensanut2023$Age >= 40 & ensanut2023$prev_cv == 1] <- 2
ensanut2023$prevent <- factor(ensanut2023$prevent, labels = c("<40 years", "Elegible for primary prevention", "Elegible for secondary prevention"))

#Hypertension treatment (regardless of hypertension diagnosis)
ensanut2023$tx_bp_all <- if_else((ensanut2023$bp_diag == "Without hypertension" & is.na(ensanut2023$TX_HBP)), "Untreated", ensanut2023$TX_HBP)
ensanut2023$tx_bp_all <- if_else((ensanut2023$bp_diag == "Undiagnosed hypertension" & is.na(ensanut2023$TX_HBP)), "Untreated", ensanut2023$tx_bp_all)
ensanut2023$tx_bp_all <- factor(ensanut2023$tx_bp_all)

#Diabetes complications
ensanut2023$diabetes_retinopatia<-if_else((ensanut2023$a0316l ==1 | ensanut2023$a0316d==1)==T,"Retinopathy", "No retinopathy")
ensanut2023$diabetes_nefropatia<-if_else(ensanut2023$a0316e==1,"Nephropathy", "No nephropathy")
ensanut2023$diabetes_neuropatia<-if_else(ensanut2023$a0316k==1,"Neuropathy", "No neuropathy")
ensanut2023$diabetes_anycomplication<-ifelse((ensanut2023$diabetes_retinopatia=="Retinopathy"|
                                                ensanut2023$diabetes_nefropatia=="Nephropathy"|
                                                ensanut2023$diabetes_neuropatia=="Neuropathy"), "Yes", "No")

ensanut2023$diabetes_allcomplications<-ifelse((ensanut2023$diabetes_retinopatia=="Retinopathy" &
                                                 ensanut2023$diabetes_nefropatia=="Nephropathy" &
                                                 ensanut2023$diabetes_neuropatia=="Neuropathy"), "Yes", "No")

#Cardiovascular outcomes and risk
ensanut2023$ascvd<-ifelse((ensanut2023$a0502a==1 | ensanut2023$a0502b==1)==T, "ASCVD", "No ASCVD")
ensanut2023$ethnicity <- rep(0, nrow(ensanut2023))
ensanut2023$eGFR<- with(ensanut2023, CKDEpi.creat(
  creatinine=Creatinine, sex=Sex, age=Age, ethnicity=ethnicity))
ensanut2023$severe_ckd<- ifelse(ensanut2023$eGFR<30, "eGFR <30mL/min/1.73m2", "eGFR ≥30mL/min/1.73m2")
ensanut2023$Risk.region<-"Moderate"
ensanut2023$total.chol<-ensanut2023$COLEST/38.67
ensanut2023$total.hdl<-ensanut2023$COL_HDL/38.67
ensanut2023$HbA1c<-10.929*(ensanut2023$HB1AC-2.15)
ensanut2023$systolic.bp<-ensanut2023$SBP
ensanut2023$diabetes.age<-as.numeric(ensanut2023$edad_diag)
ensanut2023$diabetes<-ensanut2023$diag
ensanut2023$Gender<-as.character(factor(ensanut2023$male1, labels = c("female", "male")))
ensanut2023$smoker<-ifelse(ensanut2023$Smoking=="Current smoker",1,0)
ensanut2023$classify<-TRUE
ensanut2023$SCORE2_diabetes<-score_diab(ensanut2023, FALSE)
ensanut2023$SCORE2_diabetes_cat<-score_diab(ensanut2023, TRUE)
ensanut2023$SCORE2_diabetes_cat2<-ensanut2023$SCORE2_diabetes_cat
ensanut2023$SCORE2_diabetes_cat2[ensanut2023$ascvd=="ASCVD"]<-"Very high risk"
ensanut2023$SCORE2_diabetes_cat2[ensanut2023$diabetes_allcomplications=="Yes"]<-"Very high risk"
ensanut2023$high_risk<-ifelse(ensanut2023$SCORE2_diabetes_cat2 %in% c("High risk", "Very high risk"), 1, 0)

table(ensanut2023$high_risk, ensanut2023$SCORE2_diabetes_cat2)
## LDL guidelines
ensanut2023$high_ldl<-ifelse(((ensanut2023$SCORE2_diabetes_cat2=="Very high risk" & ensanut2023$s_ldl<55)|
                                      (ensanut2023$SCORE2_diabetes_cat2=="High risk" & ensanut2023$s_ldl<70) |
                                      (ensanut2023$SCORE2_diabetes_cat2 %in% c("Low risk", "Moderate risk") & ensanut2023$s_ldl<100)), 1, 0)

ensanut2023$high_ldl<-factor(ensanut2023$high_ldl, labels = c("Uncontrolled", "coNTROLLED"))

#### Missing values####
#Year 2016
na_2016 <- ensanut2016 %>% 
  filter(!is.na(ponde_f_vv)) %>% 
  select(HBA1C, SBP, DBP, LDL, s_ldl, HDL, TG, gly_control, high_bp, high_ldl, high_sldl, high_hdl, smoke_quit, cronic_diag)

vis_dat(na_2016)
vis_miss(na_2016)
miss_var_summary(na_2016)
gg_miss_var(na_2016, show_pct = TRUE) + labs(title = paste0("ENSANUT 2016 (n = ", nrow(na_2016), ")"))

#Year 2018
na_2018 <- ensanut2018 %>% 
  filter(!is.na(ponderador_glucosa)) %>% 
  select(HBA1C, SBP, DBP, LDL, s_ldl, HDL, TG, gly_control, high_bp, high_ldl, high_sldl, high_hdl, smoke_quit)

vis_dat(na_2018)
vis_miss(na_2018)
miss_var_summary(na_2018)
gg_miss_var(na_2018, show_pct = TRUE) + labs(title = paste0("ENSANUT 2018 (n = ", nrow(na_2018), ")"))

#Year 2021
na_2021 <- ensanut2021 %>% 
  filter(!is.na(ponde_vv)) %>% 
  select(HBA1C, SBP, DBP, LDL, s_ldl, HDL, TG, gly_control, high_bp, high_ldl, high_sldl, high_hdl, smoke_quit, cronic_diag)

vis_dat(na_2021)
vis_miss(na_2021)
miss_var_summary(na_2021)
gg_miss_var(na_2021, show_pct = TRUE) + labs(title = paste0("ENSANUT 2021 (n = ", nrow(na_2021), ")"))

#Year 2022
na_2022 <- ensanut2022 %>% 
  filter(!is.na(ponde_v)) %>% 
  select(HBA1C, SBP, DBP, LDL, s_ldl, HDL, TG, gly_control, high_bp, high_ldl, high_sldl, high_hdl, smoke_quit, cronic_diag)

vis_dat(na_2022)
vis_miss(na_2022)
miss_var_summary(na_2022)
gg_miss_var(na_2022, show_pct = TRUE) + labs(title = paste0("ENSANUT 2022 (n = ", nrow(na_2022), ")"))

#Year 2023
na_2023 <- ensanut2023 %>% 
  filter(!is.na(ponde_suero)) %>% 
  select(HBA1C, SBP, DBP, LDL, s_ldl, HDL, TG, gly_control, high_bp, high_ldl, high_sldl, high_hdl, smoke_quit, cronic_diag)

vis_dat(na_2023)
vis_miss(na_2023)
miss_var_summary(na_2023)
gg_miss_var(na_2023, show_pct = TRUE) + labs(title = paste0("ENSANUT 2022 (n = ", nrow(na_2022), ")"))


#### Survey objects - Year 2016####
#Lonely PSU
options(survey.adjust.domain.lonely = TRUE, survey.lonely.psu = "adjust")

ensanut2016_pre <- ensanut2016 %>% mutate(id = row_number()) %>% 
  filter(edad.x >= 20) %>% 
  filter(!is.na(est_var), !is.na(ponde_f_vv))

#Survey design
ensanut2016_survey <- svydesign(id = ~ code_upm.x, 
                                strata= ~ est_var.x, 
                                weights = ~ ponde_f_vv,  
                                nest = TRUE, data = ensanut2016_pre)

#Subpopulations of interest
#Individuals with diagnosed diabetes
ensanut2016_diag <- subset(ensanut2016_survey, diag == 1)

#Individuals with hypertension
ensanut2016_bp <- subset(ensanut2016_survey, (diag == 1 & bp_diag == "Diagnosed hypertension"))

#Individuals elegible for primary prevention
ensanut2016_primary <- subset(ensanut2016_survey, (diag == 1 & prevent == "Elegible for primary prevention"))

#Individuals elegible for secondary prevention
ensanut2016_secondary <- subset(ensanut2016_survey, (diag == 1  & prevent == "Elegible for secondary prevention"))

#Individuals with SCORE2-Diabetes
ensanut2016_score <- subset(ensanut2016_survey, (diag == 1 & !is.na(SCORE2_diabetes)))

#### Survey objects - Year 2018####
ensanut2018_pre <- ensanut2018 %>% mutate(id = row_number()) %>% 
  filter(EDAD.x >= 20) %>% 
  filter(!is.na(ponderador_glucosa), !is.na(EST_DIS.x))

#Survey design
ensanut2018_survey <- svydesign(id = ~ UPM.x, 
                                strata = ~ EST_DIS.x, #Cuál? EST_DIS.x
                                weights = ~ ponderador_glucosa, 
                                nest = TRUE, data = ensanut2018_pre,)

#Subpopulations of interest
#Individuals with diagnosed diabetes
ensanut2018_diag <- subset(ensanut2018_survey, diag == 1)
ensanut2018_diag2 <- subset(ensanut2018_survey, diag == 1 & between(Age, 30, 69))

#Individuals with hypertension
ensanut2018_bp <- subset(ensanut2018_survey, (diag == 1 &  bp_diag == "Diagnosed hypertension"))

#Individuals elegible for primary prevention
ensanut2018_primary <- subset(ensanut2018_survey, (diag == 1 & prevent == "Elegible for primary prevention"))

#Individuals elegible for secondary prevention
ensanut2018_secondary <- subset(ensanut2018_survey, (diag == 1 & prevent == "Elegible for secondary prevention"))

ensanut2018_elegible <- subset(ensanut2018_survey, (diag == 1 & between(Age, 30, 69) & prevent %in% c("Elegible for primary prevention","Elegible for secondary prevention")))


#Individuals with SCORE2-Diabetes
ensanut2018_score <- subset(ensanut2018_survey, (diag == 1 & !is.na(SCORE2_diabetes)))

#### Survey objects - Year 2021####
ensanut2021_pre <- ensanut2021 %>% mutate(id = row_number()) %>% 
  filter(edad >= 20) %>% 
  filter(!is.na(est_final), !is.na(ponde_vv))

#Survey design
ensanut2021_survey <- svydesign(id = ~ upm.x, 
                                strata = ~ est_final,#Cuál? est_final
                                weights = ~ ponde_vv, 
                                nest = TRUE, data = ensanut2021_pre)

#Subpopulations of interest
#Individuals with diagnosed diabetes
ensanut2021_diag <- subset(ensanut2021_survey, diag == 1)

#Individuals with hypertension
ensanut2021_bp <- subset(ensanut2021_survey, (diag == 1  & bp_diag == "Diagnosed hypertension"))

#Individuals elegible for primary prevention
ensanut2021_primary <- subset(ensanut2021_survey, (diag == 1 & prevent == "Elegible for primary prevention"))

#Individuals elegible for primary prevention
ensanut2021_secondary <- subset(ensanut2021_survey, (diag == 1 & prevent == "Elegible for secondary prevention"))

#Individuals with SCORE2-Diabetes
ensanut2021_score <- subset(ensanut2021_survey, (diag == 1 & !is.na(SCORE2_diabetes)))

#### Survey objects - Year 2022####
ensanut2022_pre <- ensanut2022 %>% mutate(id = row_number()) %>% 
  filter(edad >= 20) %>% 
  filter(!is.na(ponde_v), !is.na(est_sel.x.x.x))

ensanut2022_survey <- svydesign(id = ~ upm.x, 
                                strata = ~ est_sel.x.x.x, 
                                weights = ~ ponde_v, 
                                nest = TRUE, data = ensanut2022_pre,)

#Subpopulations of interest
#Individuals with diagnosed diabetes
ensanut2022_diag <- subset(ensanut2022_survey, diag == 1)

#Individuals with hypertension
ensanut2022_bp <- subset(ensanut2022_survey, (diag == 1  & bp_diag == "Diagnosed hypertension"))

#Individuals elegible for primary prevention
ensanut2022_primary <- subset(ensanut2022_survey, (diag == 1  & prevent == "Elegible for primary prevention"))

#Individuals elegible for primary prevention
ensanut2022_secondary <- subset(ensanut2022_survey, (diag == 1  & prevent == "Elegible for secondary prevention"))

#Individuals with SCORE2-Diabetes
ensanut2022_score <- subset(ensanut2022_survey, (diag == 1 & !is.na(SCORE2_diabetes)))

#### Survey objects - Year 2023####
ensanut2023_pre <- ensanut2023 %>% mutate(id = row_number()) %>% 
  filter(h0303.x >= 20) %>% 
  filter(!is.na(ponde_suero_2), !is.na(est_sel.y.y))

ensanut2023_antro <- ensanut2023 %>% mutate(id = row_number()) %>% 
  filter(h0303.x >= 20) %>% 
  filter(!is.na(ponde_f), !is.na(est_sel.x.x))

ensanut2023_survey <- svydesign(id = ~ upm.y.y, 
                                strata = ~ est_sel.y.y, 
                                weights = ~ ponde_suero_2, 
                                nest = TRUE, data = ensanut2023_pre)

ensanut2023_survey_antro <- svydesign(id = ~ upm.x.x, 
                                strata = ~ est_sel.x.x, 
                                weights = ~ ponde_f, 
                                nest = TRUE, data = ensanut2023_antro)

sum(ensanut2023_survey_antro$variables$ponde_f)

#Subpopulations of interest
#Individuals with diagnosed diabetes
ensanut2023_diag <- subset(ensanut2023_survey, (diag == 1))

#Individuals with hypertension
ensanut2023_bp <- subset(ensanut2023_survey, (diag == 1  & bp_diag == "Diagnosed hypertension"))

#Individuals with hypertension
ensanut2023_bp3 <- subset(ensanut2023_survey, diag == 1)


#Individuals elegible for primary prevention
ensanut2023_primary <- subset(ensanut2023_survey, (diag == 1 & prevent == "Elegible for primary prevention"))

#Individuals elegible for secondary prevention
ensanut2023_secondary <- subset(ensanut2023_survey, (diag == 1 & prevent == "Elegible for secondary prevention"))

#Individuals with SCORE2-Diabetes
ensanut2023_score <- subset(ensanut2023_survey, (diag == 1 & !is.na(SCORE2_diabetes)))

#### Weights used for modeling####
#2016
ensanut2016$CVE_ENT <- as.numeric(ensanut2016$entidad.x)
ensanut2016$ponde_f_vv <- as.numeric(as.character(ensanut2016$ponde_f_vv)) #Ponderador o Factor de Expansión
ensanut2016$sampleid <- paste(ensanut2016$CVE_ENT, as.numeric(ensanut2016$est_var)) #Crear un ID con el estado y la variable de estratificacion

#Within each sampling unit, sum the weights
wts <- tapply(ensanut2016$ponde_f_vv,ensanut2016$sampleid,sum, na.rm = TRUE)

#Make a data frame from this
wts <- data.frame(id=names(unlist(wts)), wt=unlist(wts))

#Get the unique sampling location ids'
t1 <- as.data.frame(table(ensanut2016$sampleid))

#Put all of this into a data set
wts2 <- data.frame(ids=wts$id, sumwt=wts$wt, jn=t1$Freq)

#Merge all of this back to the original data file
ensanut2016 <- merge(ensanut2016, wts2, by.x="sampleid", by.y="ids", all.x=T)

#In the new data set, multiply the original weight by the fraction of the sampling unit total population each person represents
ensanut2016$swts <- ensanut2016$ponde_f_vv*(ensanut2016$jn/ensanut2016$sumwt)

#2018
ensanut2018$CVE_ENT <- as.numeric(ensanut2018$ENT.x)
ensanut2018$ponderador_glucosa <- as.numeric(as.character(ensanut2018$ponderador_glucosa)) #Ponderador o Factor de Expansión
ensanut2018$sampleid <- paste(ensanut2018$CVE_ENT, as.numeric(ensanut2018$ESTRATO.y)) #Crear un ID con el estado y la variable de estratificacion

wts <- tapply(ensanut2018$ponderador_glucosa,ensanut2018$sampleid,sum, na.rm = TRUE)
wts <- data.frame(id=names(unlist(wts)), wt=unlist(wts))
t1 <- as.data.frame(table(ensanut2018$sampleid))
wts2 <- data.frame(ids=wts$id, sumwt=wts$wt, jn=t1$Freq)
ensanut2018 <- merge(ensanut2018, wts2, by.x="sampleid", by.y="ids", all.x=T)
ensanut2018$swts <- ensanut2018$ponderador_glucosa*(ensanut2018$jn/ensanut2018$sumwt)

#2021
ensanut2021$CVE_ENT <- as.numeric(ensanut2021$entidad.x)
ensanut2021$ponde_vv <- as.numeric(as.character(ensanut2021$ponde_vv)) #Ponderador o Factor de Expansión
ensanut2021$sampleid <- paste(ensanut2021$CVE_ENT, as.numeric(ensanut2021$est_final)) #Crear un ID con el estado y la variable de estratificacion

wts <- tapply(ensanut2021$ponde_vv,ensanut2021$sampleid,sum, na.rm = TRUE)
wts <- data.frame(id=names(unlist(wts)), wt=unlist(wts))
t1 <- as.data.frame(table(ensanut2021$sampleid))
wts2 <- data.frame(ids=wts$id, sumwt=wts$wt, jn=t1$Freq)
ensanut2021 <- merge(ensanut2021, wts2, by.x="sampleid", by.y="ids", all.x=T)
ensanut2021$swts <- ensanut2021$ponde_vv*(ensanut2021$jn/ensanut2021$sumwt)

#2022
ensanut2022$CVE_ENT <- as.numeric(ensanut2022$entidad.x)
ensanut2022$ponde_v <- as.numeric(as.character(ensanut2022$ponde_v)) #Ponderador o Factor de Expansión
ensanut2022$sampleid <- paste(ensanut2022$CVE_ENT, as.numeric(ensanut2022$est_sel.x.x.x)) #Crear un ID con el estado y la variable de estratificacion

wts <- tapply(ensanut2022$ponde_v,ensanut2022$sampleid,sum, na.rm = TRUE)
wts <- data.frame(id=names(unlist(wts)), wt=unlist(wts))
t1 <- as.data.frame(table(ensanut2022$sampleid))
wts2 <- data.frame(ids=wts$id, sumwt=wts$wt, jn=t1$Freq)
ensanut2022 <- merge(ensanut2022, wts2, by.x="sampleid", by.y="ids", all.x=T)
ensanut2022$swts <- ensanut2022$ponde_v*(ensanut2022$jn/ensanut2022$sumwt)

#2023
ensanut2023$CVE_ENT <- as.numeric(ensanut2023$entidad.x)
ensanut2023$ponde_suero <- as.numeric(as.character(ensanut2023$ponde_suero)) #Ponderador o Factor de Expansión
ensanut2023$sampleid <- paste(ensanut2023$CVE_ENT, as.numeric(ensanut2023$estrato.y.y)) #Crear un ID con el estado y la variable de estratificacion

wts <- tapply(ensanut2023$ponde_suero,ensanut2023$sampleid,sum, na.rm = TRUE)
wts <- data.frame(id=names(unlist(wts)), wt=unlist(wts))
t1 <- as.data.frame(table(ensanut2023$sampleid))
wts2 <- data.frame(ids=wts$id, sumwt=wts$wt, jn=t1$Freq)
ensanut2023 <- merge(ensanut2023, wts2, by.x="sampleid", by.y="ids", all.x=T)
ensanut2023$swts <- ensanut2023$ponde_suero*(ensanut2023$jn/ensanut2023$sumwt)

#### Descriptive analysis####
#Year 2016
t_2016<-ensanut2016_diag %>% 
  tbl_svysummary(label = list(male1 ~ "Sex",
                              age_group ~ "Age group",
                              HBA1C ~ "HbA1c",
                              TX_T2D ~ "Diabetes treatment",
                              TX_HBP ~ "Hypertension treatment",
                              education ~ "Education level",
                              social_sec ~ "Social security",
                              urban ~ "Location",
                              indigenous_fct ~ "Indigenous identity",
                              bp_diag ~ "Diagnosed hypertension",
                              prevent ~ "CV risk prevention"),
                 statistic = list(all_categorical() ~ "{p}",
                                  all_continuous() ~ "{mean}"),
                 include = c(male1, age_group, HBA1C, SBP, DBP, TX_T2D, 
                             TX_HBP, education, social_sec, urban, indigenous_fct, 
                             bp_diag, prevent),
                 missing = "no") %>% 
  add_ci(pattern = "{stat} ({ci})") 

total_2016 <- svytotal(~ diag, ensanut2016_diag, na.rm = TRUE)

#Year 2018
t_2018<-ensanut2018_diag %>% 
  tbl_svysummary(label = list(male1 ~ "Sex",
                              age_group ~ "Age group",
                              HBA1C ~ "HbA1c",
                              TX_T2D ~ "Diabetes treatment",
                              TX_HBP ~ "Hypertension treatment",
                              education ~ "Education level",
                              social_sec ~ "Social security",
                              urban ~ "Location",
                              indigenous_fct ~ "Indigenous identity",
                              bp_diag ~ "Diagnosed hypertension",
                              prevent ~ "CV risk prevention"),
                 statistic = list(all_categorical() ~ "{p}",
                                  all_continuous() ~ "{mean}"),
                 include = c(male1, age_group, HBA1C, SBP, DBP, TX_T2D, 
                             TX_HBP, education, social_sec, urban, indigenous_fct,
                             bp_diag, prevent),
                 missing = "no") %>% 
  add_ci(pattern = "{stat} ({ci})") 

total_2018 <- svytotal(~ diag, ensanut2018_diag, na.rm = TRUE)

#Year 2021
t_2021<-ensanut2021_diag %>% 
  tbl_svysummary(label = list(male1 ~ "Sex",
                              age_group ~ "Age group",
                              HBA1C ~ "HbA1c",
                              TX_T2D ~ "Diabetes treatment",
                              TX_HBP ~ "Hypertension treatment",
                              education ~ "Education level",
                              social_sec ~ "Social security",
                              urban ~ "Location",
                              indigenous_fct ~ "Indigenous identity",
                              bp_diag ~ "Diagnosed hypertension",
                              prevent ~ "CV risk prevention"),
                 statistic = list(all_categorical() ~ "{p}",
                                  all_continuous() ~ "{mean}"),
                 include = c(male1, age_group, HBA1C, SBP, DBP, TX_T2D, 
                             TX_HBP, education, social_sec, urban, indigenous_fct,
                             bp_diag, prevent),
                 missing = "no") %>% 
  add_ci(pattern = "{stat} ({ci})") 

total_2021 <- svytotal(~ diag, ensanut2021_diag, na.rm = TRUE)

#Year 2022
t_2022<-ensanut2022_diag %>% 
  tbl_svysummary(label = list(male1 ~ "Sex",
                              age_group ~ "Age group",
                              HBA1C ~ "HbA1c",
                              TX_T2D ~ "Diabetes treatment",
                              TX_HBP ~ "Hypertension treatment",
                              education ~ "Education level",
                              social_sec ~ "Social security",
                              urban ~ "Location",
                              indigenous_fct ~ "Indigenous identity",
                              bp_diag ~ "Diagnosed hypertension",
                              prevent ~ "CV risk prevention"),
                 statistic = list(all_categorical() ~ "{p}",
                                  all_continuous() ~ "{mean}"),
                 include = c(male1, age_group, HBA1C, SBP, DBP, TX_T2D, 
                             TX_HBP, education, social_sec, urban, indigenous_fct,
                             bp_diag, prevent),
                 missing = "no") %>% 
  add_ci(pattern = "{stat} ({ci})") 

total_2022 <- svytotal(~ diag, ensanut2022_diag, na.rm = TRUE)

#Year 2023
t_2023<-ensanut2023_diag %>% 
  tbl_svysummary(label = list(male1 ~ "Sex",
                              age_group ~ "Age group",
                              HBA1C ~ "HbA1c",
                              TX_T2D ~ "Diabetes treatment",
                              TX_HBP ~ "Hypertension treatment",
                              education ~ "Education level",
                              social_sec ~ "Social security",
                              urban ~ "Location",
                              indigenous_fct ~ "Indigenous identity",
                              bp_diag ~ "Diagnosed hypertension",
                              prevent ~ "CV risk prevention"),
                 statistic = list(all_categorical() ~ "{p}",
                                  all_continuous() ~ "{mean}"),
                 include = c(male1, age_group, HBA1C, SBP, DBP, TX_T2D, 
                             TX_HBP, education, social_sec, urban, indigenous_fct,
                             bp_diag, prevent),
                 missing = "no") %>% 
  add_ci(pattern = "{stat} ({ci})") 

total_2023 <- svytotal(~ diag, ensanut2023_diag, na.rm = TRUE)
#gtsummary, flextable, officer

#Save table
tab2 <-tbl_merge(tbls = list(t_2016, t_2018, t_2021, t_2022, t_2023))%>% as_flex_table() %>% align(align = "center",part = "all") %>% autofit()
doc <- read_docx() %>% body_add_flextable(value = tab2, split = TRUE) %>%  body_end_section_landscape() %>%
  print(target = "Tables/tablaS1.docx"); remove(t_2016, t_2018, t_2021, t_2022, t_2023, tab2)

#### Variables of interest####
var_diag <- c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100","smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100")
var_nodiag <- c("gly_control", "high_bp", "high_ldl","high_sldl", "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "cronic_diag")
var_treat <- c("TX_T2D", "TX_HBP", "statin")

#### Prevalence estimates year 2016####
#Prevalence estimates of risk factor control among individuals with diagnosed diabetes
diab_2016_diag <- map(var_diag, wt_prop, ensanut2016_diag) %>% 
  set_names(var_diag)

#Prevalence estimates of treatment among individuals with diagnosed diabetes
diab_2016_treat <- map(var_treat, wt_prop, ensanut2016_diag) %>% 
  set_names(var_treat)

#Stratification by sex
list_sex_2016 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100",
                               "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("male1"),
                      design = list(ensanut2016_diag),
                      suffix = c("gly", "high", "high", "high", "high","smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2016)

sex_strat_2016 <- pmap(list_sex_2016, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by age
list_age_2016 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100",
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("age_group"),
                      design = list(ensanut2016_diag),
                      suffix = c("gly", "high", "high", "high", "high","smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2016)

age_strat_2016 <- pmap(list_age_2016, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by educational level
list_edu_2016 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("education"),
                      design = list(ensanut2016_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2016)

edu_strat_2016 <- pmap(list_edu_2016, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by type of settlement (rural or urban)
list_urb_2016 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("urban"),
                      design = list(ensanut2016_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2016)

urb_strat_2016 <- pmap(list_urb_2016, wt_prop_by) %>%
  set_names(var_diag)


#Stratification by indigenous identity
list_ind_2016 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("indigenous_fct"),
                      design = list(ensanut2016_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2016)

ind_strat_2016 <- pmap(list_ind_2016, wt_prop_by) %>%
  set_names(var_diag)


#Stratification by time since diagnosis
list_chr_2016 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("cronic_diag"),
                      design = list(ensanut2016_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2016)

chr_strat_2016 <- pmap(list_chr_2016, wt_prop_by) %>%
  set_names(var_diag)

#### Prevalence estimates year 2018####
#Prevalence estimates among individuals with diagnosed diabetes
diab_2018_diag <- map(var_diag, wt_prop, ensanut2018_diag) %>% 
  set_names(var_diag)

#Prevalence estimates of treatment among individuals with diagnosed diabetes
diab_2018_treat <- map(var_treat, wt_prop, ensanut2018_diag) %>% 
  set_names(var_treat)

#Stratification by sex
list_sex_2018 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("male1"),
                      design = list(ensanut2018_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2018)

sex_strat_2018 <- pmap(list_sex_2018, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by age
list_age_2018 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("age_group"),
                      design = list(ensanut2018_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2018)

age_strat_2018 <- pmap(list_age_2018, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by educational level
list_edu_2018 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("education"),
                      design = list(ensanut2018_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2018)

edu_strat_2018 <- pmap(list_edu_2018, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by type of settlement (rural or urban)
list_urb_2018 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("urban"),
                      design = list(ensanut2018_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2018)

urb_strat_2018 <- pmap(list_urb_2018, wt_prop_by) %>%
  set_names(var_diag)


#Stratification by indigenous identity
list_ind_2018 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("indigenous_fct"),
                      design = list(ensanut2018_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2018)

ind_strat_2018 <- pmap(list_ind_2018, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by time since diagnosis
list_chr_2018 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("cronic_diag"),
                      design = list(ensanut2018_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2018)

chr_strat_2018 <- pmap(list_chr_2018, wt_prop_by) %>%
  set_names(var_diag)

#### Prevalence estimates year 2021####
#Prevalence estimates among individuals with diagnosed diabetes
diab_2021_diag <- map(var_diag, wt_prop, ensanut2021_diag) %>% 
  set_names(var_diag)

#Prevalence estimates of treatment among individuals with diagnosed diabetes
diab_2021_treat <- map(var_treat, wt_prop, ensanut2021_diag) %>% 
  set_names(var_treat)

#Stratification by sex
list_sex_2021 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("male1"),
                      design = list(ensanut2021_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2021)

sex_strat_2021 <- pmap(list_sex_2021, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by age
list_age_2021 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("age_group"),
                      design = list(ensanut2021_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2021)

age_strat_2021 <- pmap(list_age_2021, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by educational level
list_edu_2021 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("education"),
                      design = list(ensanut2021_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2021)

edu_strat_2021 <- pmap(list_edu_2021, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by type of settlement (rural or urban)
list_urb_2021 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("urban"),
                      design = list(ensanut2021_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2021)

urb_strat_2021 <- pmap(list_urb_2021, wt_prop_by) %>%
  set_names(var_diag)


#Stratification by indigenous identity
list_ind_2021 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("indigenous_fct"),
                      design = list(ensanut2021_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2021)

ind_strat_2021 <- pmap(list_ind_2021, wt_prop_by) %>%
  set_names(var_diag)


#Stratification by time since diagnosis
list_chr_2021 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("cronic_diag"),
                      design = list(ensanut2021_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2021)

chr_strat_2021 <- pmap(list_chr_2021, wt_prop_by) %>%
  set_names(var_diag)

#### Prevalence estimates year 2022####
#Prevalence estimates among individuals with diagnosed diabetes
diab_2022_diag <- map(var_diag, wt_prop, ensanut2022_diag) %>% 
  set_names(var_diag)

#Prevalence estimates of treatment among individuals with diagnosed diabetes
diab_2022_treat <- map(var_treat, wt_prop, ensanut2022_diag) %>% 
  set_names(var_treat)

#Stratification by sex
list_sex_2022 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("male1"),
                      design = list(ensanut2022_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2022)

sex_strat_2022 <- pmap(list_sex_2022, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by age
list_age_2022 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("age_group"),
                      design = list(ensanut2022_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2022)

age_strat_2022 <- pmap(list_age_2022, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by educational level
list_edu_2022 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("education"),
                      design = list(ensanut2022_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2022)

edu_strat_2022 <- pmap(list_edu_2022, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by type of settlement (rural or urban)
list_urb_2022 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("urban"),
                      design = list(ensanut2022_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2022)

urb_strat_2022 <- pmap(list_urb_2022, wt_prop_by) %>%
  set_names(var_diag)


#Stratification by indigenous identity
list_ind_2022 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("indigenous_fct"),
                      design = list(ensanut2022_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2022)

ind_strat_2022 <- pmap(list_ind_2022, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by time since diagnosis
list_chr_2022 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("cronic_diag"),
                      design = list(ensanut2022_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2022)

chr_strat_2022 <- pmap(list_chr_2022, wt_prop_by) %>%
  set_names(var_diag)

#### Prevalence estimates year 2023####
#Prevalence estimates among individuals with diagnosed diabetes
diab_2023_diag <- map(var_diag, wt_prop, ensanut2023_diag) %>% 
  set_names(var_diag)

#Prevalence estimates of treatment among individuals with diagnosed diabetes
diab_2023_treat <- map(var_treat, wt_prop, ensanut2023_diag) %>% 
  set_names(var_treat)

#Stratification by sex
list_sex_2023 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("male1"),
                      design = list(ensanut2023_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2023)

sex_strat_2023 <- pmap(list_sex_2023, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by age
list_age_2023 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("age_group"),
                      design = list(ensanut2023_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2023)

age_strat_2023 <- pmap(list_age_2023, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by educational level
list_edu_2023 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("education"),
                      design = list(ensanut2023_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2023)

edu_strat_2023 <- pmap(list_edu_2023, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by type of settlement (rural or urban)
list_urb_2023 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("urban"),
                      design = list(ensanut2023_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2023)

urb_strat_2023 <- pmap(list_urb_2023, wt_prop_by) %>%
  set_names(var_diag)


#Stratification by indigenous identity
list_ind_2023 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("indigenous_fct"),
                      design = list(ensanut2023_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2023)

ind_strat_2023 <- pmap(list_ind_2023, wt_prop_by) %>%
  set_names(var_diag)

#Stratification by time since diagnosis
list_chr_2023 <- list(variable = c("gly_control", "high_bp", "high_ldl","high_sldl", "high_sldl100", 
                                   "smoke_quit", "comb_bc", "comb_abc", "comb_abcn", "comb_bc100", "comb_abc100", "comb_abcn100"),
                      by = c("cronic_diag"),
                      design = list(ensanut2023_diag),
                      suffix = c("gly", "high", "high", "high", "high", "smoke", "comb", "comb", "comb", "comb", "comb", "comb"),
                      year = 2023)

chr_strat_2023 <- pmap(list_chr_2023, wt_prop_by) %>%
  set_names(var_diag)


#### Prevalence estimates for every year####
all_years <- c(2016, 2018, 2021, 2022, 2023)

#Diagnosed diabetes
diag_glyc_total <- list_extr(years = all_years, variable = "gly_control", suffix = "_diag")
diag_high_bp_total <- list_extr(years = all_years, variable = "high_bp", suffix = "_diag")
diag_high_ldl_total <- list_extr(years = all_years, variable = "high_ldl", suffix = "_diag")
diag_shigh_ldl_total <- list_extr(years = all_years, variable = "high_sldl", suffix = "_diag")
diag_shigh_ldl100_total <- list_extr(years = all_years, variable = "high_sldl100", suffix = "_diag")
diag_smoking <- list_extr(years = all_years, variable = "smoke_quit", suffix = "_diag")
diag_shigh_ldl100_total <- list_extr(years = all_years, variable = "high_sldl100", suffix = "_diag")
diag_bc <- list_extr(years = all_years, variable = "comb_bc", suffix = "_diag")
diag_abc <- list_extr(years = all_years, variable = "comb_abc", suffix = "_diag")
diag_abcn <- list_extr(years = all_years, variable = "comb_abcn", suffix = "_diag")
diag_bc100 <- list_extr(years = all_years, variable = "comb_bc100", suffix = "_diag")
diag_abc100 <- list_extr(years = all_years, variable = "comb_abc100", suffix = "_diag")
diag_abcn100 <- list_extr(years = all_years, variable = "comb_abcn100", suffix = "_diag")
diag_tx_diab <- list_extr(years = all_years, variable = "TX_T2D", suffix = "_treat")
diag_tx_bp <- list_extr(years = all_years, variable = "TX_HBP", suffix = "_treat")
diag_tx_ldl <- list_extr(years = all_years, variable = "statin", suffix = "_treat")

#Stratified by sex
gly_sex_total <- rbind(sex_strat_2016$gly_control, sex_strat_2018$gly_control, sex_strat_2021$gly_control, sex_strat_2022$gly_control,sex_strat_2023$gly_control)
bp_sex_total <- rbind(sex_strat_2016$high_bp, sex_strat_2018$high_bp, sex_strat_2021$high_bp, sex_strat_2022$high_bp,sex_strat_2023$high_bp)
ldl_sex_total <- rbind(sex_strat_2016$high_ldl, sex_strat_2018$high_ldl, sex_strat_2021$high_ldl, sex_strat_2022$high_ldl,sex_strat_2023$high_ldl)
sldl_sex_total <- rbind(sex_strat_2016$high_sldl, sex_strat_2018$high_sldl, sex_strat_2021$high_sldl, sex_strat_2022$high_sldl,sex_strat_2023$high_sldl)
sldl100_sex_total <- rbind(sex_strat_2016$high_sldl100, sex_strat_2018$high_sldl100, sex_strat_2021$high_sldl100, sex_strat_2022$high_sldl100,sex_strat_2023$high_sldl100)
smoke_sex_total <- rbind(sex_strat_2016$smoke_quit, sex_strat_2018$smoke_quit, sex_strat_2021$smoke_quit, sex_strat_2022$smoke_quit,sex_strat_2023$smoke_quit)
bc_sex_total <- rbind(sex_strat_2016$comb_bc, sex_strat_2018$comb_bc, sex_strat_2021$comb_bc, sex_strat_2022$comb_bc,sex_strat_2023$comb_bc)
abc_sex_total <- rbind(sex_strat_2016$comb_abc, sex_strat_2018$comb_abc, sex_strat_2021$comb_abc, sex_strat_2022$comb_abc,sex_strat_2023$comb_abc)
abcn_sex_total <- rbind(sex_strat_2016$comb_abcn, sex_strat_2018$comb_abcn, sex_strat_2021$comb_abcn, sex_strat_2022$comb_abcn,sex_strat_2023$comb_abcn)
bc100_sex_total <- rbind(sex_strat_2016$comb_bc100, sex_strat_2018$comb_bc100, sex_strat_2021$comb_bc100, sex_strat_2022$comb_bc100,sex_strat_2023$comb_bc100)
abc100_sex_total <- rbind(sex_strat_2016$comb_abc100, sex_strat_2018$comb_abc100, sex_strat_2021$comb_abc100, sex_strat_2022$comb_abc100,sex_strat_2023$comb_abc100)
abcn100_sex_total <- rbind(sex_strat_2016$comb_abcn100, sex_strat_2018$comb_abcn100, sex_strat_2021$comb_abcn100, sex_strat_2022$comb_abcn100,sex_strat_2023$comb_abcn100)

#Stratified by age
gly_age_total <- rbind(age_strat_2016$gly_control, age_strat_2018$gly_control, age_strat_2021$gly_control, age_strat_2022$gly_control,age_strat_2023$gly_control)
bp_age_total <- rbind(age_strat_2016$high_bp, age_strat_2018$high_bp, age_strat_2021$high_bp, age_strat_2022$high_bp,age_strat_2023$high_bp)
ldl_age_total <- rbind(age_strat_2016$high_ldl, age_strat_2018$high_ldl, age_strat_2021$high_ldl, age_strat_2022$high_ldl,age_strat_2023$high_ldl)
sldl_age_total <- rbind(age_strat_2016$high_sldl, age_strat_2018$high_sldl, age_strat_2021$high_sldl, age_strat_2022$high_sldl,age_strat_2023$high_sldl)
sldl100_age_total <- rbind(age_strat_2016$high_sldl100, age_strat_2018$high_sldl100, age_strat_2021$high_sldl100, age_strat_2022$high_sldl100,age_strat_2023$high_sldl100)
smoke_age_total <- rbind(age_strat_2016$smoke_quit, age_strat_2018$smoke_quit, age_strat_2021$smoke_quit, age_strat_2022$smoke_quit,age_strat_2023$smoke_quit)
bc_age_total <- rbind(age_strat_2016$comb_bc, age_strat_2018$comb_bc, age_strat_2021$comb_bc, age_strat_2022$comb_bc,age_strat_2023$comb_bc)
abc_age_total <- rbind(age_strat_2016$comb_abc, age_strat_2018$comb_abc, age_strat_2021$comb_abc, age_strat_2022$comb_abc,age_strat_2023$comb_abc)
abcn_age_total <- rbind(age_strat_2016$comb_abcn, age_strat_2018$comb_abcn, age_strat_2021$comb_abcn, age_strat_2022$comb_abcn,age_strat_2023$comb_abcn)
bc100_age_total <- rbind(age_strat_2016$comb_bc100, age_strat_2018$comb_bc100, age_strat_2021$comb_bc100, age_strat_2022$comb_bc100,age_strat_2023$comb_bc100)
abc100_age_total <- rbind(age_strat_2016$comb_abc100, age_strat_2018$comb_abc100, age_strat_2021$comb_abc100, age_strat_2022$comb_abc100,age_strat_2023$comb_abc100)
abcn100_age_total <- rbind(age_strat_2016$comb_abcn100, age_strat_2018$comb_abcn100, age_strat_2021$comb_abcn100, age_strat_2022$comb_abcn100,age_strat_2023$comb_abcn100)

#Stratified by educational level
gly_edu_total <- rbind(edu_strat_2016$gly_control, edu_strat_2018$gly_control, edu_strat_2021$gly_control, edu_strat_2022$gly_control,edu_strat_2023$gly_control)
bp_edu_total <- rbind(edu_strat_2016$high_bp, edu_strat_2018$high_bp, edu_strat_2021$high_bp, edu_strat_2022$high_bp,edu_strat_2023$high_bp)
ldl_edu_total <- rbind(edu_strat_2016$high_ldl, edu_strat_2018$high_ldl, edu_strat_2021$high_ldl, edu_strat_2022$high_ldl,edu_strat_2023$high_ldl)
sldl_edu_total <- rbind(edu_strat_2016$high_sldl, edu_strat_2018$high_sldl, edu_strat_2021$high_sldl, edu_strat_2022$high_sldl,edu_strat_2023$high_sldl)
sldl100_edu_total <- rbind(edu_strat_2016$high_sldl100, edu_strat_2018$high_sldl100, edu_strat_2021$high_sldl100, edu_strat_2022$high_sldl100,edu_strat_2023$high_sldl100)
smoke_edu_total <- rbind(edu_strat_2016$smoke_quit, edu_strat_2018$smoke_quit, edu_strat_2021$smoke_quit, edu_strat_2022$smoke_quit,edu_strat_2023$smoke_quit)
bc_edu_total <- rbind(edu_strat_2016$comb_bc, edu_strat_2018$comb_bc, edu_strat_2021$comb_bc, edu_strat_2022$comb_bc,edu_strat_2023$comb_bc)
abc_edu_total <- rbind(edu_strat_2016$comb_abc, edu_strat_2018$comb_abc, edu_strat_2021$comb_abc, edu_strat_2022$comb_abc,edu_strat_2023$comb_abc)
abcn_edu_total <- rbind(edu_strat_2016$comb_abcn, edu_strat_2018$comb_abcn, edu_strat_2021$comb_abcn, edu_strat_2022$comb_abcn,edu_strat_2023$comb_abcn)
bc100_edu_total <- rbind(edu_strat_2016$comb_bc100, edu_strat_2018$comb_bc100, edu_strat_2021$comb_bc100, edu_strat_2022$comb_bc100,edu_strat_2023$comb_bc100)
abc100_edu_total <- rbind(edu_strat_2016$comb_abc100, edu_strat_2018$comb_abc100, edu_strat_2021$comb_abc100, edu_strat_2022$comb_abc100,edu_strat_2023$comb_abc100)
abcn100_edu_total <- rbind(edu_strat_2016$comb_abcn100, edu_strat_2018$comb_abcn100, edu_strat_2021$comb_abcn100, edu_strat_2022$comb_abcn100,edu_strat_2023$comb_abcn100)

#Stratified by area (rural or urban)
gly_urban_total <- rbind(urb_strat_2016$gly_control, urb_strat_2018$gly_control, urb_strat_2021$gly_control, urb_strat_2022$gly_control,urb_strat_2023$gly_control)
bp_urban_total <- rbind(urb_strat_2016$high_bp, urb_strat_2018$high_bp, urb_strat_2021$high_bp, urb_strat_2022$high_bp,urb_strat_2023$high_bp)
ldl_urban_total <- rbind(urb_strat_2016$high_ldl, urb_strat_2018$high_ldl, urb_strat_2021$high_ldl, urb_strat_2022$high_ldl,urb_strat_2023$high_ldl)
sldl_urban_total <- rbind(urb_strat_2016$high_sldl, urb_strat_2018$high_sldl, urb_strat_2021$high_sldl, urb_strat_2022$high_sldl,urb_strat_2023$high_sldl)
sldl100_urban_total <- rbind(urb_strat_2016$high_sldl100, urb_strat_2018$high_sldl100, urb_strat_2021$high_sldl100, urb_strat_2022$high_sldl100,urb_strat_2023$high_sldl100)
smoke_urban_total <- rbind(urb_strat_2016$smoke_quit, urb_strat_2018$smoke_quit, urb_strat_2021$smoke_quit, urb_strat_2022$smoke_quit,urb_strat_2023$smoke_quit)
bc_urban_total <- rbind(urb_strat_2016$comb_bc, urb_strat_2018$comb_bc, urb_strat_2021$comb_bc, urb_strat_2022$comb_bc,urb_strat_2023$comb_bc)
abc_urban_total <- rbind(urb_strat_2016$comb_abc, urb_strat_2018$comb_abc, urb_strat_2021$comb_abc, urb_strat_2022$comb_abc,urb_strat_2023$comb_abc)
abcn_urban_total <- rbind(urb_strat_2016$comb_abcn, urb_strat_2018$comb_abcn, urb_strat_2021$comb_abcn, urb_strat_2022$comb_abcn,urb_strat_2023$comb_abcn)
bc100_urban_total <- rbind(urb_strat_2016$comb_bc100, urb_strat_2018$comb_bc100, urb_strat_2021$comb_bc100, urb_strat_2022$comb_bc100,urb_strat_2023$comb_bc100)
abc100_urban_total <- rbind(urb_strat_2016$comb_abc100, urb_strat_2018$comb_abc100, urb_strat_2021$comb_abc100, urb_strat_2022$comb_abc100,urb_strat_2023$comb_abc100)
abcn100_urban_total <- rbind(urb_strat_2016$comb_abcn100, urb_strat_2018$comb_abcn100, urb_strat_2021$comb_abcn100, urb_strat_2022$comb_abcn100,urb_strat_2023$comb_abcn100)

#Stratified by indigenous identity
gly_indigenous_total <- rbind(ind_strat_2016$gly_control, ind_strat_2018$gly_control, ind_strat_2021$gly_control, ind_strat_2022$gly_control,ind_strat_2023$gly_control)
bp_indigenous_total <- rbind(ind_strat_2016$high_bp, ind_strat_2018$high_bp, ind_strat_2021$high_bp, ind_strat_2022$high_bp,ind_strat_2023$high_bp)
ldl_indigenous_total <- rbind(ind_strat_2016$high_ldl, ind_strat_2018$high_ldl, ind_strat_2021$high_ldl, ind_strat_2022$high_ldl,ind_strat_2023$high_ldl)
sldl_indigenous_total <- rbind(ind_strat_2016$high_sldl, ind_strat_2018$high_sldl, ind_strat_2021$high_sldl, ind_strat_2022$high_sldl,ind_strat_2023$high_sldl)
sldl100_indigenous_total <- rbind(ind_strat_2016$high_sldl100, ind_strat_2018$high_sldl100, ind_strat_2021$high_sldl100, ind_strat_2022$high_sldl100,ind_strat_2023$high_sldl100)
smoke_indigenous_total <- rbind(ind_strat_2016$smoke_quit, ind_strat_2018$smoke_quit, ind_strat_2021$smoke_quit, ind_strat_2022$smoke_quit,ind_strat_2023$smoke_quit)
bc_indigenous_total <- rbind(ind_strat_2016$comb_bc, ind_strat_2018$comb_bc, ind_strat_2021$comb_bc, ind_strat_2022$comb_bc,ind_strat_2023$comb_bc)
abc_indigenous_total <- rbind(ind_strat_2016$comb_abc, ind_strat_2018$comb_abc, ind_strat_2021$comb_abc, ind_strat_2022$comb_abc, ind_strat_2023$comb_abc)
abcn_indigenous_total <- rbind(ind_strat_2016$comb_abcn, ind_strat_2018$comb_abcn, ind_strat_2021$comb_abcn, ind_strat_2022$comb_abcn,ind_strat_2023$comb_abcn)
bc100_indigenous_total <- rbind(ind_strat_2016$comb_bc100, ind_strat_2018$comb_bc100, ind_strat_2021$comb_bc100, ind_strat_2022$comb_bc100,ind_strat_2023$comb_bc100)
abc100_indigenous_total <- rbind(ind_strat_2016$comb_abc100, ind_strat_2018$comb_abc100, ind_strat_2021$comb_abc100, ind_strat_2022$comb_abc100,ind_strat_2023$comb_abc100)
abcn100_indigenous_total <- rbind(ind_strat_2016$comb_abcn100, ind_strat_2018$comb_abcn100, ind_strat_2021$comb_abcn100, ind_strat_2022$comb_abcn100,ind_strat_2023$comb_abcn100)


#Stratified by time since diagnosis
gly_chronic_total <- rbind(chr_strat_2016$gly_control, chr_strat_2018$gly_control, chr_strat_2021$gly_control, chr_strat_2022$gly_control,chr_strat_2023$gly_control)
bp_chronic_total <- rbind(chr_strat_2016$high_bp, chr_strat_2018$high_bp, chr_strat_2021$high_bp, chr_strat_2022$high_bp,chr_strat_2023$high_bp)
ldl_chronic_total <- rbind(chr_strat_2016$high_ldl, chr_strat_2018$high_ldl, chr_strat_2021$high_ldl, chr_strat_2022$high_ldl,chr_strat_2023$high_ldl)
sldl_chronic_total <- rbind(chr_strat_2016$high_sldl, chr_strat_2018$high_sldl, chr_strat_2021$high_sldl, chr_strat_2022$high_sldl,chr_strat_2023$high_sldl)
sldl100_chronic_total <- rbind(chr_strat_2016$high_sldl100, chr_strat_2018$high_sldl100, chr_strat_2021$high_sldl100, chr_strat_2022$high_sldl100,chr_strat_2023$high_sldl100)
smoke_chronic_total <- rbind(chr_strat_2016$smoke_quit, chr_strat_2018$smoke_quit, chr_strat_2021$smoke_quit, chr_strat_2022$smoke_quit,chr_strat_2023$smoke_quit)
bc_chronic_total <- rbind(chr_strat_2016$comb_bc, chr_strat_2018$comb_bc, chr_strat_2021$comb_bc, chr_strat_2022$comb_bc,chr_strat_2023$comb_bc)
abc_chronic_total <- rbind(chr_strat_2016$comb_abc, chr_strat_2018$comb_abc, chr_strat_2021$comb_abc, chr_strat_2022$comb_abc,chr_strat_2023$comb_abc)
abcn_chronic_total <- rbind(chr_strat_2016$comb_abcn, chr_strat_2018$comb_abcn, chr_strat_2021$comb_abcn, chr_strat_2022$comb_abcn,chr_strat_2023$comb_abcn)
bc100_chronic_total <- rbind(chr_strat_2016$comb_bc100, chr_strat_2018$comb_bc100, chr_strat_2021$comb_bc100, chr_strat_2022$comb_bc100,chr_strat_2023$comb_bc100)
abc100_chronic_total <- rbind(chr_strat_2016$comb_abc100, chr_strat_2018$comb_abc100, chr_strat_2021$comb_abc100, chr_strat_2022$comb_abc100,chr_strat_2023$comb_abc100)
abcn100_chronic_total <- rbind(chr_strat_2016$comb_abcn100, chr_strat_2018$comb_abcn100, chr_strat_2021$comb_abcn100, chr_strat_2022$comb_abcn100,chr_strat_2023$comb_abcn100)

#### Prevalence of CVD and CV risk categories ####

diab_2016_cvd<- map("prev_cv", wt_prop, ensanut2016_diag) %>% 
  set_names("prev_cv") %>% as.data.frame()
diab_2018_cvd <- map("prev_cv", wt_prop, ensanut2018_diag) %>% 
  set_names("prev_cv") %>% as.data.frame()
diab_2021_cvd <- map("prev_cv", wt_prop, ensanut2021_diag) %>% 
  set_names("prev_cv") %>% as.data.frame()
diab_2022_cvd <- map("prev_cv", wt_prop, ensanut2022_diag) %>% 
  set_names("prev_cv") %>% as.data.frame()
diab_2023_cvd <- map("prev_cv", wt_prop, ensanut2023_diag) %>% 
  set_names("prev_cv") %>% as.data.frame()

cvd2 <- rbind(diab_2016_cvd %>% mutate(x=2016), diab_2018_cvd %>% mutate(x=2018),
               diab_2021_cvd %>% mutate(x=2021),diab_2022_cvd %>% mutate(x=2022), diab_2023_cvd %>% mutate(x=2023)) %>%
  rename(prop=prev_cv.prop, lower=prev_cv.X2.5.., upper=prev_cv.X97.5.., Year=x) %>%
  remove_rownames() %>% mutate(Year=factor(Year))


svyciprop(~high_risk, ensanut2023_score)

diab_2016_score<- map("SCORE2_diabetes_cat2", wt_prop, ensanut2016_score) %>% 
  set_names("SCORE2_diabetes_cat2") %>% as.data.frame()
diab_2018_score <- map("SCORE2_diabetes_cat2", wt_prop, ensanut2018_score) %>% 
  set_names("SCORE2_diabetes_cat2") %>% as.data.frame()
diab_2021_score <- map("SCORE2_diabetes_cat2", wt_prop, ensanut2021_score) %>% 
  set_names("SCORE2_diabetes_cat2") %>% as.data.frame()
diab_2022_score <- map("SCORE2_diabetes_cat2", wt_prop, ensanut2022_score) %>% 
  set_names("SCORE2_diabetes_cat2") %>% as.data.frame()
diab_2023_score <- map("SCORE2_diabetes_cat2", wt_prop, ensanut2023_score) %>% 
  set_names("SCORE2_diabetes_cat2") %>% as.data.frame()

risk2<-c("Low", "Moderate", "High", "Very high")
risk3<-c("Very high", "High", "Moderate", "Low")
prev2 <- rbind(diab_2016_score %>% mutate(x=2016), diab_2018_score %>% mutate(x=2018),
               diab_2021_score %>% mutate(x=2021),diab_2022_score %>% mutate(x=2022), diab_2023_score %>% mutate(x=2023)) %>%
  mutate(risk=rep(c("High", "Low", "Moderate", "Very high"),5)) %>%
  mutate(risk=factor(risk, levels=risk2, ordered = TRUE)) %>%
  rename(prop=SCORE2_diabetes_cat2.prop, lower=SCORE2_diabetes_cat2.X2.5.., upper=SCORE2_diabetes_cat2.X97.5.., Year=x) %>%
  arrange(factor(risk, levels = risk2)) %>% remove_rownames() %>% mutate(Year=factor(Year))

#### Hypertension prevalence among individuals with hypertension####
bp_2016 <- wt_prop("high_bp", ensanut2016_bp)
bp_2018 <- wt_prop("high_bp", ensanut2018_bp)
bp_2021 <- wt_prop("high_bp", ensanut2021_bp)
bp_2022 <- wt_prop("high_bp", ensanut2022_bp)
bp_2023 <- wt_prop("high_bp", ensanut2023_bp)

bp_total <- rbind(bp_2016, bp_2018, bp_2021, bp_2022,bp_2023)
bp_total <- bp_total %>% 
  mutate(type = rep(c("Controlled", "Uncontrolled"), 5),
         year = c(2016, 2016, 2018, 2018, 2021, 2021, 2022, 2022,2023,2023)) %>% 
  rename(lower = X2.5..,
         upper = X97.5..) %>% 
  filter(type == "Controlled")


bp_sex_2016 <- wt_prop_by("high_bp", "male1", ensanut2016_bp, "high", 2016)
bp_sex_2018 <- wt_prop_by("high_bp", "male1", ensanut2018_bp, "high", 2018)
bp_sex_2021 <- wt_prop_by("high_bp", "male1", ensanut2021_bp, "high", 2021)
bp_sex_2022 <- wt_prop_by("high_bp", "male1", ensanut2022_bp, "high", 2022)
bp_sex_2023 <- wt_prop_by("high_bp", "male1", ensanut2023_bp, "high", 2023)

bp_sex_total <- rbind(bp_sex_2016, bp_sex_2018, bp_sex_2021, bp_sex_2022, bp_sex_2023)
bp_sex_total <- bp_sex_total %>%
  mutate(var = rep(c("Controlled", "Uncontrolled"), 10))

bp_total %>% ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = factor(year))) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Hypertension control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

#### Total estimations of control for 2023####
#Glycemic control
gly_tot <- as.data.frame(svytotal(~gly_control, ensanut2023_diag, na.rm = TRUE))
gly_tot <- cbind(gly_tot, confint(svytotal(~gly_control, ensanut2023_diag, na.rm = TRUE)))

#Blood pressure control
bp_tot <- as.data.frame(svytotal(~high_bp, ensanut2023_diag, na.rm = TRUE))
bp_tot <- cbind(bp_tot, confint(svytotal(~high_bp, ensanut2023_diag, na.rm = TRUE)))

#Cholesterol control
ldl_tot <- as.data.frame(svytotal(~high_sldl, ensanut2023_diag, na.rm = TRUE))
ldl_tot <- cbind(ldl_tot, confint(svytotal(~high_sldl, ensanut2023_diag, na.rm = TRUE)))

#Smoking control
smok_tot <- as.data.frame(svytotal(~smoke_quit, ensanut2023_diag, na.rm = TRUE))
smok_tot <- cbind(smok_tot, confint(svytotal(~smoke_quit, ensanut2023_diag, na.rm = TRUE)))

#Stratification by risk category

statin_risk<-svyby(~statin, by = ~high_risk, design = ensanut2023_diag, svyciprop, na.rm = TRUE)
diabetestx_risk<-svyby(~TX_T2D, by = ~high_risk, design = ensanut2023_diag, svymean, na.rm = TRUE)
gly_risk<-svyby(~gly_control, by = ~high_risk, design = ensanut2023_diag, svymean, na.rm = TRUE)
highbp_risk<-svyby(~bp_diag, by = ~high_risk, design = ensanut2023_antro, svymean, na.rm = TRUE)
ldl_risk<-svyby(~high_ldl, by = ~SCORE2_diabetes_cat2, design = ensanut2023_diag, svymean, na.rm = TRUE)

svyby(~high_ldl, by = ~high_risk, design = ensanut2023_diag, svyciprop, na.rm = TRUE) %>% confint()

svyciprop(~high_risk, design = ensanut2023_score, na.rm = TRUE)
svyciprop(~high_risk, design = ensanut2022_score, na.rm = TRUE)
svyciprop(~high_risk, design = ensanut2021_score, na.rm = TRUE)
svyciprop(~high_risk, design = ensanut2018_score, na.rm = TRUE)
svyciprop(~high_risk, design = ensanut2016_score, na.rm = TRUE)

svytotal(~high_risk, design = ensanut2023_score, na.rm = TRUE)
svytotal(~high_risk, design = ensanut2023_score, na.rm = TRUE) %>% confint()


svyciprop(~prev_cv, design = ensanut2023_diag, na.rm = TRUE)

svytotal(~SCORE2_diabetes_cat2, design = ensanut2023_diag, na.rm = TRUE) %>% confint()

svymean(~bp_diag, design = ensanut2023_bp3, na.rm = TRUE)

svyby(~gly_control, by = "SCORE2_diabetes_cat2", svymean,design = ensanut2023_score, na.rm = TRUE )

#BC control
bc_tot <- as.data.frame(svytotal(~comb_bc, ensanut2023_diag, na.rm = TRUE))
bc_tot <- cbind(bc_tot, confint(svytotal(~comb_bc, ensanut2023_diag, na.rm = TRUE)))

#ABC control
abc_tot <- as.data.frame(svytotal(~comb_abc, ensanut2023_diag, na.rm = TRUE))
abc_tot <- cbind(abc_tot, confint(svytotal(~comb_abc, ensanut2023_diag, na.rm = TRUE)))

#ABCN control
abcn_tot <- as.data.frame(svytotal(~comb_abcn, ensanut2023_diag, na.rm = TRUE))
abcn_tot <- cbind(abcn_tot, confint(svytotal(~comb_abcn, ensanut2023_diag, na.rm = TRUE)))

tot <- rbind(gly_tot, bp_tot, ldl_tot, smok_tot, bc_tot, abc_tot, abcn_tot)
tot <- tot %>% select(-SE)
write.table(tot, file = "~/Desktop/tab.txt", col.names = TRUE, row.names = TRUE)

#### Table 1 - Logistic regression models####
#Data set
m_2016 <- ensanut2016 %>% select(gly_control, high_bp, high_ldl, high_sldl, high_sldl100, 
                                 smoke_quit, comb_bc100, comb_abc100, comb_abcn100, 
                                 education, age_group, indigenous_fct, male1, 
                                 Year, DISLI_cat, social_sec, high_risk,swts)

m_2018 <- ensanut2018 %>% select(gly_control, high_bp, high_ldl, high_sldl, high_sldl100, 
                                 smoke_quit, comb_bc100, comb_abc100, comb_abcn100, 
                                 education, age_group, indigenous_fct, male1, 
                                 Year, DISLI_cat, social_sec, high_risk,swts)

m_2021 <- ensanut2021 %>% select(gly_control, high_bp, high_ldl, high_sldl, high_sldl100,
                                 smoke_quit, comb_bc100, comb_abc100, comb_abcn100, 
                                 education, age_group, indigenous_fct, male1, 
                                 Year, DISLI_cat, social_sec, high_risk,swts)

m_2022 <- ensanut2022 %>% select(gly_control, high_bp, high_ldl, high_sldl, high_sldl100,
                                 smoke_quit, comb_bc100, comb_abc100, comb_abcn100, 
                                 education, age_group, indigenous_fct, male1, 
                                 Year, DISLI_cat, social_sec, high_risk,swts)

m_2023 <- ensanut2023 %>% select(gly_control, high_bp, high_ldl, high_sldl, high_sldl100,
                                 smoke_quit, comb_bc100, comb_abc100, comb_abcn100, 
                                 education, age_group, indigenous_fct, male1, 
                                 Year, DISLI_cat, social_sec, high_risk,swts)

m_total <- rbind(m_2016, m_2018, m_2021, m_2022,m_2023)
m_total$high_bp <- relevel(m_total$high_bp, ref = "Uncontrolled")
m_total$high_sldl <- relevel(m_total$high_sldl, ref = "≥70 mg/dL")
m_total$high_sldl100 <- relevel(m_total$high_sldl100, ref = "≥100 mg/dL")
m_total$smoke_quit <- relevel(m_total$smoke_quit, ref = "Smoker")
m_total$Year <- factor(m_total$Year, levels = c(2016, 2018, 2021, 2022, 2023), labels = c(c("2016", "2018", "2021", "2022","2023")))
m_total$high_risk<-factor(m_total$high_risk, labels = c("Low/Moderate risk", "High/Very high risk"))
m_total$high_risk<-relevel(m_total$high_risk, ref = c("High/Very high risk"))

#Models
m1 <- glm(gly_control ~ male1 + age_group + education + indigenous_fct + Year + DISLI_cat + social_sec, family = binomial, data = m_total, weights = swts)
gly_model <- tidy(m1, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(across(where(is.numeric), round, digits = 2)) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(variable = "Glycemic control")

m2 <- glm(high_bp ~ male1 + age_group + education + indigenous_fct + Year + DISLI_cat + social_sec, family = binomial, data = m_total, weights = swts)
bp_model <- tidy(m2, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(across(where(is.numeric), round, digits = 2)) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(variable = "BP control")

m3 <- glm(high_sldl100 ~ male1 + age_group + education + indigenous_fct + Year + DISLI_cat + social_sec, family = binomial, data = m_total, weights = swts)
ldl_model <- tidy(m3, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(across(where(is.numeric), round, digits = 2)) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(variable = "LDL control")

m4 <- glm(smoke_quit ~ male1 + age_group + education + indigenous_fct + Year + DISLI_cat + social_sec, family = binomial, data = m_total, weights = swts)
smoke_model <- tidy(m4, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(across(where(is.numeric), round, digits = 2)) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(variable = "Smoke control")

m5 <- glm(comb_bc100 ~ male1 + age_group + education + indigenous_fct + Year + DISLI_cat + social_sec, family = binomial, data = m_total, weights = swts)
bc_model <- tidy(m5, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(across(where(is.numeric), round, digits = 2)) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(variable = "BC control")

m6 <- glm(comb_abcn100 ~ male1 + age_group + education + indigenous_fct + Year + DISLI_cat + social_sec, family = binomial, data = m_total, weights = swts)
abcn_model <- tidy(m6, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(across(where(is.numeric), round, digits = 2)) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(variable = "ABCN control")

m7 <- glm(high_risk ~ male1 + age_group + education + indigenous_fct + Year + DISLI_cat + social_sec, family = binomial, data = m_total, weights = swts)
abcn_model <- tidy(m6, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(across(where(is.numeric), round, digits = 2)) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(variable = "ABCN control")


tab1<-m1 %>% tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)
tab2<-m2 %>% tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)
tab3<-m3 %>% tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)
tab4<-m4 %>% tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)
tab5<-m5 %>% tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)
tab6<-m6 %>% tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)
tab7<-m7 %>% tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)

table1 <- tbl_merge(list(tab1, tab2, tab3, tab4, tab5, tab6,tab7), tab_spanner = c("Glycemic control", "BP control", "Cholesterol control", "Smoking control", "Low/Moderate CVD risk","BC control", "ABCN control"))%>% 
  as_flex_table() %>% align(align = "center",part = "all") %>% autofit()

#Save table
doc <- read_docx() %>% body_add_flextable(value = table1, split = TRUE) %>%  body_end_section_landscape() %>%
  print(target = "Tables/Tabla 1.docx"); remove(tab1, tab2, tab3, tab4, tab5, tab6)

#### Figure 1 - Prevalence of control in individuals with diagnosed diabetes####
fig_gen_1 <- diag_glyc_total %>% filter(var == "Controlled") %>% 
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Glycemic control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_gen_2 <- bp_total %>% ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = factor(year))) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Hypertension control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_gen_3 <- diag_high_ldl_total %>% filter(var == "coNTROLLED") %>% 
  mutate(lower = if_else(lower < 0, 0, lower)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "LDL-C control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_gen_4 <- diag_shigh_ldl100_total %>% filter(var == "<100 mg/dL")  %>% 
  mutate(lower = if_else(lower < 0, 0, lower)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<100 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_gen_5 <- diag_smoking %>% filter(var == "Non-smoker") %>% 
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Non-current smoking",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_1 <- ggarrange(ggarrange(fig_gen_1, fig_gen_2, fig_gen_5, nrow = 1, ncol=3, labels="AUTO", common.legend = T,legend = "bottom"))

ggsave(fig_1, file="Figures/Figure 1.jpg", bg="transparent",
       width = 25, height = 10, units = c("cm"), dpi = 600, limitsize = FALSE)

#### Figure 2 - CV risk categories and management ####

leg<-prev2 %>% 
  ggplot(aes(x=Year, y=prop, fill=risk), xLabels=NA) +
  geom_bar(stat="identity", color="black", linetype=1) + labs(fill="") +
  ggpubr::theme_pubclean() +
  scale_fill_manual(values =c("#51B364","#F0BD27","red","darkred"))


f1a<-ggarrange(prev2 %>% group_by(Year) %>%
                 mutate(risk=factor(risk, levels=risk3, ordered = TRUE)) %>%
                 mutate(pos = cumsum(prop), IC95=upper-lower,
                        upper1=pos+IC95, lower1=pos-IC95) %>% ungroup() %>% mutate(
                          "lab1"=paste0(round(prop*100,1),"%","  \n(",round(lower*100,1),"-", round(upper*100,1),")"),
                          "lower1"=ifelse(risk=="Very high risk", NA, lower1),
                          "upper1"=ifelse(risk=="Very high risk", NA, upper1)) %>% 
                 mutate(Year=factor(Year, levels=c("2023", "2022", "2021", "2018", "2016"), ordered = T)) %>%
                 ggplot(aes(x=Year, y=prop, fill=risk), xLabels=NA) +
                 geom_bar(stat="identity", color="black", linetype=1) + labs(fill="") +
                 ggpubr::theme_pubclean() +
                 scale_x_discrete(breaks = c(2023, 2022, 2021, 2018, 2016))+
                 xlab("ENSANUT cycle") + ylab ("Weighted prevalence (%)") +
                 scale_fill_manual(values =c("darkred", "red", "#F0BD27","#51B364")) +
                 ggtitle("SCORE2-Diabetes risk categories") + theme(legend.position = "bottom") +
                 theme(plot.title = element_text(size=15, face="bold", hjust=0.5, vjust=0)) +
                 geom_text(aes(label = lab1),
                           size=2.7, col="white", position = position_stack(vjust = .5), fontface="bold.italic") +
                 scale_y_continuous(labels = scales::percent)+coord_flip()+
                 theme(legend.position = "none")) %>%
  ggarrange(legend.grob = get_legend(leg), legend = "bottom")

f1b<-svyby(~TX_T2D, ~high_risk, ensanut2023_diag, svymean, na.rm = TRUE) %>% 
  pivot_longer(-1, names_to = c(".value", "type"), names_sep = "_") %>% 
  rename(prop = TX, se = se.TX) %>% 
  mutate(lower = prop - (se * 1.96),upper = prop + (se * 1.96)) %>%
  select(high_risk, type, prop, lower, upper) %>% 
  mutate(type = factor(type, labels = c("Both", "Insulin", "None", "Pills"))) %>% 
  mutate(high_risk=factor(high_risk, labels=c("Low/Moderate", "High/Very high"))) %>%
  ggplot(aes(high_risk, prop, fill = type)) +
  geom_col(position = "stack", width = 0.65) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.01)) +
  labs(title = "Diabetes treatment",
       x = "SCORE2-Diabetes risk category",
       y = "Proportion of diabetes treatment",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = c("#3171BC", "#E8C241", "#7E7E7E","#665590")) +
  geom_text(aes(label = paste0(round(prop*100, 1), "%")),
            color = "white", 
            position = position_stack(vjust = 0.5))

f1c<-svyby(~TX_HBP, by = ~high_risk, design =  ensanut2023_diag, svymean, na.rm = TRUE)%>% 
  pivot_longer(-1, names_to = c(".value", "type"), names_sep = "_") %>% 
  rename(prop = TX, se = se.TX) %>% filter(type=="HBPTreated") %>%
  mutate(lower = prop - (se * 1.96),upper = c(1.0, 1.0)) %>%
  select(high_risk, type, prop, lower, upper) %>% 
  mutate(high_risk=factor(high_risk, labels=c("Low/Moderate", "High/Very high"))) %>%
  ggplot(aes(high_risk, prop)) +
  geom_col(position = "stack", width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)+
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.01)) +
  labs(title = "Blood pressure treatment",
       x = "SCORE2-Diabetes risk category",
       y = "Proportion of hypertension treatment",
       fill = "") +
  theme_classic() +
  scale_fill_jama() +
  geom_text(aes(label = paste0(round(prop*100, 1), "%", " \n(",paste0(round(lower*100,1)), "-",paste0(round(upper*100,1)),")")),
            color = "white", 
            position = position_stack(vjust = 0.5))+theme(legend.position = "none")

f1d<-svyby(~statin, by = ~high_risk, design = ensanut2023_diag, svymean, na.rm = TRUE)%>% 
  pivot_longer(-1, names_to = c(".value", "type"), names_sep = "_") %>% 
  select(high_risk, statinStatin, se.statinStatin) %>%
  rename(prop = statinStatin, se = se.statinStatin) %>% 
  mutate(lower = prop - (se * 1.96),upper = prop + (se * 1.96)) %>%
  mutate(upper=c(0.811, 1.00))%>%
  select(high_risk, prop, lower, upper) %>% 
  mutate(high_risk=factor(high_risk, labels=c("Low/Moderate", "High/Very high"))) %>%
  ggplot(aes(high_risk, prop)) +
  geom_col(position = "stack", width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)+
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.01)) +
  labs(title = "Statin use",
       x = "SCORE2-Diabetes risk category",
       y = "Proportion of statin use",
       fill = "") +
  theme_classic() +
  scale_fill_jama() +
  geom_text(aes(label = paste0(round(prop*100, 1), "%", " \n(",paste0(round(lower*100,1)), "-",paste0(round(upper*100,1)),")")),
            color = "white", 
            position = position_stack(vjust = 0.5))+theme(legend.position = "none")


f1e<-svyby(~bp_diag, by = ~high_risk, design = ensanut2023_diag, svymean, na.rm = TRUE)%>% 
  pivot_longer(-1, names_to = c(".value", "type"), names_sep = "_") %>% 
  rename(prop = bp, se = se.bp) %>% 
  mutate(lower = prop - (se * 1.96),upper = prop + (se * 1.96)) %>%
  select(high_risk, type, prop, lower, upper) %>% 
  mutate(high_risk=factor(high_risk, labels=c("Low/Moderate", "High/Very high"))) %>%
  mutate(type = factor(type, labels = c("Diagnosed ", "Undiagnosed","No hypertension"))) %>%
  ggplot(aes(high_risk, prop, fill = type)) +
  geom_col(position = "stack", width = 0.65) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.01)) +
  labs(title = "Hypertension diagnosis",
       x = "SCORE2-Diabetes risk category",
       y = "Proportion of hypertension diagnosis",
       fill = "") +
  theme_classic() +
  scale_fill_jama() +
  geom_text(aes(label = paste0(round(prop*100, 1), "%")),
            color = "white", 
            position = position_stack(vjust = 0.5))

f1f<-ldl_risk %>% select(high_ldlcoNTROLLED, se.high_ldlcoNTROLLED) %>%
            rename(prop = high_ldlcoNTROLLED, se = se.high_ldlcoNTROLLED) %>%
            mutate(lower = prop - (se * 1.96),upper = prop + (se * 1.96)) %>%
            mutate(risk=factor(c("High", "Low", "Moderate", "Very high"), levels=risk3)) %>%
            mutate(lower=c(0, 0.013857259, 0.020564279, 0)) %>%
            mutate("lab1"=paste0(round(prop*100,1),"%"," (",round(lower*100,1),"-", round(upper*100,1),")")) %>%
            select(risk, prop, lower, upper, lab1) %>% 
  ggplot(aes(x=risk, y=prop, fill=risk), xLabels=NA) +
  geom_bar(stat="identity", color="black", linetype=1) + labs(fill="") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)+
  theme_pubclean() +
  xlab("SCORE2-Diabetes risk category cycle") + ylab ("Weighted prevalence (%)") +
  scale_fill_manual(values =c("darkred", "red", "#F0BD27","#51B364")) +
  ggtitle("LDL-C control in 2023") + theme(legend.position = "bottom") +
    geom_text(aes(label = lab1),
            size=2.7, col="white", vjust = 1.5,
            fontface="bold.italic") +
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = "none")+
  scale_y_break(c(0.035,0.15), scales = c(1.5,5))
  

f_1a<-ggarrange(f1b, print(f1f), labels = c("B", "D"), ncol=1, nrow=2)
f_1b<-ggarrange(f1c,f1d, labels = c("C", "E"), ncol=1, nrow=2)
figure1<-ggarrange(f1a, f_1a,f_1b, labels = c("A", "", ""), nrow=1, ncol=3, widths = c(0.5, 0.3, 0.2))

ggsave(figure1, file="Figures/Figure 2.jpg", bg="transparent",
       width = 45, height = 20, units = c("cm"), dpi = 600, limitsize = FALSE)


#### Figure 3 - Prevalence of combined control####
fig_comb_1 <- diag_bc100 %>% filter(var == "Control")  %>% 
  mutate(lower = if_else(lower < 0, 0, lower)) %>% 
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.2)) +
  labs(title = "BC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_comb_2 <- diag_abc100 %>% filter(var == "Control")  %>% 
  mutate(lower = if_else(lower < 0, 0, lower)) %>% 
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.2)) +
  labs(title = "ABC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_comb_3 <- diag_abcn100 %>% filter(var == "Control")  %>% 
  mutate(lower = if_else(lower < 0, 0, lower)) %>% 
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.2)) +
  labs(title = "ABCN control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_2 <- ggarrange(fig_comb_1, fig_comb_2, fig_comb_3,
                   ncol = 3, nrow = 1, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(fig_2, file="Figures/Figure 3.jpg", bg="transparent",
       width = 30, height = 10, units = c("cm"), dpi = 600, limitsize = FALSE)

#### Figure 4 - Prevalence of treatment in individuals with diagnosed diabetes####
fig_tx_1 <- diag_tx_diab %>% 
  mutate(var = factor(var, levels = c("Both", "Insulin", "Pills", "None"))) %>% 
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = var)) +
  geom_col(position = "stack", width = 0.65) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.01)) +
  labs(title = "Diabetes treatment",
       x = "Year",
       y = "Proportion of treatment",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = c("#3171BC", "#E8C241", "#665590", "#7E7E7E")) +
  geom_text(aes(label = paste0(round(prop*100, 1), "%")),
            color = "white", 
            position = position_stack(vjust = 0.5))

fig_tx_2 <- diag_tx_bp %>% filter(var == "Treated") %>% 
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Blood pressure treatment",
       x = "Year",
       y = "Proportion of treatment",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_tx_3 <- diag_tx_ldl %>% filter(var == "Statin") %>% 
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Overall statin use",
       x = "Year",
       y = "Proportion of treatment",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_3 <- ggarrange(ggarrange(fig_tx_1, nrow=1, ncol=1, labels=c("A"), common.legend = T,legend = "bottom"), 
                   ggarrange(fig_tx_2, fig_tx_3, nrow = 1, ncol = 2, labels =c("B", "C"), common.legend = T, legend = "bottom"),
                   nrow=1, ncol=2, widths = c(0.4, 0.6))

ggsave(fig_3, file="Figures/Figure 4.jpg", bg="transparent",
       width = 30, height = 10, units = c("cm"), dpi = 600, limitsize = FALSE)

#### Supplementary Figure 1 - LDL-C and combined control using <70 mg/dL target####
fig_gen_3 <- diag_shigh_ldl_total %>% filter(var == "<70 mg/dL") %>% 
  mutate(lower = if_else(lower < 0, 0, lower)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_comb_4 <- diag_bc %>% filter(var == "Control")  %>% 
  mutate(lower = if_else(lower < 0, 0, lower)) %>%  
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.1)) +
  labs(title = "BC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_comb_5 <- diag_abc %>% filter(var == "Control")  %>% 
  mutate(lower = if_else(lower < 0, 0, lower)) %>% 
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.1)) +
  labs(title = "ABC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_comb_6 <- diag_abcn %>% filter(var == "Control")  %>% 
  mutate(lower = if_else(lower < 0, 0, lower)) %>%  
  ggplot(aes(year, prop, ymin = lower, ymax = upper, fill = year)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.1)) +
  labs(title = "ABCN control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

supp_fig_1 <- ggarrange(fig_gen_3, fig_comb_4, fig_comb_5, fig_comb_6,
                        ncol = 4, nrow = 1, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(supp_fig_1, file="Figures/supp_fig1.jpg", bg="transparent",
       width = 35, height = 10, units = c("cm"), dpi = 600, limitsize = FALSE)

#### Supplementary Figure 2 - Control prevalence stratified by sex####
fig_sex_1 <- gly_sex_total %>% filter(type == "controlControlled") %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = male1)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Glycemic control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = c("#9467bd", "#aec7e8"))

fig_sex_2 <- bp_sex_total %>% filter(type == "bpControlled") %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = male1)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Blood pressure control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = c("#9467bd", "#aec7e8"))

fig_sex_3 <- ldl_sex_total %>% filter(type == "ldl<70 mg/dL") %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = male1)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<70 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = c("#9467bd", "#aec7e8"))

fig_sex_4 <- sldl100_sex_total %>% filter(type == "sldl100<100 mg/dL") %>%
  mutate(lower = if_else(lower < 0, 0, lower)) %>%  
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = male1)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<100 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = c("#9467bd", "#aec7e8"))

fig_sex_5 <- smoke_sex_total %>% filter(type == "quitNon-smoker") %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = male1)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Non-current smoking",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = c("#9467bd", "#aec7e8"))

fig_sex_6 <- bc100_sex_total %>% filter(type == "bc100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = male1)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "BC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = c("#9467bd", "#aec7e8"))

fig_sex_7 <- abc100_sex_total %>% filter(type == "abc100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = male1)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "ABC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = c("#9467bd", "#aec7e8"))

fig_sex_8 <- abcn100_sex_total %>% filter(type == "abcn100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = male1)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "ABCN control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = c("#9467bd", "#aec7e8"))

fig_sex1 <- ggarrange(fig_sex_1, fig_sex_2, fig_sex_3, fig_sex_4, fig_sex_5, fig_sex_6, fig_sex_7, fig_sex_8, 
                      ncol=4, nrow=2, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(fig_sex1, file="Figures/supp_fig2.jpg", bg="transparent",
       width = 40, height = 20, units = c("cm"), dpi = 600, limitsize = FALSE)

#### Supplementary Figure 3 - Control prevalence stratified by age group####
fig_age_1 <- gly_age_total %>% filter(type == "controlControlled") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = age_group)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Glycemic control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_age_2 <- bp_age_total %>% filter(type == "bpControlled") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = age_group)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Blood pressure control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_age_3 <- sldl_age_total %>% filter(type == "sldl<70 mg/dL") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = age_group)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<70 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_age_4 <- sldl100_age_total %>% filter(type == "sldl100<100 mg/dL")   %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = age_group)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<100 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_age_5 <- smoke_age_total %>% filter(type == "quitNon-smoker") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = age_group)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Non-current smoking",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_age_6 <- bc100_age_total %>% filter(type == "bc100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = age_group)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "BC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_age_7 <- abc100_age_total %>% filter(type == "abc100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = age_group)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "ABC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_age_8 <- abcn100_age_total %>% filter(type == "abcn100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = age_group)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "ABCN control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_jco()

fig_age1 <- ggarrange(fig_age_1, fig_age_2, fig_age_3, fig_age_4, fig_age_5, fig_age_6, fig_age_7, fig_age_8, 
                      ncol=4, nrow=2, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(fig_age1, file="Figures/supp_fig3.jpg", bg="transparent",
       width = 40, height = 20, units = c("cm"), dpi = 600, limitsize = FALSE)

#### Supplementary Figure 4 - Control prevalence stratified by area (rural or urban)####
fig_urb_1 <- gly_urban_total %>% filter(type == "controlControlled") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = urban)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Glycemic control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_urb_2 <- bp_urban_total %>% filter(type == "bpControlled") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = urban)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Blood pressure control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_urb_3 <- sldl_urban_total %>% filter(type == "sldl<70 mg/dL") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = urban)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<70 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_urb_4 <- sldl100_urban_total %>% filter(type == "sldl100<100 mg/dL") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = urban)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<100 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_urb_5 <- smoke_urban_total %>% filter(type == "quitNon-smoker") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = urban)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Non-current smoking",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_urb_6 <- bc100_urban_total %>% filter(type == "bc100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
          upper = if_else(upper > 0.25, 0.25, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = urban)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "BC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_urb_7 <- abc100_urban_total %>% filter(type == "abc100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = urban)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "ABC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_urb_8 <- abcn100_urban_total %>% filter(type == "abcn100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = urban)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "ABCN control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_urb1 <- ggarrange(fig_urb_1, fig_urb_2, fig_urb_3, fig_urb_4, fig_urb_5, fig_urb_6, fig_urb_7, fig_urb_8, 
                      ncol=4, nrow=2, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(fig_urb1, file="Figures/supp_fig4.jpg", bg="transparent",
       width = 40, height = 20, units = c("cm"), dpi = 600, limitsize = FALSE)

#### Supplementary Figure 5 - Control prevalence stratified by by indigenous identity####
fig_ind_1 <- gly_indigenous_total %>% filter(type == "controlControlled") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = indigenous_fct)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Glycemic control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_ind_2 <- bp_indigenous_total %>% filter(type == "bpControlled") %>%
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = indigenous_fct)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Blood pressure control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_ind_3 <- sldl_indigenous_total %>% filter(type == "sldl<70 mg/dL") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = indigenous_fct)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<70 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_ind_4 <- sldl100_indigenous_total %>% filter(type == "sldl100<100 mg/dL") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = indigenous_fct)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<100 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_ind_5 <- smoke_indigenous_total %>% filter(type == "quitNon-smoker") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = indigenous_fct)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Non-current smoking",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_ind_6 <- bc100_indigenous_total %>% filter(type == "bc100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 0.3, 0.3, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = indigenous_fct)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.3)) +
  labs(title = "BC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_ind_7 <- abc100_indigenous_total %>% filter(type == "abc100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = indigenous_fct)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.3)) +
  labs(title = "ABC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_ind_8 <- abcn100_indigenous_total %>% filter(type == "abcn100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = indigenous_fct)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.3)) +
  labs(title = "ABCN control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_ind1 <- ggarrange(fig_ind_1, fig_ind_2, fig_ind_3, fig_ind_4, fig_ind_5, fig_ind_6, fig_ind_7, fig_ind_8, 
                      ncol=4, nrow=2, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(fig_ind1, file="Figures/supp_fig5.jpg", bg="transparent",
       width = 40, height = 20, units = c("cm"), dpi = 600, limitsize = FALSE)

#### Supplementary Figure 6 - Prevalence of statin use####
statin_cv_prev_2016 <- svyby(~ statin, by = ~ prevent, ensanut2016_diag, svymean, na.rm = TRUE) %>% 
  pivot_longer(-1, names_to = c(".value", "type"), names_sep = c(6,6,9,9)) %>% 
  mutate(type = rep(c("No statin", "Statin"), 3)) %>% 
  rename(prop = statin,
         se = se.sta) %>% 
  mutate(lower = prop - (se * 1.96),
         upper = prop + (se * 1.96),
         year = 2016)

statin_cv_prev_2018 <- svyby(~ statin, by = ~ prevent, ensanut2018_diag, svymean, na.rm = TRUE) %>% 
  pivot_longer(-1, names_to = c(".value", "type"), names_sep = c(6,6,9,9)) %>% 
  mutate(type = rep(c("No statin", "Statin"), 3)) %>% 
  rename(prop = statin,
         se = se.sta) %>% 
  mutate(lower = prop - (se * 1.96),
         upper = prop + (se * 1.96),
         year = 2018)

statin_cv_prev_2021 <- svyby(~ statin, by = ~ prevent, ensanut2021_diag, svymean, na.rm = TRUE) %>% 
  pivot_longer(-1, names_to = c(".value", "type"), names_sep = c(6,6,9,9)) %>% 
  mutate(type = rep(c("No statin", "Statin"), 3)) %>% 
  rename(prop = statin,
         se = se.sta) %>% 
  mutate(lower = prop - (se * 1.96),
         upper = prop + (se * 1.96),
         year = 2021)

statin_cv_prev_2022 <- svyby(~ statin, by = ~ prevent, ensanut2022_diag, svymean, na.rm = TRUE) %>% 
  pivot_longer(-1, names_to = c(".value", "type"), names_sep = c(6,6,9,9)) %>% 
  mutate(type = rep(c("No statin", "Statin"), 3)) %>% 
  rename(prop = statin,
         se = se.sta) %>% 
  mutate(lower = prop - (se * 1.96),
         upper = prop + (se * 1.96),
         year = 2022)

statin_cv_prev_2023 <- svyby(~ statin, by = ~ prevent, ensanut2023_diag, svymean, na.rm = TRUE) %>% 
  pivot_longer(-1, names_to = c(".value", "type"), names_sep = c(6,6,9,9)) %>% 
  mutate(type = rep(c("No statin", "Statin"), 3)) %>% 
  rename(prop = statin,
         se = se.sta) %>% 
  mutate(lower = prop - (se * 1.96),
         upper = prop + (se * 1.96),
         year = 2023)

statin_cv_total <- rbind(statin_cv_prev_2016, statin_cv_prev_2018, statin_cv_prev_2021, statin_cv_prev_2022,statin_cv_prev_2023)

fig_statin <- statin_cv_total %>% filter(type == "Statin") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = prevent)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_pointrange(show.legend = FALSE) +
  labs(x = "Year",
       y = "Proportion of statin use") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_wrap(~ prevent) +
  theme_classic()

ggsave(fig_statin, file="Figures/supp_fig6.jpg", bg="transparent",
       width = 18, height = 10, units = c("cm"), dpi = 600, limitsize = FALSE)


#### Other Figures - Control prevalence stratified by educational level####
fig_edu_1 <- gly_edu_total %>% filter(type == "controlControlled") %>%
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>%
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = education)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  facet_wrap(~ education, nrow = 1) +
  ylim(0, 1) +
  labs(title = "Glycemic control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_edu_2 <- bp_edu_total %>% filter(type == "bpControlled") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>%
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = education)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  facet_wrap(~ education, nrow = 1) +
  ylim(0, 1) +
  labs(title = "Blood pressure control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_edu_3 <- sldl_edu_total %>% filter(type == "sldl<70 mg/dL") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = education)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  facet_wrap(~ education, nrow = 1) +
  ylim(0, 1) +
  labs(title = "Cholesterol control (<70 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_edu_4 <- sldl100_edu_total %>% filter(type == "sldl100<100 mg/dL") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = education)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  facet_wrap(~ education, nrow = 1) +
  ylim(0, 1) +
  labs(title = "Cholesterol control (<100 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_edu_5 <- smoke_edu_total %>% filter(type == "quitNon-smoker") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = education)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  facet_wrap(~ education, nrow = 1) +
  ylim(0, 1) +
  labs(title = "Smoking control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_edu_6 <- bc_edu_total %>% filter(type == "bcControl") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = education)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  facet_wrap(~ education, nrow = 1) +
  ylim(0, 1) +
  labs(title = "BC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_edu_7 <- abc_edu_total %>% filter(type == "abcControl") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = education)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  facet_wrap(~ education, nrow = 1) +
  ylim(0, 1) +
  labs(title = "ABC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_edu_8 <- abcn_edu_total %>% filter(type == "abcnControl") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = education)) +
  geom_col(width = 0.65) +
  geom_pointrange(show.legend = FALSE) +
  facet_wrap(~ education, nrow = 1) +
  ylim(0, 1) +
  labs(title = "ABCN control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_edu <- ggarrange(fig_edu_1, fig_edu_2, fig_edu_3, fig_edu_4, fig_edu_5, nrow = 5, ncol = 1, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(fig_edu, file="~/Desktop/supp_fig_edu.jpg", bg="transparent",
       width = 50, height = 60, units = c("cm"), dpi = 600, limitsize = FALSE)

#### Other Figures - Control prevalence stratified by time since diagnosis####
fig_chr_1 <- gly_chronic_total %>% filter(type == "controlControlled") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = cronic_diag)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Glycemic control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_chr_2 <- bp_chronic_total %>% filter(type == "bpControlled") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = cronic_diag)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Blood pressure control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_chr_3 <- sldl_chronic_total %>% filter(type == "sldl<70 mg/dL") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = cronic_diag)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<70 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_chr_4 <- sldl100_chronic_total %>% filter(type == "sldl100<100 mg/dL") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = cronic_diag)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Cholesterol control (<100 mg/dL)",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_chr_5 <- smoke_chronic_total %>% filter(type == "quitNon-smoker") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = cronic_diag)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Smoking control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_chr_6 <- bc100_chronic_total %>% filter(type == "bc100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = cronic_diag)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "BC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_chr_7 <- abc100_chronic_total %>% filter(type == "abc100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = cronic_diag)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "ABC control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_chr_8 <- abcn100_chronic_total %>% filter(type == "abcn100Control") %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         upper = if_else(upper > 1, 1, upper)) %>% 
  ggplot(aes(factor(year), prop, ymin = lower, ymax = upper, fill = cronic_diag)) +
  geom_col(width = 0.65, position = position_dodge(width = 0.7)) +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.7)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  labs(title = "ABCN control",
       x = "Year",
       y = "Proportion of control",
       fill = "") +
  theme_classic() +
  scale_fill_npg()

fig_chr1 <- ggarrange(fig_chr_1, fig_chr_2, fig_chr_4, fig_chr_5, fig_chr_6, fig_chr_7, fig_chr_8, 
                      ncol=4, nrow=2, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(fig_chr1, file="~/Desktop/supp_fig5.jpg", bg="transparent",
       width = 40, height = 20, units = c("cm"), dpi = 600, limitsize = FALSE)


#### Figures - Logistic regression models####
model_total <- rbind(gly_model, bp_model, ldl_model, smoke_model, bc_model, abcn_model)
model_total$term_name <-  rep(c("Male", "45-64 years", "≥65 years", "Elementary school", "Middle/High scool", "University", "Other",
                              "Indigenous", "Year 2018", "Year 2021", "Year 2022", "DISLI low", "DISLI Moderate", "DISLI Very low", "Public insurance", "Private insurance"), 6)
model_total$term_name <- factor(model_total$term_name, levels = c("Male", "45-64 years", "≥65 years", "Elementary school", "Middle/High scool", "University", "Other",
                                                                  "Indigenous", "Year 2018", "Year 2021", "Year 2022", "DISLI low", "DISLI Moderate", "DISLI Very low","Public insurance", "Private insurance"))
model_total$variable <- factor(model_total$variable, levels = c("Glycemic control", "BP control", "LDL control", "Smoke control", "BC control", "ABCN control"))

fig_models <- model_total %>% ggplot(aes(y = fct_rev(term_name), x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_pointrange() +
  labs(x = "Estimate",
       y = "Variable") +
  scale_x_log10() +
  geom_vline(xintercept = 1, linetype = 2) +
  facet_wrap(~ variable) +
  theme_bw()

ggsave(fig_models, file="~/Proyectos/Cardiovascular risk factors in diabetes/Figuras_v3/fig_models.jpg", bg="transparent",
       width = 30, height = 18, units = c("cm"), dpi = 600, limitsize = FALSE)


#### Geographic disrtibution 2023####
#Shapes
setwd("~/Mi unidad (obello@facmed.unam.mx)/Datasets/")
geom_mx <-  sf::st_read(dsn="shapes", layer="areas_geoestadisticas_estatales")  %>% sf::st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",stringsAsFactors=FALSE)
geom_mx<-ms_simplify(geom_mx, keep = 0.01, keep_shapes = T)
geom_mx$CVE_ENT_NUM<-as.numeric(geom_mx$CVE_ENT)
geom_mx$SUBREGION[geom_mx$CVE_ENT_NUM %in% c(2,3,18,25,26)]<-1
geom_mx$SUBREGION[geom_mx$CVE_ENT_NUM %in% c(5,8,19,28)]<-2
geom_mx$SUBREGION[geom_mx$CVE_ENT_NUM %in% c(6,14,16)]<-3
geom_mx$SUBREGION[geom_mx$CVE_ENT_NUM %in% c(1,10,11,22,24,32)]<-4
geom_mx$SUBREGION[geom_mx$CVE_ENT_NUM %in% c(13,29,30)]<-5
geom_mx$SUBREGION[geom_mx$CVE_ENT_NUM %in% c(9)]<-6
geom_mx$SUBREGION[geom_mx$CVE_ENT_NUM %in% c(15)]<-7
geom_mx$SUBREGION[geom_mx$CVE_ENT_NUM %in% c(12,17,20,21)]<-8
geom_mx$SUBREGION[geom_mx$CVE_ENT_NUM %in% c(4,7,27,23,31)]<-9
geom_df.mun <-  sf::st_read(dsn="shapes", layer="areas_geoestadisticas_municipales")  %>% sf::st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",stringsAsFactors=FALSE)
setwd("~/Mi unidad (obello@facmed.unam.mx)/CV Risk Diabetes")
geom_df.mun$id<-paste0(str_pad(geom_df.mun$CVE_ENT, 2,pad = "0"),str_pad(geom_df.mun$CVE_MUN,3, pad="0"))
geo_data<-poly2nb(geom_mx)

#Maps by region
geom_gly_control <- svyby(~gly_control, by = ~SUBREGION, ensanut2023_diag, svymean, na.rm = TRUE)
geom_gly_control <- geom_gly_control %>% 
  select(1, control = 3) %>% 
  mutate(type = "Glycemic control") %>% 
  mutate(control=control*100) %>%
  left_join(geom_mx, by = "SUBREGION")

geom_bp_control <-  svyby(~high_bp, by = ~SUBREGION, ensanut2023_bp, svymean, na.rm = TRUE)
geom_bp_control <- geom_bp_control %>% 
  select(1, control = 2) %>% 
  mutate(type = "Blood pressure control") %>% 
  mutate(control=control*100) %>%
  left_join(geom_mx, by = "SUBREGION")

geom_ldl_control <- svyby(~high_sldl100, by = ~SUBREGION, ensanut2023_diag, svymean, na.rm = TRUE)
geom_ldl_control <- geom_ldl_control %>% 
  select(1, control = 2) %>% 
  mutate(type = "Cholesterol control") %>%
  mutate(control=control*100) %>%
  left_join(geom_mx, by = "SUBREGION")

geom_smoke <-  svyby(~smoke_quit, by = ~SUBREGION, ensanut2023_diag, svymean, na.rm = TRUE)
geom_smoke <- geom_smoke %>% 
  select(1, control = 2) %>% 
  mutate(type = "Non-current smoking") %>% 
  mutate(control=control*100) %>%
  left_join(geom_mx, by = "SUBREGION")

geom_bc <- svyby(~comb_bc100, by = ~SUBREGION, ensanut2023_diag, svymean, na.rm = TRUE)
geom_bc <- geom_bc %>% 
  select(1, control = 3) %>% 
  mutate(type = "BC control") %>% 
  mutate(control=control*100) %>%
  left_join(geom_mx, by = "SUBREGION")

geom_abcn <- svyby(~comb_abcn100, by = ~SUBREGION, ensanut2023_diag, svymean, na.rm = TRUE)
geom_abcn <- geom_abcn %>% 
  select(1, control = 3) %>% 
  mutate(type = "ABCN control") %>% 
  mutate(control=control*100) %>%
  left_join(geom_mx, by = "SUBREGION")

geom_risk <- svyby(~high_risk, by = ~SUBREGION, ensanut2023_diag, svymean, na.rm = TRUE)
geom_risk <- geom_risk %>% 
  select(1, control = 3) %>% 
  mutate(type = "High/Very high CVD risk") %>% 
  mutate(control=control*100) %>%
  left_join(geom_mx, by = "SUBREGION")

geom_total <- rbind(geom_gly_control, geom_bp_control, geom_ldl_control, geom_smoke, geom_bc, geom_abcn,geom_risk)
geom_total$type <- factor(geom_total$type, levels = c("Glycemic control", "Blood pressure control", "Cholesterol control", "Non-current smoking", "BC control", "ABCN control","High/Very high CVD risk"))

map1 <- geom_total %>% filter(type %in% c("Glycemic control", "Blood pressure control", "Cholesterol control")) %>% ggplot() +
  geom_sf(mapping = aes(fill = control, geometry=geometry), color = "black", size = 0.1, show.legend = TRUE) +
  scale_fill_distiller(type = "seq", palette = 3, direction = 1) +
  labs(fill = "% Control") +
  facet_wrap(~ type) +
  ggthemes::theme_map()+theme(legend.position = "bottom")

map2 <- geom_total %>% filter(type %in% c("Non-current smoking")) %>% ggplot() +
  geom_sf(mapping = aes(fill = control, geometry=geometry), color = "black", size = 0.1, show.legend = TRUE) +
  scale_fill_distiller(type = "seq", palette = 1, direction = 1) +
  labs(fill = "% Control") +
  facet_wrap(~ type) +
  ggthemes::theme_map()+theme(legend.position = "bottom")

map3 <- geom_total %>% filter(type %in% c("ABCN control")) %>% ggplot() +
  geom_sf(mapping = aes(fill = control, geometry=geometry), color = "black", size = 0.1, show.legend = TRUE) +
  scale_fill_distiller(type = "seq", palette = 4, direction = 1) +
  labs(fill = "% Control") +
  facet_wrap(~ type) +
  ggthemes::theme_map()+theme(legend.position = "bottom")

map4 <- geom_total %>% filter(type %in% c("High/Very high CVD risk")) %>% ggplot() +
  geom_sf(mapping = aes(fill = control, geometry=geometry), color = "black", size = 0.1, show.legend = TRUE) +
  scale_fill_distiller(type = "seq", palette = 7, direction = 1) +
  labs(fill = "% Prevalence") +
  facet_wrap(~ type) +
  ggthemes::theme_map()+theme(legend.position = "bottom")

map_pre <- ggarrange(map2, map3,map4, nrow = 1, ncol = 3)
map_total <- ggarrange(map1, map_pre, nrow = 2, ncol = 1)

ggsave(map_total, file="Figures/Figure 5.jpg", bg="transparent",
       width = 30, height = 18, units = c("cm"), dpi = 600, limitsize = FALSE)
