*************************************************************
*		METABOLIC SYNDROME AND PERINATAL DIET PAPER			*
*			  WRITTEN BY QIANHUI JIN						*
*				LAST EDITED JUNE 2024						*
*		J Am Heart Assoc. 2024;13:e035555. 					*
*			DOI: 10.1161/JAHA.124.035555 					*
*************************************************************


***THIS DO FILE USES DATA FROM THE NUMOM2B STUDY AND THE NUMOM2B HH1 STUDY***
***https://numom2b.org

***WE DO NOT HAVE PERMISSION TO POST THE DATA ONLINE***
***RESEARCHERS WHO WANT TO USE THE PUBLICLY AVAILABLE DATA CAN ACCESS IT ON THE NICHD DASH HERE: ***
***https://dash.nichd.nih.gov/ 

***NOTES:
*SOME VARIABLE NAMES IN THE DASH DATASET DIFFER FROM THOSE USED IN THIS CODE, AS WE ASSIGNED NEW, MORE DESCRIPTIVE VARIABLES NAMES***
*THE N OF THE PUBLICLY AVAIALABLE DASH DATASET IS SLIGHTLY LOWER THAN THE NUMOM2B STUDY INVESTIGATOR DATASET THAT WE HAVE USED FOR THIS ANALYSIS***


*LOAD IN IMPUTED DATA
clear all

*SET WORKING DIRECTORY
cd "I:\Bodnar Qianhui\HEI MetSx\Datasets"

*IMPORT THE IMPUTED DATASET
use "HEI and MetS Imputed_20qj_FINAL_mid_wc.dta", clear
set type double


*****************************************************************************
*Table 1																	*
*Characteristics of Participants and Healthy Eating Index-2020 Total Scores *
*****************************************************************************

*CREATE COLLAPSED VERSIONS OF SOME VARIABLES
gen agecat3=0
replace agecat3=1 if momage<25
replace agecat3=2 if momage>=25 & momage<=35
replace agecat3=3 if momage>35

gen povperc_v1cat2=0 if povperc_v1<0.2  
replace povperc_v1cat2=1 if povperc_v1>=0.2


*CREATING SQUARED VERSIONS OF CONTINUOUS VARIABLES TO USE BELOW WHEN CALCULATING SD'S
mi xeq: gen momage_sq=momage*momage
mi xeq: gen iliac_v5_sq=iliac_v5*iliac_v5
mi xeq: gen sbp_v5_sq=sbp_v5*sbp_v5
mi xeq: gen dbp_v5_sq=dbp_v5*dbp_v5
mi xeq: gen trig_v5_sq=trig_v5*trig_v5
mi xeq: gen hdl_v5_sq=hdl_v5*hdl_v5
mi xeq: gen glucose_v5_sq=glucose_v5*glucose_v5
mi xeq: gen povperc_v1_sq=povperc_v1*povperc_v1
mi xeq: gen heiy1_total_v1_sq=heiy1_total_v1*heiy1_total_v1


*CREATE A VARIABLE THAT WILL TAG EACH NUMOMID
*NEED THIS FOR CALCULATING N'S FOR THE TABLE
mi xeq: bysort numomid: generate tag =_n==1
mi xeq: count if tag==1

mi estimate: total tag, over(agecat3)
mi estimate, esampvaryok: total tag, over(momrace4)
mi estimate: total tag, over(povperccat2_v1)
mi estimate: total tag, over(bmiprepregcat6)
mi estimate: total tag, over(momrace4)
mi estimate: total tag, over(momeduc4)
mi estimate: total tag, over(married)
mi estimate: total tag, over(gravcat)
mi estimate: total tag, over(insurpub)
mi estimate: total tag, over(smokerpre)
mi estimate: total tag, over(pa_totmetwk3_v1)
mi estimate: total tag, over(prehtn)
mi estimate: total tag, over(prediab)
mi estimate: total tag, over(sleepsat3)
mi estimate: total tag, over(epds_dich_v1)
mi estimate: total tag, over(anx_dich)
mi estimate: total tag, over(povperc_v1cat2)
mi estimate: total tag, over(high_waist_il_v5)
mi estimate: total tag, over(high_bp_v5)
mi estimate: total tag, over(high_trig_v5)
mi estimate: total tag, over(low_hdl_v5)
mi estimate: total tag, over(high_glucose_v5)


mi estimate: proportion agecat3 
mi estimate: proportion povperccat2_v1 
mi estimate: proportion bmiprepregcat6 
mi estimate: proportion momrace4 
mi estimate: proportion momeduc4  
mi estimate: proportion married  
mi estimate: proportion insurpub 
mi estimate: proportion smokerpre  
mi estimate: proportion pa_totmetwk3_v1  
mi estimate: proportion prehtn  
mi estimate: proportion prediab  
mi estimate: proportion sleepsat3  
mi estimate: proportion epds_dich_v1  
mi estimate: proportion anx_dich  
mi estimate: proportion povperc_v1cat2 
mi estimate: proportion high_waist_il_v5  
mi estimate: proportion high_bp_v5  
mi estimate: proportion high_trig_v5  
mi estimate: proportion low_hdl_v5  
mi estimate: proportion high_glucose_v5  


*NOTE: FOLLOWING CODES MIGHT GENERATE WARNING MESSAGES BUT THEY STILL PRODUCE THE DESIRED AVERAGED OUTPUT 
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ): mean heiy1_total_v1 heiy1_total_v1_sq  if gravcat==1
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ): mean heiy1_total_v1 heiy1_total_v1_sq  if gravcat==2
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ): mean heiy1_total_v1 heiy1_total_v1_sq  if gravcat==3


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ): mean heiy1_total_v1 heiy1_total_v1_sq  if momagecat3==1
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if momagecat3==2
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if momagecat3==3

mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if bmiprepregcat6==1
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if bmiprepregcat6==2
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if bmiprepregcat6==3
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if bmiprepregcat6==4
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if bmiprepregcat6==5
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if bmiprepregcat6==6


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if momrace4==1
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if momrace4==2
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if momrace4==3
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if momrace4==4

mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if momeduc4==1
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if momeduc4==2
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if momeduc4==3
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if momeduc4==4

mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if married==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if married==1


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if insurpub==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if insurpub==1


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if smokerpre==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if smokerpre==1


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if pa_totmetwk3_v1==1
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if pa_totmetwk3_v1==2
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if pa_totmetwk3_v1==3

mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if prehtn==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if prehtn==1


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if prediab==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if prediab==1


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if sleepsat3==1
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if sleepsat3==2
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if sleepsat3==3

mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if epds_dich_v1==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if epds_dich_v1==1

mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if anx_dich==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if anx_dich==1


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if povperc_v1cat2==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if povperc_v1cat2==1


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if high_waist_il_v1==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if high_waist_il_v1==1


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if high_bp_v1==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if high_bp_v1==1


mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if high_trig_v1==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if high_trig_v1==1

mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if low_hdl_v1==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if low_hdl_v1==1

mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if high_glucose_v1==0
mi estimate (sd : sqrt( _b[heiy1_total_v1_sq] - _b[heiy1_total_v1]*_b[heiy1_total_v1] ) ), esampvaryok: mean heiy1_total_v1 heiy1_total_v1_sq  if high_glucose_v1==1

*OUTPUT FROM THE ABOVE COMMANDS WAS MANUALLY ENTERED INTO THE TABLE*
  
  
  
*************************************************************************
*TABLE 2 																*
*Characteristics of participants according to the presence of 			*
*metabolic syndrome at 2- to 7-years postpartum  						*
*************************************************************************

*NOTE: FOLLOWING CODES MIGHT GENERATE WARNING MESSAGES BUT THEY STILL PRODUCE THE DESIRED AVERAGED OUTPUT 
collect clear
mi estimate: total tag, over(metsx_v5)
mi estimate: total tag, over(metsx_v5_midwc)
mi estimate: total tag, over(povperc_v1cat2)
mi estimate: total tag, over(high_waist_il_v5)
mi estimate: total tag, over(high_bp_v5)
mi estimate: total tag, over(high_trig_v5)
mi estimate: total tag, over(low_hdl_v5)
mi estimate: total tag, over(high_glucose_v5)


bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(bmiprepregcat6)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(momrace4)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(momeduc4)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(married)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(insurpub)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(smokerpre)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(pa_totmetwk3_v1)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(prehtn)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(prediab)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(sleepsat3)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(epds_dich_v1)
bysort metsx_v5: collect col_2 = r(table)["b",.]:mi estimate: total tag, over(anx_dich)


bysort metsx_v5: collect col_1 = r(table)["b",.]: mi estimate : mean momage 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[momage_sq] - _b[momage]*_b[momage] ) ): mean momage momage_sq
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion  bmiprepregcat6
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion  momrace4
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion momeduc4 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion married 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion insurpub
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion smokerpre 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion pa_totmetwk3_v1 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion prehtn 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion prediab 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion sleepsat3 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion epds_dich_v1 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion anx_dich 
bysort metsx_v5: collect col_1 = r(table)["b",.]: mi estimate : mean povperc_v1 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[povperc_v1_sq] - _b[povperc_v1]*_b[povperc_v1] ) ): mean povperc_v1 povperc_v1_sq
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion pree_acog
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion gdm 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion ptb37 
bysort metsx_v5: collect col_1 = r(table)["b",.]: mi estimate : mean iliac_v5
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[iliac_v5_sq] - _b[iliac_v5]*_b[iliac_v5] ) ): mean iliac_v5 iliac_v5_sq 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion  high_waist_il_v5  if tag ==1
bysort metsx_v5: collect col_1 = r(table)["b",.]: mi estimate : mean sbp_v5
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[sbp_v5_sq] - _b[sbp_v5]*_b[sbp_v5] ) ): mean sbp_v5 sbp_v5_sq
bysort metsx_v5: collect col_1 = r(table)["b",.]: mi estimate : mean  dbp_v5
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[dbp_v5_sq] - _b[dbp_v5]*_b[dbp_v5] ) ): mean dbp_v5 dbp_v5_sq 
bysort metsx_v5: collect col_1 = r(table)["b",.]: mi estimate : mean trig_v5 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[trig_v5_sq] - _b[trig_v5 ]*_b[trig_v5 ] ) ): mean trig_v5  trig_v5_sq
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion  high_trig_v5 
bysort metsx_v5: collect col_1 = r(table)["b",.]: mi estimate : mean hdl_v5 
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[hdl_v5_sq] - _b[hdl_v5 ]*_b[hdl_v5 ] ) ): mean hdl_v5 hdl_v5_sq
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion  low_hdl_v5 
bysort metsx_v5: collect col_1 = r(table)["b",.]: mi estimate : mean glucose_v5
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[glucose_v5_sq] - _b[glucose_v5 ]*_b[glucose_v5 ] ) ): mean glucose_v5  glucose_v5_sq
bysort metsx_v5: collect col_2 = r(table)["b",.]: mi estimate: proportion  high_glucose_v5 

*OUTPUT FROM THE ABOVE COMMANDS WAS MANUALLY ENTERED INTO THE TABLE



********************************************************************************
*TABLE 3 Association between Healthy Eating Index-2020 total score and postpartum metabolic syndrome in the Nulliparous Pregnancy Outcomes Study
********************************************************************************
* RESULTS SHOWN IN TABLE 3 WERE GENERATED FROM REGRESSION MODELS IN DO FILE "PeriDiet_MetS_Modeling_FINAL"



	
********************************************************************************
*TABLE S2 Healthy Eating Index-2020 component scores according to postpartum metabolic syndrome in the Nulliparous Pregnancy Outcomes Study: monitoring mothers-to-be Heart Health Study (nuMoM2b-HHS), n=4423.
********************************************************************************
	  
use  "I:\Bodnar Qianhui\HEI MetSx\Datasets\HEI and MetS Imputed_20qj_FINAL.dta"
collect clear

*CREATING SQUARED VERSIONS OF CONTINUOUS VARIABLES TO USE BELOW WHEN CALCULATING SD'S
gen heiy1_totalveg_v1_sq=heiy1_totalveg_v1*heiy1_totalveg_v1
gen heiy2_green_and_bean_v1_sq=heiy2_green_and_bean_v1*heiy2_green_and_bean_v1
gen heiy3_totalfruit_v1_sq=heiy3_totalfruit_v1*heiy3_totalfruit_v1
gen heiy4_wholefruit_v1_sq=heiy4_wholefruit_v1*heiy4_wholefruit_v1
gen heiy5_wholegrain_v1_sq=heiy5_wholegrain_v1*heiy5_wholegrain_v1
gen heiy6_totaldairy_v1_sq=heiy6_totaldairy_v1*heiy6_totaldairy_v1
gen heiy7_totprot_v1_sq=heiy7_totprot_v1*heiy7_totprot_v1
gen heiy8_seaplant_prot_v1_sq=heiy8_seaplant_prot_v1*heiy8_seaplant_prot_v1
gen heiy9_fattyacid_v1_sq=heiy9_fattyacid_v1*heiy9_fattyacid_v1
gen heiy10_sodium_v1_sq=heiy10_sodium_v1*heiy10_sodium_v1
gen heiy11_refinedgrain_v1_sq=heiy11_refinedgrain_v1*heiy11_refinedgrain_v1
gen heiy12_addsug_v1_sq=heiy12_addsug_v1*heiy12_addsug_v1
gen heiy13_sfa_v1_sq=heiy13_sfa_v1*heiy13_sfa_v1


bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy3_totalfruit_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy3_totalfruit_v1_sq] - _b[heiy3_totalfruit_v1]*_b[heiy3_totalfruit_v1] ) ): mean heiy3_totalfruit_v1 heiy3_totalfruit_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy4_wholefruit_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy4_wholefruit_v1_sq] - _b[heiy4_wholefruit_v1]*_b[heiy4_wholefruit_v1] ) ): mean heiy4_wholefruit_v1 heiy4_wholefruit_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy1_totalveg_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy1_totalveg_v1_sq] - _b[heiy1_totalveg_v1]*_b[heiy1_totalveg_v1] ) ): mean heiy1_totalveg_v1 heiy1_totalveg_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy2_green_and_bean_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy2_green_and_bean_v1_sq] - _b[heiy2_green_and_bean_v1]*_b[heiy2_green_and_bean_v1] ) ): mean heiy2_green_and_bean_v1 heiy2_green_and_bean_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy5_wholegrain_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy5_wholegrain_v1_sq] - _b[heiy5_wholegrain_v1]*_b[heiy5_wholegrain_v1] ) ): mean heiy5_wholegrain_v1 heiy5_wholegrain_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy6_totaldairy_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy6_totaldairy_v1_sq] - _b[heiy6_totaldairy_v1]*_b[heiy6_totaldairy_v1] ) ): mean heiy6_totaldairy_v1 heiy6_totaldairy_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy7_totprot_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy7_totprot_v1_sq] - _b[heiy7_totprot_v1]*_b[heiy7_totprot_v1] ) ): mean heiy7_totprot_v1 heiy7_totprot_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy8_seaplant_prot_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy8_seaplant_prot_v1_sq] - _b[heiy8_seaplant_prot_v1]*_b[heiy8_seaplant_prot_v1] ) ): mean heiy8_seaplant_prot_v1 heiy8_seaplant_prot_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy9_fattyacid_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy9_fattyacid_v1_sq] - _b[heiy9_fattyacid_v1]*_b[heiy9_fattyacid_v1] ) ): mean heiy9_fattyacid_v1 heiy9_fattyacid_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy11_refinedgrain_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy11_refinedgrain_v1_sq] - _b[heiy11_refinedgrain_v1]*_b[heiy11_refinedgrain_v1] ) ): mean heiy11_refinedgrain_v1 heiy11_refinedgrain_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy10_sodium_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy10_sodium_v1_sq] - _b[heiy10_sodium_v1]*_b[heiy10_sodium_v1] ) ): mean heiy10_sodium_v1 heiy10_sodium_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy12_addsug_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy12_addsug_v1_sq] - _b[heiy12_addsug_v1]*_b[heiy12_addsug_v1] ) ): mean heiy12_addsug_v1 heiy12_addsug_v1_sq

bysort metsx_v5:collect col_1 = r(table)["b",.]: mi estimate : mean heiy13_sfa_v1 
bysort metsx_v5:collect col_2 = r(table)["b",.]: mi estimate (sd : sqrt( _b[heiy13_sfa_v1_sq] - _b[heiy13_sfa_v1]*_b[heiy13_sfa_v1] ) ): mean heiy13_sfa_v1 heiy13_sfa_v1_sq

*OUTPUT FROM THE ABOVE COMMANDS WAS MANUALLY ENTERED INTO THE TABLE



********************************************************************************
*TABLE S3 Associations between Healthy Eating Index-2020 total score and 5 health conditions that make up metabolic syndrome in the Nulliparous Pregnancy Outcomes Study: monitoring mothers-to-be Heart Health Study (nuMoM2b-HHS), n=4423.
********************************************************************************
* RESULTS SHOWN IN TABLE S3 WERE GENERATED FROM REGRESSION MODELS IN DO FILE "PeriDiet_MetS_Modeling_FINAL"


********************************************************************************
*TABLE S4 Sensitivity analysis comparing the characteristics of individuals with postpartum metabolic syndrome when defined using the iliac crest waist circumference measure versus the midpoint waist circumference measure, Nulliparous Pregnancy Outcomes Study: monitoring mothers-to-be Heart Health Study (nuMoM2b-HHS), n=4,423.
********************************************************************************

*NOTE: FOLLOWING CODES MIGHT GENERATE WARNING MESSAGES BUT THEY STILL PRODUCE THE DESIRED AVERAGED OUTPUT 

bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion   gravcat 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion  bmiprepregcat6
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion  bmiprepregcat6
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion  momrace4
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion momeduc4 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion married 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion insurpub
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion smokerpre 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion pa_totmetwk3_v1 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion prehtn 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion prediab 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion sleepsat3 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion epds_dich_v1 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion anx_dich 
bysort metsx_v5_midwc: collect col_1 = r(table)["b",.]: mi estimate: mean povperc_v1 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion pree_acog
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion gdm 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion ptb37 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion  high_waist_mid_v5 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion   high_bp_v5 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion  high_trig_v5 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion  low_hdl_v5 
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]: mi estimate: proportion  high_glucose_v5 


bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg momage, quantile(50)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg momage, quantile(25)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg momage, quantile(75)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg waistmid_v5, quantile(50)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg waistmid_v5, quantile(25)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg waistmid_v5, quantile(75)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg sbp_v5, quantile(50)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg sbp_v5, quantile(25)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg sbp_v5, quantile(75)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg dbp_v5, quantile(50)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg dbp_v5, quantile(25)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg dbp_v5, quantile(75)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg trig_v5, quantile(50)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg trig_v5, quantile(25)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg trig_v5, quantile(75)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg hdl_v5, quantile(50)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg hdl_v5, quantile(25)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg hdl_v5, quantile(75)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg glucose_v5, quantile(50)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg glucose_v5, quantile(25)
bysort metsx_v5_midwc: collect col_2 = r(table)["b",.]:mi estimate: qreg glucose_v5, quantile(75)

*OUTPUT FROM THE ABOVE COMMANDS WAS MANUALLY ENTERED INTO THE TABLE