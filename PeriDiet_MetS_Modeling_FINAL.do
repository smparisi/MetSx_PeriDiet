*************************************************************
*		METABOLIC SYNDROME AND PERINATAL DIET PAPER			*
*			WRITTEN BY SARA PARISI & QIANHUI JIN			*
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


***************************************************************************************************************************************
*THIS DO FILE CONTAINS THREE PARTS:

*PART #1 - CREATE SPLINE TERMS OF HEI-2020 TOTAL SCORE

*PART #2 - MULTIVARIABLE LOGISTIC REGRESSION MODELS FOR THE FOLLOWING:
*	-ASSOCIATION BETWEEN HEI-2020 TOTAL SCORE AND POSTPARTUM METABOLIC SYNDROME
*	-ASSOCIATION BETWEEN THE HEI-202 TOTAL SCORE AND THE METS 5 COMPONENT HEALTH CONDITIONS
*	(ELEVATED WAIST CIRCUMFERENCE, ELEVATED BLOOD PRESSURE, REDUCED HDL CHOLESTEROL, ELEVATED TRIGLYCERIDES,ELEVATED BLOOD GLUCOSE)	

*PART #3 - SENSITIVITY ANALYSIS OF USING A DIFFERENT WAIST CIRCUMFERENCE MEASURE

****************************************************************************************************************************************

clear all

*SET WORKING DIRECTORY
cd "I:\Bodnar Qianhui\HEI MetSx\Datasets"

*IMPORTING THE IMPUTED DATASET
use "HEI and MetS Imputed_20qj_FINAL_mid_wc.dta", clear
set type double

*INSTALL THE MIMRGNS COMMAND*
*Daniel Klein, 2014. "MIMRGNS: Stata module to run margins after mi estimate," Statistical Software Components S457795, Boston College Department of Economics, revised 25 Jul 2022.
*ssc install mimrgns


***************************************************************
*PART #1 - CREATE SPLINE TERMS 
***************************************************************

*CREATE THE MAIN SPLINE TERMS
gen _S_hei_total1=.
gen _S_hei_total2=.

*CREATE RANGE VARIABLE THAT THE MCP COMMAND NEEDS LATER
gen _S_hei_total1_r=.
gen _S_hei_total2_r=.

centile heiy1_total_v1, centile(1)
summarize heiy1_total_v1, detail

*REMOVE THE 1ST PERCENTILE OUTLIERS
*WE ARE ASKING IT TO MAKE EVERY 1 VALUE BETWEN 38 AND 96, 61 TERMS
mi xeq: range hei_total_range 38 96 61


*MAKE SPLINES SEPEATELY FOR EACH OF THE MI DATASETS (N=20)
forvalues k=0(1)20 {
				
		*MAKE SPLINES FROM THE HEI TOTAL SCORE VARIABLE
		mkspline _S`k'hei_total=heiy1_total_v1 if _mi_m==`k', cubic nknots(3) display
		matrix knots`k' = r(knots)
		scalar knot`k'1 = knots`k'[1,1]
		scalar knot`k'2 = knots`k'[1,2]
		scalar knot`k'3 = knots`k'[1,3]
		
		*NOW TAKE NEW RANGE VARIABLE AND CREATE A SPLINE USING THE SAME KNOTS AS ABOVE AND THE SAME SPLINE COMMAND*
		mkspline  _S`k'hei_total_r=hei_total_range if _mi_m==`k', cubic knots(`=knot`k'1' `=knot`k'2' `=knot`k'3') di
		
		*NOW COMBINING THEM ALL INTO ONE SET OF _Szscore SPLINE VARIABLES
		replace _S_hei_total1=_S`k'hei_total1  if _mi_m==`k'
		replace _S_hei_total2=_S`k'hei_total2  if _mi_m==`k'
		replace _S_hei_total1_r=_S`k'hei_total_r1  if _mi_m==`k'
		replace _S_hei_total2_r=_S`k'hei_total_r2 if _mi_m==`k'
}


*DROPPING THE NON-COMBINED SPLINE TERMS THAT WE DON'T NEED ANYMORE
drop _S0hei_total1- _S20hei_total2

*SAVE NEW TEMP DATASET THAT  WILL USE FOR MODELING 
tempfile new
save "`new'", replace



***********************************************************
*PART #2 - MULTIVARIABLE LOGISTIC REGRESSION MODELS 
***********************************************************


*****MET SX OVERALL*****

*RUNNING A LOGISTIC REGRESSION MODEL ON ONE IMPUTATION SET
mi extract 1	
		logit  metsx_v5 _S_hei_total1 _S_hei_total2 momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*USING THE MCP COMMAND TO GENERATE THE PREDICTION COMMAND WE WILL USE ON THE MI DATASET BELOW
mcp heiy1_total_v1 (_S_hei_total1- _S_hei_total2), var1(hei_total_range (_S_hei_total1_r-_S_hei_total2_r)) show ci

*NOTE: COPY AND PASTE THE LONG MARGINS COMMAND THAT DISPLAYS FROM RUNNING THIS INTO THE CODE BELOW
*THIS WILL ENSURE THAT YOU'VE GOT A SMOOTH LINE OF PREDICTIONS BELOW
*AFTER PASTING, BE SURE TO ADD "POST" TO THE MIMRGNS COMMAND IN ORDER TO BE ABLE TO SAVE THE ESTIMATES

*RUNNING A LOGISTIC REGRESSION MODEL ACROSS ALL IMPUTATIONS
clear all
use "`new'"

mi estimate, or saving("Model Estimates\miestfile_total", replace) esample(esample) dots: logit  metsx_v5 _S_hei_total1 _S_hei_total2  momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust) 


*CALCULATING PREDICTIONS AT SPECIFIC VALUES OF HEI-2020 TOTAL SCORE
*PASTE VALUES FROM THE MCP COMMAND ABOVE INTO THE CODE BELOW, ADD "POST"
mimrgns, predict(pr) post  at( _S_hei_total1=38 _S_hei_total2=0 ) at( _S_hei_total1=38.96666666666667 _S_hei_total2=0 ) at( _S_hei_total1=39.93333333333333 _S_hei_total2=0 ) at( _S_hei_total1=40.9 _S_hei_total2=0 ) at( _S_hei_total1=41.86666666666667 _S_hei_total2=0 ) at( _S_hei_total1=42.83333333333334 _S_hei_total2=0 ) at( _S_hei_total1=43.8 _S_hei_total2=0 ) at( _S_hei_total1=44.76666666666667 _S_hei_total2=0 ) at( _S_hei_total1=45.73333333333333 _S_hei_total2=0 ) at( _S_hei_total1=46.7 _S_hei_total2=0 ) at( _S_hei_total1=47.66666666666666 _S_hei_total2=0 ) at( _S_hei_total1=48.63333333333333 _S_hei_total2=0 ) at( _S_hei_total1=49.6 _S_hei_total2=.000626487988512 ) at( _S_hei_total1=50.56666666666666 _S_hei_total2=.0063462042598175 ) at( _S_hei_total1=51.53333333333333 _S_hei_total2=.023080276213218 ) at( _S_hei_total1=52.5 _S_hei_total2=.0567525803642306 ) at( _S_hei_total1=53.46666666666667 _S_hei_total2=.1132869932283724 ) at( _S_hei_total1=54.43333333333334 _S_hei_total2=.1986073913211604 ) at( _S_hei_total1=55.4 _S_hei_total2=.3186376511581107 ) at( _S_hei_total1=56.36666666666667 _S_hei_total2=.4793016492547422 ) at( _S_hei_total1=57.33333333333333 _S_hei_total2=.6865232621265694 ) at( _S_hei_total1=58.3 _S_hei_total2=.9462263662891126 ) at( _S_hei_total1=59.26666666666667 _S_hei_total2=1.264334838257887 ) at( _S_hei_total1=60.23333333333333 _S_hei_total2=1.646772554548411 ) at( _S_hei_total1=61.2 _S_hei_total2=2.0994633916762 ) at( _S_hei_total1=62.16666666666667 _S_hei_total2=2.628331226156772 ) at( _S_hei_total1=63.13333333333333 _S_hei_total2=3.239299934505639 ) at( _S_hei_total1=64.09999999999999 _S_hei_total2=3.938293393238322 ) at( _S_hei_total1=65.06666666666666 _S_hei_total2=4.731205761142141 ) at( _S_hei_total1=66.03333333333333 _S_hei_total2=5.620019450435214 ) at( _S_hei_total1=67 _S_hei_total2=6.598943550142404 ) at( _S_hei_total1=67.96666666666667 _S_hei_total2=7.661268961895122 ) at( _S_hei_total1=68.93333333333334 _S_hei_total2=8.800286587324782 ) at( _S_hei_total1=69.90000000000001 _S_hei_total2=10.00928732806279 ) at( _S_hei_total1=70.86666666666667 _S_hei_total2=11.28156208574056 ) at( _S_hei_total1=71.83333333333334 _S_hei_total2=12.6104017619895 ) at( _S_hei_total1=72.8 _S_hei_total2=13.98909725844099 ) at( _S_hei_total1=73.76666666666667 _S_hei_total2=15.4109394767265 ) at( _S_hei_total1=74.73333333333333 _S_hei_total2=16.86921931847741 ) at( _S_hei_total1=75.7 _S_hei_total2=18.35722768532513 ) at( _S_hei_total1=76.66666666666666 _S_hei_total2=19.86825547890104 ) at( _S_hei_total1=77.63333333333333 _S_hei_total2=21.39559360083662 ) at( _S_hei_total1=78.59999999999999 _S_hei_total2=22.93253295276323 ) at( _S_hei_total1=79.56666666666666 _S_hei_total2=24.47257034542139 ) at( _S_hei_total1=80.53333333333333 _S_hei_total2=26.01269731055832 ) at( _S_hei_total1=81.5 _S_hei_total2=27.55282427569524 ) at( _S_hei_total1=82.46666666666667 _S_hei_total2=29.09295124083217 ) at( _S_hei_total1=83.43333333333334 _S_hei_total2=30.63307820596908 ) at( _S_hei_total1=84.40000000000001 _S_hei_total2=32.17320517110601 ) at( _S_hei_total1=85.36666666666667 _S_hei_total2=33.71333213624293 ) at( _S_hei_total1=86.33333333333334 _S_hei_total2=35.25345910137985 ) at( _S_hei_total1=87.3 _S_hei_total2=36.79358606651675 ) at( _S_hei_total1=88.26666666666667 _S_hei_total2=38.33371303165369 ) at( _S_hei_total1=89.23333333333333 _S_hei_total2=39.8738399967906 ) at( _S_hei_total1=90.2 _S_hei_total2=41.41396696192753 ) at( _S_hei_total1=91.16666666666666 _S_hei_total2=42.95409392706442 ) at( _S_hei_total1=92.13333333333333 _S_hei_total2=44.49422089220133 ) at( _S_hei_total1=93.09999999999999 _S_hei_total2=46.03434785733828 ) at( _S_hei_total1=94.06666666666666 _S_hei_total2=47.57447482247517 ) at( _S_hei_total1=95.03333333333333 _S_hei_total2=49.11460178761212 ) at( _S_hei_total1=96 _S_hei_total2=50.65472875274904 ) 


*SAVING THE PREDICTIONS FOR LATER ANALYSIS OR GRAPHING
parmest, saving("Post Estimation Predictions\adjusted_mimrgns_HEI_total_bp80.dta",replace) 


*CALCULATING PREVALENCE DIFFERENCES AT DIFFERENT HEI-2020 TOTAL SCORES
*EACH LINE WITHIN nlcom REPRESENTS THE DIFFERENCE BETWEEN THE PREDICTED PROBABILITY AT A GIVEN SCORE AND THE REFERENCE SCORE ( SCORE 60).
nlcom ///
 ((_b[1._at] - _b[24._at])) ///
 ((_b[2._at] - _b[24._at])) ///
 ((_b[3._at] - _b[24._at])) ///
 ((_b[4._at] - _b[24._at])) ///
 ((_b[5._at] - _b[24._at])) ///
 ((_b[6._at] - _b[24._at])) ///
 ((_b[7._at] - _b[24._at])) ///
 ((_b[8._at] - _b[24._at])) ///
 ((_b[9._at] - _b[24._at])) ///
 ((_b[10._at] - _b[24._at])) ///
 ((_b[11._at] - _b[24._at])) ///
 ((_b[12._at] - _b[24._at])) ///
 ((_b[13._at] - _b[24._at])) ///
 ((_b[14._at] - _b[24._at])) ///
 ((_b[15._at] - _b[24._at])) ///
 ((_b[16._at] - _b[24._at])) ///
 ((_b[17._at] - _b[24._at])) ///
 ((_b[18._at] - _b[24._at])) ///
 ((_b[19._at] - _b[24._at])) ///
 ((_b[20._at] - _b[24._at])) ///
 ((_b[21._at] - _b[24._at])) ///
 ((_b[22._at] - _b[24._at])) ///
 ((_b[23._at] - _b[24._at])) ///
 ((_b[24._at] - _b[24._at])) ///
 ((_b[25._at] - _b[24._at])) ///
 ((_b[26._at] - _b[24._at])) ///
 ((_b[27._at] - _b[24._at])) ///
 ((_b[28._at] - _b[24._at])) ///
 ((_b[29._at] - _b[24._at])) ///
 ((_b[30._at] - _b[24._at])) ///
 ((_b[31._at] - _b[24._at])) ///
 ((_b[32._at] - _b[24._at])) ///
 ((_b[33._at] - _b[24._at])) ///
 ((_b[34._at] - _b[24._at])) ///
 ((_b[35._at] - _b[24._at])) ///
 ((_b[36._at] - _b[24._at])) ///
 ((_b[37._at] - _b[24._at])) ///
 ((_b[38._at] - _b[24._at])) ///
 ((_b[39._at] - _b[24._at])) ///
 ((_b[40._at] - _b[24._at])) ///
 ((_b[41._at] - _b[24._at])) ///
 ((_b[42._at] - _b[24._at])) ///
 ((_b[43._at] - _b[24._at])) ///
 ((_b[44._at] - _b[24._at])) ///
 ((_b[45._at] - _b[24._at])) ///
 ((_b[46._at] - _b[24._at])) ///
 ((_b[47._at] - _b[24._at])) ///
 ((_b[48._at] - _b[24._at])) ///
 ((_b[49._at] - _b[24._at])) ///
 ((_b[50._at] - _b[24._at])) ///
 ((_b[51._at] - _b[24._at])) ///
 ((_b[52._at] - _b[24._at])) ///
 ((_b[53._at] - _b[24._at])) ///
 ((_b[54._at] - _b[24._at])) ///
 ((_b[55._at] - _b[24._at])) ///
 ((_b[56._at] - _b[24._at])) ///
 ((_b[57._at] - _b[24._at])) ///
 ((_b[58._at] - _b[24._at])) ///
 ((_b[59._at] - _b[24._at])) ///
 ((_b[60._at] - _b[24._at])) ///
 ((_b[61._at] - _b[24._at])),post
 
*THE RESULTS OF THE nlcom COMMAND (I.E., THE PREVALENCE DIFFERENCES) ARE SAVED INTO A .dta FILE FOR LATER ANALYSIS OR GRAPHING.
parmest, saving("Post Estimation Predictions\RD_HEI_total_bp80.dta",replace)
 
 
 
****COMPONENT #1 - HIGH WAIST CIRCUMFERENCE****

clear all
use "`new'"

*RUNNING A LOGISTIC REGRESSION MODEL ON ONE IMPUTATION SET
mi extract 1	
		logit  high_waist_il_v5 _S_hei_total1 _S_hei_total2 momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*USING THE MCP COMMAND TO GENERATE THE PREDICTION COMMAND WE WILL USE ON THE MI DATASET BELOW	
mcp heiy1_total_v1 (_S_hei_total1- _S_hei_total2), var1(hei_total_range (_S_hei_total1_r-_S_hei_total2_r)) show ci


*NOTE: COPY AND PASTE THE LONG MARGINS COMMAND THAT DISPLAYS FROM RUNNING THIS INTO THE CODE BELOW
*THIS WILL ENSURE THAT YOU'VE GOT A SMOOTH LINE OF PREDICTIONS BELOW
*AFTER PASTING, BE SURE TO ADD "POST" TO THE MIMRGNS COMMAND IN ORDER TO BE ABLE TO SAVE THE ESTIMATES

*RUNNING A LOGISTIC REGRESSION MODEL ACROSS ALL IMPUTATIONS
clear all
use "`new'"

mi estimate, or saving("Model Estimates\WC_miestfile_total", replace) esample(esample) dots: logit  high_waist_il_v5 _S_hei_total1 _S_hei_total2  momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)


*CALCULATING PREDICTIONS AT SPECIFIC VALUES OF HEI-2020 TOTAL SCORE
*PASTE VALUES FROM THE MCP COMMAND ABOVE INTO THE CODE BELOW, ADD "POST"
mimrgns, predict( pr) post   at( _S_hei_total1=38 _S_hei_total2=0 ) at( _S_hei_total1=38.96666666666667 _S_hei_total2=0 ) at( _S_hei_total1=39.93333333333333 _S_hei_total2=0 ) at( _S_hei_total1=40.9 _S_hei_total2=0 ) at( _S_hei_total1=41.86666666666667 _S_hei_total2=0 ) at( _S_hei_total1=42.83333333333334 _S_hei_total2=0 ) at( _S_hei_total1=43.8 _S_hei_total2=0 ) at( _S_hei_total1=44.76666666666667 _S_hei_total2=0 ) at( _S_hei_total1=45.73333333333333 _S_hei_total2=0 ) at( _S_hei_total1=46.7 _S_hei_total2=0 ) at( _S_hei_total1=47.66666666666666 _S_hei_total2=0 ) at( _S_hei_total1=48.63333333333333 _S_hei_total2=0 ) at( _S_hei_total1=49.6 _S_hei_total2=.000626487988512 ) at( _S_hei_total1=50.56666666666666 _S_hei_total2=.0063462042598175 ) at( _S_hei_total1=51.53333333333333 _S_hei_total2=.023080276213218 ) at( _S_hei_total1=52.5 _S_hei_total2=.0567525803642306 ) at( _S_hei_total1=53.46666666666667 _S_hei_total2=.1132869932283724 ) at( _S_hei_total1=54.43333333333334 _S_hei_total2=.1986073913211604 ) at( _S_hei_total1=55.4 _S_hei_total2=.3186376511581107 ) at( _S_hei_total1=56.36666666666667 _S_hei_total2=.4793016492547422 ) at( _S_hei_total1=57.33333333333333 _S_hei_total2=.6865232621265694 ) at( _S_hei_total1=58.3 _S_hei_total2=.9462263662891126 ) at( _S_hei_total1=59.26666666666667 _S_hei_total2=1.264334838257887 ) at( _S_hei_total1=60.23333333333333 _S_hei_total2=1.646772554548411 ) at( _S_hei_total1=61.2 _S_hei_total2=2.0994633916762 ) at( _S_hei_total1=62.16666666666667 _S_hei_total2=2.628331226156772 ) at( _S_hei_total1=63.13333333333333 _S_hei_total2=3.239299934505639 ) at( _S_hei_total1=64.09999999999999 _S_hei_total2=3.938293393238322 ) at( _S_hei_total1=65.06666666666666 _S_hei_total2=4.731205761142141 ) at( _S_hei_total1=66.03333333333333 _S_hei_total2=5.620019450435214 ) at( _S_hei_total1=67 _S_hei_total2=6.598943550142404 ) at( _S_hei_total1=67.96666666666667 _S_hei_total2=7.661268961895122 ) at( _S_hei_total1=68.93333333333334 _S_hei_total2=8.800286587324782 ) at( _S_hei_total1=69.90000000000001 _S_hei_total2=10.00928732806279 ) at( _S_hei_total1=70.86666666666667 _S_hei_total2=11.28156208574056 ) at( _S_hei_total1=71.83333333333334 _S_hei_total2=12.6104017619895 ) at( _S_hei_total1=72.8 _S_hei_total2=13.98909725844099 ) at( _S_hei_total1=73.76666666666667 _S_hei_total2=15.4109394767265 ) at( _S_hei_total1=74.73333333333333 _S_hei_total2=16.86921931847741 ) at( _S_hei_total1=75.7 _S_hei_total2=18.35722768532513 ) at( _S_hei_total1=76.66666666666666 _S_hei_total2=19.86825547890104 ) at( _S_hei_total1=77.63333333333333 _S_hei_total2=21.39559360083662 ) at( _S_hei_total1=78.59999999999999 _S_hei_total2=22.93253295276323 ) at( _S_hei_total1=79.56666666666666 _S_hei_total2=24.47257034542139 ) at( _S_hei_total1=80.53333333333333 _S_hei_total2=26.01269731055832 ) at( _S_hei_total1=81.5 _S_hei_total2=27.55282427569524 ) at( _S_hei_total1=82.46666666666667 _S_hei_total2=29.09295124083217 ) at( _S_hei_total1=83.43333333333334 _S_hei_total2=30.63307820596908 ) at( _S_hei_total1=84.40000000000001 _S_hei_total2=32.17320517110601 ) at( _S_hei_total1=85.36666666666667 _S_hei_total2=33.71333213624293 ) at( _S_hei_total1=86.33333333333334 _S_hei_total2=35.25345910137985 ) at( _S_hei_total1=87.3 _S_hei_total2=36.79358606651675 ) at( _S_hei_total1=88.26666666666667 _S_hei_total2=38.33371303165369 ) at( _S_hei_total1=89.23333333333333 _S_hei_total2=39.8738399967906 ) at( _S_hei_total1=90.2 _S_hei_total2=41.41396696192753 ) at( _S_hei_total1=91.16666666666666 _S_hei_total2=42.95409392706442 ) at( _S_hei_total1=92.13333333333333 _S_hei_total2=44.49422089220133 ) at( _S_hei_total1=93.09999999999999 _S_hei_total2=46.03434785733828 ) at( _S_hei_total1=94.06666666666666 _S_hei_total2=47.57447482247517 ) at( _S_hei_total1=95.03333333333333 _S_hei_total2=49.11460178761212 ) at( _S_hei_total1=96 _S_hei_total2=50.65472875274904 ) 


*SAVING THE PREDICTIONS FOR LATER ANALYSIS OR GRAPHING
parmest, saving("Post Estimation Predictions\WC_adjusted_mimrgns_HEI_total.dta",replace) 


*CALCULATING PREVALENCE DIFFERENCES AT DIFFERENT HEI-2020 TOTAL SCORES
*EACH LINE WITHIN nlcom REPRESENTS THE DIFFERENCE BETWEEN THE PREDICTED PROBABILITY AT A GIVEN SCORE AND THE REFERENCE SCORE ( SCORE 60). 
nlcom ///
 ((_b[1._at] - _b[24._at])) ///
 ((_b[2._at] - _b[24._at])) ///
 ((_b[3._at] - _b[24._at])) ///
 ((_b[4._at] - _b[24._at])) ///
 ((_b[5._at] - _b[24._at])) ///
 ((_b[6._at] - _b[24._at])) ///
 ((_b[7._at] - _b[24._at])) ///
 ((_b[8._at] - _b[24._at])) ///
 ((_b[9._at] - _b[24._at])) ///
 ((_b[10._at] - _b[24._at])) ///
 ((_b[11._at] - _b[24._at])) ///
 ((_b[12._at] - _b[24._at])) ///
 ((_b[13._at] - _b[24._at])) ///
 ((_b[14._at] - _b[24._at])) ///
 ((_b[15._at] - _b[24._at])) ///
 ((_b[16._at] - _b[24._at])) ///
 ((_b[17._at] - _b[24._at])) ///
 ((_b[18._at] - _b[24._at])) ///
 ((_b[19._at] - _b[24._at])) ///
 ((_b[20._at] - _b[24._at])) ///
 ((_b[21._at] - _b[24._at])) ///
 ((_b[22._at] - _b[24._at])) ///
 ((_b[23._at] - _b[24._at])) ///
 ((_b[24._at] - _b[24._at])) ///
 ((_b[25._at] - _b[24._at])) ///
 ((_b[26._at] - _b[24._at])) ///
 ((_b[27._at] - _b[24._at])) ///
 ((_b[28._at] - _b[24._at])) ///
 ((_b[29._at] - _b[24._at])) ///
 ((_b[30._at] - _b[24._at])) ///
 ((_b[31._at] - _b[24._at])) ///
 ((_b[32._at] - _b[24._at])) ///
 ((_b[33._at] - _b[24._at])) ///
 ((_b[34._at] - _b[24._at])) ///
 ((_b[35._at] - _b[24._at])) ///
 ((_b[36._at] - _b[24._at])) ///
 ((_b[37._at] - _b[24._at])) ///
 ((_b[38._at] - _b[24._at])) ///
 ((_b[39._at] - _b[24._at])) ///
 ((_b[40._at] - _b[24._at])) ///
 ((_b[41._at] - _b[24._at])) ///
 ((_b[42._at] - _b[24._at])) ///
 ((_b[43._at] - _b[24._at])) ///
 ((_b[44._at] - _b[24._at])) ///
 ((_b[45._at] - _b[24._at])) ///
 ((_b[46._at] - _b[24._at])) ///
 ((_b[47._at] - _b[24._at])) ///
 ((_b[48._at] - _b[24._at])) ///
 ((_b[49._at] - _b[24._at])) ///
 ((_b[50._at] - _b[24._at])) ///
 ((_b[51._at] - _b[24._at])) ///
 ((_b[52._at] - _b[24._at])) ///
 ((_b[53._at] - _b[24._at])) ///
 ((_b[54._at] - _b[24._at])) ///
 ((_b[55._at] - _b[24._at])) ///
 ((_b[56._at] - _b[24._at])) ///
 ((_b[57._at] - _b[24._at])) ///
 ((_b[58._at] - _b[24._at])) ///
 ((_b[59._at] - _b[24._at])) ///
 ((_b[60._at] - _b[24._at])) ///
 ((_b[61._at] - _b[24._at])),post
 
*THE RESULTS OF THE nlcom COMMAND (I.E., THE PREVALENCE DIFFERENCES) ARE SAVED INTO A .dta FILE FOR LATER ANALYSIS OR GRAPHING.
parmest, saving("Post Estimation Predictions\WC_RD_HEI_total.dta",replace) 
 
 

 
 

****COMPONENT #2 - HIGH BLOOD PRESSURE****

clear all
use "`new'"
 
*RUNNING A LOGISTIC REGRESSION MODEL ON ONE IMPUTATION SET
mi extract 1	
		logit  high_bp_v5 _S_hei_total1 _S_hei_total2 momage  povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*USING THE MCP COMMAND TO GENERATE THE PREDICTION COMMAND WE WILL USE ON THE MI DATASET BELOW
mcp heiy1_total_v1 (_S_hei_total1- _S_hei_total2), var1(hei_total_range (_S_hei_total1_r-_S_hei_total2_r)) show ci

*NOTE: COPY AND PASTE THE LONG MARGINS COMMAND THAT DISPLAYS FROM RUNNING THIS INTO THE CODE BELOW
*THIS WILL ENSURE THAT YOU'VE GOT A SMOOTH LINE OF PREDICTIONS BELOW
*AFTER PASTING, BE SURE TO ADD "POST" TO THE MIMRGNS COMMAND IN ORDER TO BE ABLE TO SAVE THE ESTIMATES

*RUNNING A LOGISTIC REGRESSION MODEL ACROSS ALL IMPUTATIONS
clear all
use "`new'"

mi estimate, or saving("Model Estimates\BP_miestfile_total", replace) esample(esample) dots: logit  high_bp_v5 _S_hei_total1 _S_hei_total2  momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*CALCULATING PREDICTIONS AT SPECIFIC VALUES OF HEI-2020 TOTAL SCORE
*PASTE VALUES FROM THE MCP COMMAND ABOVE INTO THE CODE BELOW, ADD "POST"
mimrgns, predict( pr) post at( _S_hei_total1=38 _S_hei_total2=0 ) at( _S_hei_total1=38.96666666666667 _S_hei_total2=0 ) at( _S_hei_total1=39.93333333333333 _S_hei_total2=0 ) at( _S_hei_total1=40.9 _S_hei_total2=0 ) at( _S_hei_total1=41.86666666666667 _S_hei_total2=0 ) at( _S_hei_total1=42.83333333333334 _S_hei_total2=0 ) at( _S_hei_total1=43.8 _S_hei_total2=0 ) at( _S_hei_total1=44.76666666666667 _S_hei_total2=0 ) at( _S_hei_total1=45.73333333333333 _S_hei_total2=0 ) at( _S_hei_total1=46.7 _S_hei_total2=0 ) at( _S_hei_total1=47.66666666666666 _S_hei_total2=0 ) at( _S_hei_total1=48.63333333333333 _S_hei_total2=0 ) at( _S_hei_total1=49.6 _S_hei_total2=.000626487988512 ) at( _S_hei_total1=50.56666666666666 _S_hei_total2=.0063462042598175 ) at( _S_hei_total1=51.53333333333333 _S_hei_total2=.023080276213218 ) at( _S_hei_total1=52.5 _S_hei_total2=.0567525803642306 ) at( _S_hei_total1=53.46666666666667 _S_hei_total2=.1132869932283724 ) at( _S_hei_total1=54.43333333333334 _S_hei_total2=.1986073913211604 ) at( _S_hei_total1=55.4 _S_hei_total2=.3186376511581107 ) at( _S_hei_total1=56.36666666666667 _S_hei_total2=.4793016492547422 ) at( _S_hei_total1=57.33333333333333 _S_hei_total2=.6865232621265694 ) at( _S_hei_total1=58.3 _S_hei_total2=.9462263662891126 ) at( _S_hei_total1=59.26666666666667 _S_hei_total2=1.264334838257887 ) at( _S_hei_total1=60.23333333333333 _S_hei_total2=1.646772554548411 ) at( _S_hei_total1=61.2 _S_hei_total2=2.0994633916762 ) at( _S_hei_total1=62.16666666666667 _S_hei_total2=2.628331226156772 ) at( _S_hei_total1=63.13333333333333 _S_hei_total2=3.239299934505639 ) at( _S_hei_total1=64.09999999999999 _S_hei_total2=3.938293393238322 ) at( _S_hei_total1=65.06666666666666 _S_hei_total2=4.731205761142141 ) at( _S_hei_total1=66.03333333333333 _S_hei_total2=5.620019450435214 ) at( _S_hei_total1=67 _S_hei_total2=6.598943550142404 ) at( _S_hei_total1=67.96666666666667 _S_hei_total2=7.661268961895122 ) at( _S_hei_total1=68.93333333333334 _S_hei_total2=8.800286587324782 ) at( _S_hei_total1=69.90000000000001 _S_hei_total2=10.00928732806279 ) at( _S_hei_total1=70.86666666666667 _S_hei_total2=11.28156208574056 ) at( _S_hei_total1=71.83333333333334 _S_hei_total2=12.6104017619895 ) at( _S_hei_total1=72.8 _S_hei_total2=13.98909725844099 ) at( _S_hei_total1=73.76666666666667 _S_hei_total2=15.4109394767265 ) at( _S_hei_total1=74.73333333333333 _S_hei_total2=16.86921931847741 ) at( _S_hei_total1=75.7 _S_hei_total2=18.35722768532513 ) at( _S_hei_total1=76.66666666666666 _S_hei_total2=19.86825547890104 ) at( _S_hei_total1=77.63333333333333 _S_hei_total2=21.39559360083662 ) at( _S_hei_total1=78.59999999999999 _S_hei_total2=22.93253295276323 ) at( _S_hei_total1=79.56666666666666 _S_hei_total2=24.47257034542139 ) at( _S_hei_total1=80.53333333333333 _S_hei_total2=26.01269731055832 ) at( _S_hei_total1=81.5 _S_hei_total2=27.55282427569524 ) at( _S_hei_total1=82.46666666666667 _S_hei_total2=29.09295124083217 ) at( _S_hei_total1=83.43333333333334 _S_hei_total2=30.63307820596908 ) at( _S_hei_total1=84.40000000000001 _S_hei_total2=32.17320517110601 ) at( _S_hei_total1=85.36666666666667 _S_hei_total2=33.71333213624293 ) at( _S_hei_total1=86.33333333333334 _S_hei_total2=35.25345910137985 ) at( _S_hei_total1=87.3 _S_hei_total2=36.79358606651675 ) at( _S_hei_total1=88.26666666666667 _S_hei_total2=38.33371303165369 ) at( _S_hei_total1=89.23333333333333 _S_hei_total2=39.8738399967906 ) at( _S_hei_total1=90.2 _S_hei_total2=41.41396696192753 ) at( _S_hei_total1=91.16666666666666 _S_hei_total2=42.95409392706442 ) at( _S_hei_total1=92.13333333333333 _S_hei_total2=44.49422089220133 ) at( _S_hei_total1=93.09999999999999 _S_hei_total2=46.03434785733828 ) at( _S_hei_total1=94.06666666666666 _S_hei_total2=47.57447482247517 ) at( _S_hei_total1=95.03333333333333 _S_hei_total2=49.11460178761212 ) at( _S_hei_total1=96 _S_hei_total2=50.65472875274904 )


*SAVING THE PREDICTIONS
parmest, saving("Post Estimation Predictions\BP_adjusted_mimrgns_HEI_total.dta",replace) 


* THIS COMMAND CALCULATES THE PREVALENCE DIFFERENCES BETWEEN THE SPECIFIED PREDICTED PROBABILITIES AT DIFFERENT HEI-2020 TOTAL SCORES. EACH LINE WITHIN nlcom REPRESENTS THE DIFFERENCE BETWEEN THE PREDICTED PROBABILITY AT A GIVEN SCORE AND THE REFERENCE SCORE (SCORE 60).
nlcom ///
 ((_b[1._at] - _b[24._at])) ///
 ((_b[2._at] - _b[24._at])) ///
 ((_b[3._at] - _b[24._at])) ///
 ((_b[4._at] - _b[24._at])) ///
 ((_b[5._at] - _b[24._at])) ///
 ((_b[6._at] - _b[24._at])) ///
 ((_b[7._at] - _b[24._at])) ///
 ((_b[8._at] - _b[24._at])) ///
 ((_b[9._at] - _b[24._at])) ///
 ((_b[10._at] - _b[24._at])) ///
 ((_b[11._at] - _b[24._at])) ///
 ((_b[12._at] - _b[24._at])) ///
 ((_b[13._at] - _b[24._at])) ///
 ((_b[14._at] - _b[24._at])) ///
 ((_b[15._at] - _b[24._at])) ///
 ((_b[16._at] - _b[24._at])) ///
 ((_b[17._at] - _b[24._at])) ///
 ((_b[18._at] - _b[24._at])) ///
 ((_b[19._at] - _b[24._at])) ///
 ((_b[20._at] - _b[24._at])) ///
 ((_b[21._at] - _b[24._at])) ///
 ((_b[22._at] - _b[24._at])) ///
 ((_b[23._at] - _b[24._at])) ///
 ((_b[24._at] - _b[24._at])) ///
 ((_b[25._at] - _b[24._at])) ///
 ((_b[26._at] - _b[24._at])) ///
 ((_b[27._at] - _b[24._at])) ///
 ((_b[28._at] - _b[24._at])) ///
 ((_b[29._at] - _b[24._at])) ///
 ((_b[30._at] - _b[24._at])) ///
 ((_b[31._at] - _b[24._at])) ///
 ((_b[32._at] - _b[24._at])) ///
 ((_b[33._at] - _b[24._at])) ///
 ((_b[34._at] - _b[24._at])) ///
 ((_b[35._at] - _b[24._at])) ///
 ((_b[36._at] - _b[24._at])) ///
 ((_b[37._at] - _b[24._at])) ///
 ((_b[38._at] - _b[24._at])) ///
 ((_b[39._at] - _b[24._at])) ///
 ((_b[40._at] - _b[24._at])) ///
 ((_b[41._at] - _b[24._at])) ///
 ((_b[42._at] - _b[24._at])) ///
 ((_b[43._at] - _b[24._at])) ///
 ((_b[44._at] - _b[24._at])) ///
 ((_b[45._at] - _b[24._at])) ///
 ((_b[46._at] - _b[24._at])) ///
 ((_b[47._at] - _b[24._at])) ///
 ((_b[48._at] - _b[24._at])) ///
 ((_b[49._at] - _b[24._at])) ///
 ((_b[50._at] - _b[24._at])) ///
 ((_b[51._at] - _b[24._at])) ///
 ((_b[52._at] - _b[24._at])) ///
 ((_b[53._at] - _b[24._at])) ///
 ((_b[54._at] - _b[24._at])) ///
 ((_b[55._at] - _b[24._at])) ///
 ((_b[56._at] - _b[24._at])) ///
 ((_b[57._at] - _b[24._at])) ///
 ((_b[58._at] - _b[24._at])) ///
 ((_b[59._at] - _b[24._at])) ///
 ((_b[60._at] - _b[24._at])) ///
 ((_b[61._at] - _b[24._at])),post
 
*THE RESULTS OF THE nlcom COMMAND (I.E., THE PREVALENCE DIFFERENCES) ARE SAVED INTO A .dta FILE FOR LATER ANALYSIS OR GRAPHING.
parmest, saving("Post Estimation Predictions\BP_RD_HEI_total.dta",replace) 
 

 
 
 
 
****COMPONENT #3 - LOW HDL****

clear all
use "`new'"

*RUNNING A LOGISTIC REGRESSION MODEL ON ONE IMPUTATION SET
mi extract 1	
		logit  low_hdl_v5 _S_hei_total1 _S_hei_total2 momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*USING THE MCP COMMAND TO GENERATE THE PREDICTION COMMAND WE WILL USE ON THE MI DATASET BELOW
mcp heiy1_total_v1 (_S_hei_total1- _S_hei_total2), var1(hei_total_range (_S_hei_total1_r-_S_hei_total2_r)) show ci

*NOTE: COPY AND PASTE THE LONG MARGINS COMMAND THAT DISPLAYS FROM RUNNING THIS INTO THE CODE BELOW
*THIS WILL ENSURE THAT YOU'VE GOT A SMOOTH LINE OF PREDICTIONS BELOW
*AFTER PASTING, BE SURE TO ADD "POST" TO THE MIMRGNS COMMAND IN ORDER TO BE ABLE TO SAVE THE ESTIMATES

*RUNNING A LOGISTIC REGRESSION MODEL ACROSS ALL IMPUTATIONS
clear all
use "`new'"

mi estimate, or saving("Model Estimates\HDL_miestfile_total", replace) esample(esample) dots:  logit  low_hdl_v5 _S_hei_total1 _S_hei_total2  momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)


*CALCULATING PREDICTIONS AT SPECIFIC VALUES OF HEI-2020 TOTAL SCORE
*PASTE VALUES FROM THE MCP COMMAND ABOVE INTO THE CODE BELOW, ADD "POST"
mimrgns, predict( pr) post at( _S_hei_total1=38 _S_hei_total2=0 ) at( _S_hei_total1=38.96666666666667 _S_hei_total2=0 ) at( _S_hei_total1=39.93333333333333 _S_hei_total2=0 ) at( _S_hei_total1=40.9 _S_hei_total2=0 ) at( _S_hei_total1=41.86666666666667 _S_hei_total2=0 ) at( _S_hei_total1=42.83333333333334 _S_hei_total2=0 ) at( _S_hei_total1=43.8 _S_hei_total2=0 ) at( _S_hei_total1=44.76666666666667 _S_hei_total2=0 ) at( _S_hei_total1=45.73333333333333 _S_hei_total2=0 ) at( _S_hei_total1=46.7 _S_hei_total2=0 ) at( _S_hei_total1=47.66666666666666 _S_hei_total2=0 ) at( _S_hei_total1=48.63333333333333 _S_hei_total2=0 ) at( _S_hei_total1=49.6 _S_hei_total2=.000626487988512 ) at( _S_hei_total1=50.56666666666666 _S_hei_total2=.0063462042598175 ) at( _S_hei_total1=51.53333333333333 _S_hei_total2=.023080276213218 ) at( _S_hei_total1=52.5 _S_hei_total2=.0567525803642306 ) at( _S_hei_total1=53.46666666666667 _S_hei_total2=.1132869932283724 ) at( _S_hei_total1=54.43333333333334 _S_hei_total2=.1986073913211604 ) at( _S_hei_total1=55.4 _S_hei_total2=.3186376511581107 ) at( _S_hei_total1=56.36666666666667 _S_hei_total2=.4793016492547422 ) at( _S_hei_total1=57.33333333333333 _S_hei_total2=.6865232621265694 ) at( _S_hei_total1=58.3 _S_hei_total2=.9462263662891126 ) at( _S_hei_total1=59.26666666666667 _S_hei_total2=1.264334838257887 ) at( _S_hei_total1=60.23333333333333 _S_hei_total2=1.646772554548411 ) at( _S_hei_total1=61.2 _S_hei_total2=2.0994633916762 ) at( _S_hei_total1=62.16666666666667 _S_hei_total2=2.628331226156772 ) at( _S_hei_total1=63.13333333333333 _S_hei_total2=3.239299934505639 ) at( _S_hei_total1=64.09999999999999 _S_hei_total2=3.938293393238322 ) at( _S_hei_total1=65.06666666666666 _S_hei_total2=4.731205761142141 ) at( _S_hei_total1=66.03333333333333 _S_hei_total2=5.620019450435214 ) at( _S_hei_total1=67 _S_hei_total2=6.598943550142404 ) at( _S_hei_total1=67.96666666666667 _S_hei_total2=7.661268961895122 ) at( _S_hei_total1=68.93333333333334 _S_hei_total2=8.800286587324782 ) at( _S_hei_total1=69.90000000000001 _S_hei_total2=10.00928732806279 ) at( _S_hei_total1=70.86666666666667 _S_hei_total2=11.28156208574056 ) at( _S_hei_total1=71.83333333333334 _S_hei_total2=12.6104017619895 ) at( _S_hei_total1=72.8 _S_hei_total2=13.98909725844099 ) at( _S_hei_total1=73.76666666666667 _S_hei_total2=15.4109394767265 ) at( _S_hei_total1=74.73333333333333 _S_hei_total2=16.86921931847741 ) at( _S_hei_total1=75.7 _S_hei_total2=18.35722768532513 ) at( _S_hei_total1=76.66666666666666 _S_hei_total2=19.86825547890104 ) at( _S_hei_total1=77.63333333333333 _S_hei_total2=21.39559360083662 ) at( _S_hei_total1=78.59999999999999 _S_hei_total2=22.93253295276323 ) at( _S_hei_total1=79.56666666666666 _S_hei_total2=24.47257034542139 ) at( _S_hei_total1=80.53333333333333 _S_hei_total2=26.01269731055832 ) at( _S_hei_total1=81.5 _S_hei_total2=27.55282427569524 ) at( _S_hei_total1=82.46666666666667 _S_hei_total2=29.09295124083217 ) at( _S_hei_total1=83.43333333333334 _S_hei_total2=30.63307820596908 ) at( _S_hei_total1=84.40000000000001 _S_hei_total2=32.17320517110601 ) at( _S_hei_total1=85.36666666666667 _S_hei_total2=33.71333213624293 ) at( _S_hei_total1=86.33333333333334 _S_hei_total2=35.25345910137985 ) at( _S_hei_total1=87.3 _S_hei_total2=36.79358606651675 ) at( _S_hei_total1=88.26666666666667 _S_hei_total2=38.33371303165369 ) at( _S_hei_total1=89.23333333333333 _S_hei_total2=39.8738399967906 ) at( _S_hei_total1=90.2 _S_hei_total2=41.41396696192753 ) at( _S_hei_total1=91.16666666666666 _S_hei_total2=42.95409392706442 ) at( _S_hei_total1=92.13333333333333 _S_hei_total2=44.49422089220133 ) at( _S_hei_total1=93.09999999999999 _S_hei_total2=46.03434785733828 ) at( _S_hei_total1=94.06666666666666 _S_hei_total2=47.57447482247517 ) at( _S_hei_total1=95.03333333333333 _S_hei_total2=49.11460178761212 ) at( _S_hei_total1=96 _S_hei_total2=50.65472875274904 ) 


*SAVING THE PREDICTIONS FOR LATER ANALYSIS OR GRAPHING
parmest, saving("Post Estimation Predictions\HDL_adjusted_mimrgns_HEI_total.dta",replace) 


*CALCULATING PREVALENCE DIFFERENCES AT DIFFERENT HEI-2020 TOTAL SCORES
*EACH LINE WITHIN nlcom REPRESENTS THE DIFFERENCE BETWEEN THE PREDICTED PROBABILITY AT A GIVEN SCORE AND THE REFERENCE SCORE ( SCORE 60).
nlcom ///
 ((_b[1._at] - _b[24._at])) ///
 ((_b[2._at] - _b[24._at])) ///
 ((_b[3._at] - _b[24._at])) ///
 ((_b[4._at] - _b[24._at])) ///
 ((_b[5._at] - _b[24._at])) ///
 ((_b[6._at] - _b[24._at])) ///
 ((_b[7._at] - _b[24._at])) ///
 ((_b[8._at] - _b[24._at])) ///
 ((_b[9._at] - _b[24._at])) ///
 ((_b[10._at] - _b[24._at])) ///
 ((_b[11._at] - _b[24._at])) ///
 ((_b[12._at] - _b[24._at])) ///
 ((_b[13._at] - _b[24._at])) ///
 ((_b[14._at] - _b[24._at])) ///
 ((_b[15._at] - _b[24._at])) ///
 ((_b[16._at] - _b[24._at])) ///
 ((_b[17._at] - _b[24._at])) ///
 ((_b[18._at] - _b[24._at])) ///
 ((_b[19._at] - _b[24._at])) ///
 ((_b[20._at] - _b[24._at])) ///
 ((_b[21._at] - _b[24._at])) ///
 ((_b[22._at] - _b[24._at])) ///
 ((_b[23._at] - _b[24._at])) ///
 ((_b[24._at] - _b[24._at])) ///
 ((_b[25._at] - _b[24._at])) ///
 ((_b[26._at] - _b[24._at])) ///
 ((_b[27._at] - _b[24._at])) ///
 ((_b[28._at] - _b[24._at])) ///
 ((_b[29._at] - _b[24._at])) ///
 ((_b[30._at] - _b[24._at])) ///
 ((_b[31._at] - _b[24._at])) ///
 ((_b[32._at] - _b[24._at])) ///
 ((_b[33._at] - _b[24._at])) ///
 ((_b[34._at] - _b[24._at])) ///
 ((_b[35._at] - _b[24._at])) ///
 ((_b[36._at] - _b[24._at])) ///
 ((_b[37._at] - _b[24._at])) ///
 ((_b[38._at] - _b[24._at])) ///
 ((_b[39._at] - _b[24._at])) ///
 ((_b[40._at] - _b[24._at])) ///
 ((_b[41._at] - _b[24._at])) ///
 ((_b[42._at] - _b[24._at])) ///
 ((_b[43._at] - _b[24._at])) ///
 ((_b[44._at] - _b[24._at])) ///
 ((_b[45._at] - _b[24._at])) ///
 ((_b[46._at] - _b[24._at])) ///
 ((_b[47._at] - _b[24._at])) ///
 ((_b[48._at] - _b[24._at])) ///
 ((_b[49._at] - _b[24._at])) ///
 ((_b[50._at] - _b[24._at])) ///
 ((_b[51._at] - _b[24._at])) ///
 ((_b[52._at] - _b[24._at])) ///
 ((_b[53._at] - _b[24._at])) ///
 ((_b[54._at] - _b[24._at])) ///
 ((_b[55._at] - _b[24._at])) ///
 ((_b[56._at] - _b[24._at])) ///
 ((_b[57._at] - _b[24._at])) ///
 ((_b[58._at] - _b[24._at])) ///
 ((_b[59._at] - _b[24._at])) ///
 ((_b[60._at] - _b[24._at])) ///
 ((_b[61._at] - _b[24._at])),post
 
*THE RESULTS OF THE nlcom COMMAND (I.E., THE PREVALENCE DIFFERENCES) ARE SAVED INTO A .dta FILE FOR LATER ANALYSIS OR GRAPHING.
parmest, saving("Post Estimation Predictions\HDL_RD_HEI_total.dta",replace) 
 
 


****COMPONENT #4 - HIGH TRIGLYCERIDES****
clear all
use "`new'"

*RUNNING A LOGISTIC REGRESSION MODEL ON ONE IMPUTATION SET
mi extract 1	
		logit  high_trig_v5 _S_hei_total1 _S_hei_total2 momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*USING THE MCP COMMAND TO GENERATE THE PREDICTION COMMAND WE WILL USE ON THE MI DATASET BELOW
mcp heiy1_total_v1 (_S_hei_total1- _S_hei_total2), var1(hei_total_range (_S_hei_total1_r-_S_hei_total2_r)) show ci

		
*NOTE: COPY AND PASTE THE LONG MARGINS COMMAND THAT DISPLAYS FROM RUNNING THIS INTO THE CODE BELOW
*THIS WILL ENSURE THAT YOU'VE GOT A SMOOTH LINE OF PREDICTIONS BELOW
*AFTER PASTING, BE SURE TO ADD "POST" TO THE MIMRGNS COMMAND IN ORDER TO BE ABLE TO SAVE THE ESTIMATES

*RUNNING A LOGISTIC REGRESSION MODEL ACROSS ALL IMPUTATIONS
clear all
use "`new'"

mi estimate, or saving("Model Estimates\TG_miestfile_total", replace) esample(esample) dots: logit  high_trig_v5 _S_hei_total1 _S_hei_total2  momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*CALCULATING PREDICTIONS AT SPECIFIC VALUES OF HEI-2020 TOTAL SCORE
*PASTE VALUES FROM THE MCP COMMAND ABOVE INTO THE CODE BELOW, ADD "POST"
mimrgns, predict( pr) post  at( _S_hei_total1=38 _S_hei_total2=0 ) at( _S_hei_total1=38.96666666666667 _S_hei_total2=0 ) at( _S_hei_total1=39.93333333333333 _S_hei_total2=0  ) at( _S_hei_total1=40.9 _S_hei_total2=0 ) at( _S_hei_total1=41.86666666666667 _S_hei_total2=0 ) at( _S_hei_total1=42.83333333333334 _S_hei_total2=0 ) at( _S_hei_total1=43.8 _S_hei_total2=0 ) at( _S_hei_total1=44.76666666666667 _S_hei_total2=0 ) at( _S_hei_total1=45.73333333333333 _S_hei_total2=0 ) at( _S_hei_total1=46.7 _S_hei_total2=0 ) at( _S_hei_total1=47.66666666666666 _S_hei_total2=0 ) at( _S_hei_total1=48.63333333333333 _S_hei_total2=0 ) at( _S_hei_total1=49.6 _S_hei_total2=.000626487988512 ) at( _S_hei_total1=50.56666666666666 _S_hei_total2=.0063462042598175 ) at( _S_hei_total1=51.53333333333333 _S_hei_total2=.023080276213218 ) at( _S_hei_total1=52.5 _S_hei_total2=.0567525803642306 ) at( _S_hei_total1=53.46666666666667 _S_hei_total2=.1132869932283724 ) at( _S_hei_total1=54.43333333333334 _S_hei_total2=.1986073913211604 ) at( _S_hei_total1=55.4 _S_hei_total2=.3186376511581107 ) at( _S_hei_total1=56.36666666666667 _S_hei_total2=.4793016492547422 ) at( _S_hei_total1=57.33333333333333 _S_hei_total2=.6865232621265694 ) at( _S_hei_total1=58.3 _S_hei_total2=.9462263662891126 ) at( _S_hei_total1=59.26666666666667 _S_hei_total2=1.264334838257887 ) at( _S_hei_total1=60.23333333333333 _S_hei_total2=1.646772554548411 ) at( _S_hei_total1=61.2 _S_hei_total2=2.0994633916762 ) at( _S_hei_total1=62.16666666666667 _S_hei_total2=2.628331226156772 ) at( _S_hei_total1=63.13333333333333 _S_hei_total2=3.239299934505639 ) at( _S_hei_total1=64.09999999999999 _S_hei_total2=3.938293393238322 ) at( _S_hei_total1=65.06666666666666 _S_hei_total2=4.731205761142141 ) at( _S_hei_total1=66.03333333333333 _S_hei_total2=5.620019450435214 ) at( _S_hei_total1=67 _S_hei_total2=6.598943550142404 ) at( _S_hei_total1=67.96666666666667 _S_hei_total2=7.661268961895122 ) at( _S_hei_total1=68.93333333333334 _S_hei_total2=8.800286587324782 ) at( _S_hei_total1=69.90000000000001 _S_hei_total2=10.00928732806279 ) at( _S_hei_total1=70.86666666666667 _S_hei_total2=11.28156208574056 ) at( _S_hei_total1=71.83333333333334 _S_hei_total2=12.6104017619895 ) at( _S_hei_total1=72.8 _S_hei_total2=13.98909725844099 ) at( _S_hei_total1=73.76666666666667 _S_hei_total2=15.4109394767265 ) at( _S_hei_total1=74.73333333333333 _S_hei_total2=16.86921931847741 ) at( _S_hei_total1=75.7 _S_hei_total2=18.35722768532513 ) at( _S_hei_total1=76.66666666666666 _S_hei_total2=19.86825547890104 ) at( _S_hei_total1=77.63333333333333 _S_hei_total2=21.39559360083662 ) at( _S_hei_total1=78.59999999999999 _S_hei_total2=22.93253295276323 ) at( _S_hei_total1=79.56666666666666 _S_hei_total2=24.47257034542139 ) at( _S_hei_total1=80.53333333333333 _S_hei_total2=26.01269731055832 ) at( _S_hei_total1=81.5 _S_hei_total2=27.55282427569524 ) at( _S_hei_total1=82.46666666666667 _S_hei_total2=29.09295124083217 ) at( _S_hei_total1=83.43333333333334 _S_hei_total2=30.63307820596908 ) at( _S_hei_total1=84.40000000000001 _S_hei_total2=32.17320517110601 ) at( _S_hei_total1=85.36666666666667 _S_hei_total2=33.71333213624293 ) at( _S_hei_total1=86.33333333333334 _S_hei_total2=35.25345910137985 ) at( _S_hei_total1=87.3 _S_hei_total2=36.79358606651675 ) at( _S_hei_total1=88.26666666666667 _S_hei_total2=38.33371303165369 ) at( _S_hei_total1=89.23333333333333 _S_hei_total2=39.8738399967906 ) at( _S_hei_total1=90.2 _S_hei_total2=41.41396696192753 ) at( _S_hei_total1=91.16666666666666 _S_hei_total2=42.95409392706442 ) at( _S_hei_total1=92.13333333333333 _S_hei_total2=44.49422089220133 ) at( _S_hei_total1=93.09999999999999 _S_hei_total2=46.03434785733828 ) at( _S_hei_total1=94.06666666666666 _S_hei_total2=47.57447482247517 ) at( _S_hei_total1=95.03333333333333 _S_hei_total2=49.11460178761212 ) at( _S_hei_total1=96 _S_hei_total2=50.65472875274904 ) 


*SAVING THE PREDICTIONS FOR LATER ANALYSIS OR GRAPHING
parmest, saving("Post Estimation Predictions\TG_adjusted_mimrgns_HEI_total.dta",replace) 

*THIS COMMAND CALCULATES THE PREVALENCE DIFFERENCES BETWEEN THE SPECIFIED PREDICTED PROBABILITIES AT DIFFERENT HEI-2020 TOTAL SCORES. EACH LINE WITHIN nlcom REPRESENTS THE DIFFERENCE BETWEEN THE PREDICTED PROBABILITY AT A GIVEN SCORE AND THE REFERENCE SCORE ( SCORE 60).
nlcom ///
 ((_b[1._at] - _b[24._at])) ///
 ((_b[2._at] - _b[24._at])) ///
 ((_b[3._at] - _b[24._at])) ///
 ((_b[4._at] - _b[24._at])) ///
 ((_b[5._at] - _b[24._at])) ///
 ((_b[6._at] - _b[24._at])) ///
 ((_b[7._at] - _b[24._at])) ///
 ((_b[8._at] - _b[24._at])) ///
 ((_b[9._at] - _b[24._at])) ///
 ((_b[10._at] - _b[24._at])) ///
 ((_b[11._at] - _b[24._at])) ///
 ((_b[12._at] - _b[24._at])) ///
 ((_b[13._at] - _b[24._at])) ///
 ((_b[14._at] - _b[24._at])) ///
 ((_b[15._at] - _b[24._at])) ///
 ((_b[16._at] - _b[24._at])) ///
 ((_b[17._at] - _b[24._at])) ///
 ((_b[18._at] - _b[24._at])) ///
 ((_b[19._at] - _b[24._at])) ///
 ((_b[20._at] - _b[24._at])) ///
 ((_b[21._at] - _b[24._at])) ///
 ((_b[22._at] - _b[24._at])) ///
 ((_b[23._at] - _b[24._at])) ///
 ((_b[24._at] - _b[24._at])) ///
 ((_b[25._at] - _b[24._at])) ///
 ((_b[26._at] - _b[24._at])) ///
 ((_b[27._at] - _b[24._at])) ///
 ((_b[28._at] - _b[24._at])) ///
 ((_b[29._at] - _b[24._at])) ///
 ((_b[30._at] - _b[24._at])) ///
 ((_b[31._at] - _b[24._at])) ///
 ((_b[32._at] - _b[24._at])) ///
 ((_b[33._at] - _b[24._at])) ///
 ((_b[34._at] - _b[24._at])) ///
 ((_b[35._at] - _b[24._at])) ///
 ((_b[36._at] - _b[24._at])) ///
 ((_b[37._at] - _b[24._at])) ///
 ((_b[38._at] - _b[24._at])) ///
 ((_b[39._at] - _b[24._at])) ///
 ((_b[40._at] - _b[24._at])) ///
 ((_b[41._at] - _b[24._at])) ///
 ((_b[42._at] - _b[24._at])) ///
 ((_b[43._at] - _b[24._at])) ///
 ((_b[44._at] - _b[24._at])) ///
 ((_b[45._at] - _b[24._at])) ///
 ((_b[46._at] - _b[24._at])) ///
 ((_b[47._at] - _b[24._at])) ///
 ((_b[48._at] - _b[24._at])) ///
 ((_b[49._at] - _b[24._at])) ///
 ((_b[50._at] - _b[24._at])) ///
 ((_b[51._at] - _b[24._at])) ///
 ((_b[52._at] - _b[24._at])) ///
 ((_b[53._at] - _b[24._at])) ///
 ((_b[54._at] - _b[24._at])) ///
 ((_b[55._at] - _b[24._at])) ///
 ((_b[56._at] - _b[24._at])) ///
 ((_b[57._at] - _b[24._at])) ///
 ((_b[58._at] - _b[24._at])) ///
 ((_b[59._at] - _b[24._at])) ///
 ((_b[60._at] - _b[24._at])) ///
 ((_b[61._at] - _b[24._at])),post
 
*THE RESULTS OF THE nlcom COMMAND (I.E., THE PREVALENCE DIFFERENCES) ARE SAVED INTO A .dta FILE FOR LATER ANALYSIS OR GRAPHING.
parmest, saving("Post Estimation Predictions\TG_RD_HEI_total HERE.dta",replace) 

 
 
 
 
 
****COMPONENT #5 - HIGH GLUCOSE****

clear all
use "`new'"

*RUNNING A LOGISTIC REGRESSION MODEL ON ONE IMPUTATION SET
mi extract 1	
		logit  high_glucose_v5 _S_hei_total1 _S_hei_total2 momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*USING THE MCP COMMAND TO GENERATE THE PREDICTION COMMAND WE WILL USE ON THE MI DATASET BELOW	
mcp heiy1_total_v1 (_S_hei_total1- _S_hei_total2), var1(hei_total_range (_S_hei_total1_r-_S_hei_total2_r)) show ci

	
*NOTE: COPY AND PASTE THE LONG MARGINS COMMAND THAT DISPLAYS FROM RUNNING THIS INTO THE CODE BELOW
*THIS WILL ENSURE THAT YOU'VE GOT A SMOOTH LINE OF PREDICTIONS BELOW
*AFTER PASTING, BE SURE TO ADD "POST" TO THE MIMRGNS COMMAND IN ORDER TO BE ABLE TO SAVE THE ESTIMATES

*RUNNING A LOGISTIC REGRESSION MODEL ACROSS ALL IMPUTATIONS
clear all
use "`new'"

mi estimate, or saving("Model Estimates\GLU_miestfile_total", replace) esample(esample) dots: logit  high_glucose_v5 _S_hei_total1 _S_hei_total2  momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*CALCULATING PREDICTIONS AT SPECIFIC VALUES OF HEI-2020 TOTAL SCORE
*PASTE VALUES FROM THE MCP COMMAND ABOVE INTO THE CODE BELOW, ADD "POST"
mimrgns, predict( pr) post  at( _S_hei_total1=38 _S_hei_total2=0 ) at( _S_hei_total1=38.96666666666667 _S_hei_total2=0 ) at( _S_hei_total1=39.93333333333333 _S_hei_total2=0 ) at( _S_hei_total1=40.9 _S_hei_total2=0 ) at( _S_hei_total1=41.86666666666667 _S_hei_total2=0 ) at( _S_hei_total1=42.83333333333334 _S_hei_total2=0 ) at( _S_hei_total1=43.8 _S_hei_total2=0 ) at( _S_hei_total1=44.76666666666667 _S_hei_total2=0 ) at( _S_hei_total1=45.73333333333333 _S_hei_total2=0 ) at( _S_hei_total1=46.7 _S_hei_total2=0 ) at( _S_hei_total1=47.66666666666666 _S_hei_total2=0 ) at( _S_hei_total1=48.63333333333333 _S_hei_total2=0 ) at( _S_hei_total1=49.6 _S_hei_total2=.000626487988512 ) at( _S_hei_total1=50.56666666666666 _S_hei_total2=.0063462042598175 ) at( _S_hei_total1=51.53333333333333 _S_hei_total2=.023080276213218 ) at( _S_hei_total1=52.5 _S_hei_total2=.0567525803642306 ) at( _S_hei_total1=53.46666666666667 _S_hei_total2=.1132869932283724 ) at( _S_hei_total1=54.43333333333334 _S_hei_total2=.1986073913211604 ) at( _S_hei_total1=55.4 _S_hei_total2=.3186376511581107 ) at( _S_hei_total1=56.36666666666667 _S_hei_total2=.4793016492547422 ) at( _S_hei_total1=57.33333333333333 _S_hei_total2=.6865232621265694 ) at( _S_hei_total1=58.3 _S_hei_total2=.9462263662891126 ) at( _S_hei_total1=59.26666666666667 _S_hei_total2=1.264334838257887 ) at( _S_hei_total1=60.23333333333333 _S_hei_total2=1.646772554548411 ) at( _S_hei_total1=61.2 _S_hei_total2=2.0994633916762 ) at( _S_hei_total1=62.16666666666667 _S_hei_total2=2.628331226156772 ) at( _S_hei_total1=63.13333333333333 _S_hei_total2=3.239299934505639 ) at( _S_hei_total1=64.09999999999999 _S_hei_total2=3.938293393238322 ) at( _S_hei_total1=65.06666666666666 _S_hei_total2=4.731205761142141 ) at( _S_hei_total1=66.03333333333333 _S_hei_total2=5.620019450435214 ) at( _S_hei_total1=67 _S_hei_total2=6.598943550142404 ) at( _S_hei_total1=67.96666666666667 _S_hei_total2=7.661268961895122 ) at( _S_hei_total1=68.93333333333334 _S_hei_total2=8.800286587324782 ) at( _S_hei_total1=69.90000000000001 _S_hei_total2=10.00928732806279 ) at( _S_hei_total1=70.86666666666667 _S_hei_total2=11.28156208574056 ) at( _S_hei_total1=71.83333333333334 _S_hei_total2=12.6104017619895 ) at( _S_hei_total1=72.8 _S_hei_total2=13.98909725844099 ) at( _S_hei_total1=73.76666666666667 _S_hei_total2=15.4109394767265 ) at( _S_hei_total1=74.73333333333333 _S_hei_total2=16.86921931847741 ) at( _S_hei_total1=75.7 _S_hei_total2=18.35722768532513 ) at( _S_hei_total1=76.66666666666666 _S_hei_total2=19.86825547890104 ) at( _S_hei_total1=77.63333333333333 _S_hei_total2=21.39559360083662 ) at( _S_hei_total1=78.59999999999999 _S_hei_total2=22.93253295276323 ) at( _S_hei_total1=79.56666666666666 _S_hei_total2=24.47257034542139 ) at( _S_hei_total1=80.53333333333333 _S_hei_total2=26.01269731055832 ) at( _S_hei_total1=81.5 _S_hei_total2=27.55282427569524 ) at( _S_hei_total1=82.46666666666667 _S_hei_total2=29.09295124083217 ) at( _S_hei_total1=83.43333333333334 _S_hei_total2=30.63307820596908 ) at( _S_hei_total1=84.40000000000001 _S_hei_total2=32.17320517110601 ) at( _S_hei_total1=85.36666666666667 _S_hei_total2=33.71333213624293 ) at( _S_hei_total1=86.33333333333334 _S_hei_total2=35.25345910137985 ) at( _S_hei_total1=87.3 _S_hei_total2=36.79358606651675 ) at( _S_hei_total1=88.26666666666667 _S_hei_total2=38.33371303165369 ) at( _S_hei_total1=89.23333333333333 _S_hei_total2=39.8738399967906 ) at( _S_hei_total1=90.2 _S_hei_total2=41.41396696192753 ) at( _S_hei_total1=91.16666666666666 _S_hei_total2=42.95409392706442 ) at( _S_hei_total1=92.13333333333333 _S_hei_total2=44.49422089220133 ) at( _S_hei_total1=93.09999999999999 _S_hei_total2=46.03434785733828 ) at( _S_hei_total1=94.06666666666666 _S_hei_total2=47.57447482247517 ) at( _S_hei_total1=95.03333333333333 _S_hei_total2=49.11460178761212 ) at( _S_hei_total1=96 _S_hei_total2=50.65472875274904 )
  

*SAVING THE PREDICTIONS FOR LATER ANALYSIS OR GRAPHING
parmest, saving("Post Estimation Predictions\GLU_adjusted_mimrgns_HEI_total.dta",replace) 


*THIS COMMAND CALCULATES THE PREVALENCE DIFFERENCES BETWEEN THE SPECIFIED PREDICTED PROBABILITIES AT DIFFERENT HEI-2020 TOTAL SCORES. EACH LINE WITHIN nlcom REPRESENTS THE DIFFERENCE BETWEEN THE PREDICTED PROBABILITY AT A GIVEN SCORE AND THE REFERENCE SCORE ( SCORE 60)
nlcom ///
 ((_b[1._at] - _b[24._at])) ///
 ((_b[2._at] - _b[24._at])) ///
 ((_b[3._at] - _b[24._at])) ///
 ((_b[4._at] - _b[24._at])) ///
 ((_b[5._at] - _b[24._at])) ///
 ((_b[6._at] - _b[24._at])) ///
 ((_b[7._at] - _b[24._at])) ///
 ((_b[8._at] - _b[24._at])) ///
 ((_b[9._at] - _b[24._at])) ///
 ((_b[10._at] - _b[24._at])) ///
 ((_b[11._at] - _b[24._at])) ///
 ((_b[12._at] - _b[24._at])) ///
 ((_b[13._at] - _b[24._at])) ///
 ((_b[14._at] - _b[24._at])) ///
 ((_b[15._at] - _b[24._at])) ///
 ((_b[16._at] - _b[24._at])) ///
 ((_b[17._at] - _b[24._at])) ///
 ((_b[18._at] - _b[24._at])) ///
 ((_b[19._at] - _b[24._at])) ///
 ((_b[20._at] - _b[24._at])) ///
 ((_b[21._at] - _b[24._at])) ///
 ((_b[22._at] - _b[24._at])) ///
 ((_b[23._at] - _b[24._at])) ///
 ((_b[24._at] - _b[24._at])) ///
 ((_b[25._at] - _b[24._at])) ///
 ((_b[26._at] - _b[24._at])) ///
 ((_b[27._at] - _b[24._at])) ///
 ((_b[28._at] - _b[24._at])) ///
 ((_b[29._at] - _b[24._at])) ///
 ((_b[30._at] - _b[24._at])) ///
 ((_b[31._at] - _b[24._at])) ///
 ((_b[32._at] - _b[24._at])) ///
 ((_b[33._at] - _b[24._at])) ///
 ((_b[34._at] - _b[24._at])) ///
 ((_b[35._at] - _b[24._at])) ///
 ((_b[36._at] - _b[24._at])) ///
 ((_b[37._at] - _b[24._at])) ///
 ((_b[38._at] - _b[24._at])) ///
 ((_b[39._at] - _b[24._at])) ///
 ((_b[40._at] - _b[24._at])) ///
 ((_b[41._at] - _b[24._at])) ///
 ((_b[42._at] - _b[24._at])) ///
 ((_b[43._at] - _b[24._at])) ///
 ((_b[44._at] - _b[24._at])) ///
 ((_b[45._at] - _b[24._at])) ///
 ((_b[46._at] - _b[24._at])) ///
 ((_b[47._at] - _b[24._at])) ///
 ((_b[48._at] - _b[24._at])) ///
 ((_b[49._at] - _b[24._at])) ///
 ((_b[50._at] - _b[24._at])) ///
 ((_b[51._at] - _b[24._at])) ///
 ((_b[52._at] - _b[24._at])) ///
 ((_b[53._at] - _b[24._at])) ///
 ((_b[54._at] - _b[24._at])) ///
 ((_b[55._at] - _b[24._at])) ///
 ((_b[56._at] - _b[24._at])) ///
 ((_b[57._at] - _b[24._at])) ///
 ((_b[58._at] - _b[24._at])) ///
 ((_b[59._at] - _b[24._at])) ///
 ((_b[60._at] - _b[24._at])) ///
 ((_b[61._at] - _b[24._at])),post
 
*THE RESULTS OF THE nlcom COMMAND (I.E., THE PREVALENCE DIFFERENCES) ARE SAVED INTO A .dta FILE FOR LATER ANALYSIS OR GRAPHING.
parmest, saving("Post Estimation Predictions\GLU_RD_HEI_total.dta",replace) 




**************************************************************************************************************************************************
*PART #3 - SENSITIVITY ANALYSIS: REPLACING THE ILIAC WAIST CIRCUMFERENCE MEASURE WITH MIDPOINT WASIT CIRCUMFERENCE MEASURE FOR METABOLIC SYNDROME
**************************************************************************************************************************************************


clear all
use "`new'"

*RUNNING A LOGISTIC REGRESSION MODEL ON ONE IMPUTATION SET
mi extract 1	
		logit  metsx_v5_midwc _S_hei_total1 _S_hei_total2 momage povperc_v1 i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*USING THE MCP COMMAND TO GENERATE THE PREDICTION COMMAND WE WILL USE ON THE MI DATASET BELOW
mcp heiy1_total_v1 (_S_hei_total1- _S_hei_total2), var1(hei_total_range (_S_hei_total1_r-_S_hei_total2_r)) show ci


*NOTE: COPY AND PASTE THE LONG MARGINS COMMAND THAT DISPLAYS FROM RUNNING THIS INTO THE CODE BELOW
*THIS WILL ENSURE THAT YOU'VE GOT A SMOOTH LINE OF PREDICTIONS BELOW
*AFTER PASTING, BE SURE TO ADD "POST" TO THE MIMRGNS COMMAND IN ORDER TO BE ABLE TO SAVE THE ESTIMATES

*RUNNING A LOGISTIC REGRESSION MODEL ACROSS ALL IMPUTATIONS
use "`new'", clear

mi estimate, or saving("Model Estimates\miestfile_total_waistmid", replace) esample(esample) dots:logit  metsx_v5_midwc _S_hei_total1 _S_hei_total2  momage povperc_v1  i.momeduc4 i.momrace4 i.married i.smokerpre i.insurpub i.gravcat  i.puqe_dich i.prediab i.prehtn i.pa_totmetwk3_v1 i.sleepsat i.epds_dich_v1 i.bmiprepregcat6,  vce(robust)

*CALCULATING PREDICTIONS AT SPECIFIC VALUES OF HEI-2020 TOTAL SCORE
*PASTE VALUES FROM THE MCP COMMAND ABOVE INTO THE CODE BELOW, ADD "POST"
mimrgns, predict( pr) post  at( _S_hei_total1=38 _S_hei_total2=0 ) at( _S_hei_total1=38.96666666666667 _S_hei_total2=0 ) at( _S_hei_total1=39.93333333333333 _S_hei_total2=0 ) at( _S_hei_total1=40.9 _S_hei_total2=0 ) at( _S_hei_total1=41.86666666666667 _S_hei_total2=0 ) at( _S_hei_total1=42.83333333333334 _S_hei_total2=0 ) at( _S_hei_total1=43.8 _S_hei_total2=0 ) at( _S_hei_total1=44.76666666666667 _S_hei_total2=0 ) at( _S_hei_total1=45.73333333333333 _S_hei_total2=0 ) at( _S_hei_total1=46.7 _S_hei_total2=0 ) at( _S_hei_total1=47.66666666666666 _S_hei_total2=0 ) at( _S_hei_total1=48.63333333333333 _S_hei_total2=.0000257423252395 ) at( _S_hei_total1=49.6 _S_hei_total2=.0021086127088838 ) at( _S_hei_total1=50.56666666666666 _S_hei_total2=.0116876127210559 ) at( _S_hei_total1=51.53333333333333 _S_hei_total2=.0345328192067307 ) at( _S_hei_total1=52.5 _S_hei_total2=.0764143090108828 ) at( _S_hei_total1=53.46666666666667 _S_hei_total2=.1431021589784871 ) at( _S_hei_total1=54.43333333333334 _S_hei_total2=.2403664459545183 ) at( _S_hei_total1=55.4 _S_hei_total2=.3739772467839501 ) at( _S_hei_total1=56.36666666666667 _S_hei_total2=.5497046383117591 ) at( _S_hei_total1=57.33333333333333 _S_hei_total2=.7733186973829176 ) at( _S_hei_total1=58.3 _S_hei_total2=1.050589500842403 ) at( _S_hei_total1=59.26666666666667 _S_hei_total2=1.38728712553519 ) at( _S_hei_total1=60.23333333333333 _S_hei_total2=1.789181648306252 ) at( _S_hei_total1=61.2 _S_hei_total2=2.262043146000564 ) at( _S_hei_total1=62.16666666666667 _S_hei_total2=2.811641695463102 ) at( _S_hei_total1=63.13333333333333 _S_hei_total2=3.443747373538834 ) at( _S_hei_total1=64.09999999999999 _S_hei_total2=4.164130257072739 ) at( _S_hei_total1=65.06666666666666 _S_hei_total2=4.978305212961531 ) at( _S_hei_total1=66.03333333333333 _S_hei_total2=5.88590851635488 ) at( _S_hei_total1=67 _S_hei_total2=6.880692978232025 ) at( _S_hei_total1=67.96666666666667 _S_hei_total2=7.956155225139035 ) at( _S_hei_total1=68.93333333333334 _S_hei_total2=9.10579188362199 ) at( _S_hei_total1=69.90000000000001 _S_hei_total2=10.32309958022696 ) at( _S_hei_total1=70.86666666666667 _S_hei_total2=11.60157494150001 ) at( _S_hei_total1=71.83333333333334 _S_hei_total2=12.93471459398722 ) at( _S_hei_total1=72.8 _S_hei_total2=14.31601516423465 ) at( _S_hei_total1=73.76666666666667 _S_hei_total2=15.73897327878841 ) at( _S_hei_total1=74.73333333333333 _S_hei_total2=17.19708556419455 ) at( _S_hei_total1=75.7 _S_hei_total2=18.68384864699915 ) at( _S_hei_total1=76.66666666666666 _S_hei_total2=20.19275915374825 ) at( _S_hei_total1=77.63333333333333 _S_hei_total2=21.71731371098797 ) at( _S_hei_total1=78.59999999999999 _S_hei_total2=23.25100894526438 ) at( _S_hei_total1=79.56666666666666 _S_hei_total2=24.78756919226084 ) at( _S_hei_total1=80.53333333333333 _S_hei_total2=26.32420172484829 ) at( _S_hei_total1=81.5 _S_hei_total2=27.86083425743574 ) at( _S_hei_total1=82.46666666666667 _S_hei_total2=29.39746679002319 ) at( _S_hei_total1=83.43333333333334 _S_hei_total2=30.93409932261064 ) at( _S_hei_total1=84.40000000000001 _S_hei_total2=32.47073185519809 ) at( _S_hei_total1=85.36666666666667 _S_hei_total2=34.00736438778555 ) at( _S_hei_total1=86.33333333333334 _S_hei_total2=35.54399692037299 ) at( _S_hei_total1=87.3 _S_hei_total2=37.08062945296042 ) at( _S_hei_total1=88.26666666666667 _S_hei_total2=38.61726198554787 ) at( _S_hei_total1=89.23333333333333 _S_hei_total2=40.15389451813533 ) at( _S_hei_total1=90.2 _S_hei_total2=41.69052705072277 ) at( _S_hei_total1=91.16666666666666 _S_hei_total2=43.22715958331019 ) at( _S_hei_total1=92.13333333333333 _S_hei_total2=44.76379211589766 ) at( _S_hei_total1=93.09999999999999 _S_hei_total2=46.3004246484851 ) at( _S_hei_total1=94.06666666666666 _S_hei_total2=47.83705718107257 ) at( _S_hei_total1=95.03333333333333 _S_hei_total2=49.37368971366001 ) at( _S_hei_total1=96 _S_hei_total2=50.91032224624745 ) 


*SAVING THE PREDICTIONS FOR LATER ANALYSIS OR GRAPHING
parmest, saving("Post Estimation Predictions\RD_HEI_total_waistmid.dta",replace) 
 
 

