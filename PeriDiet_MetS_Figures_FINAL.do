*************************************************************
*		METABOLIC SYNDROME AND PERINATAL DIET PAPER			*
*		WRITTEN BY SARA PARISI AND QIANHUI JIN				*
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


*SET WORKING DIRECTORY
clear all
cd "I:\Bodnar Qianhui\HEI MetSx\Datasets"


**********************************************************************************
*FIGURE 1 RADAR PLOT															 *
*Periconceptional Healthy Eating Index-2020 component scores as a percent of the *
*maximum according to the presence of postpartum metabolic syndrome				 *
**********************************************************************************

ssc install radar
*help radar

*CALL IN THE IMPUTED DATASET*
use "HEI and MetS Imputed_20qj_FINAL_mid_wc.dta",clear

*CREATE "PERCENTAGE OF TOTAL POSSIBLE POINTS" VARIABLES*
gen totalfruitpct = (heiy3_totalfruit_v1/5)*100
gen wholefruitpct = (heiy4_wholefruit_v1/5)*100
gen totalvegpct = (heiy1_totalveg_v1/5)*100
gen greenbeanpct = (heiy2_green_and_bean_v1/5)*100
gen wholegrainpct = (heiy5_wholegrain_v1/10)*100
gen dairypct = (heiy6_totaldairy_v1/10)*100
gen totprotpct = (heiy7_totprot_v1/5)*100
gen seaplantpct = (heiy8_seaplant_prot_v1/5)*100
gen fattyacidpct = (heiy9_fattyacid_v1/10)*100
gen sodiumpct = (heiy10_sodium_v1/10)*100
gen refinedgrainpct = (heiy11_refinedgrain_v1/10)*100
gen addsugarpct = (heiy12_addsug_v1/10)*100
gen satfapct = ( heiy13_sfa_v1/10)*100

rename metsx_v5 metsyn

***SAVE THIS BIG DATASET AS A TEMP DATASET SO YOU CAN GO BACK TO IT AFTER YOU'VE COLLAPSED IT***
tempfile big
save "`big'"

*REMOVE THE MI SETTING SO YOU CAN COLLAPSE IT*
mi unset

*GET THE MEANS OF EACH COMPONENT - BY MET SX - JUST TAKING THE AVERAGE ACROSS ALL 20 IMPUTED DATASETS*
collapse (mean) totalfruitpct-satfapct if mi_m!=0, by(metsyn)


*RENAME THE VARIABLES SO IT CAN BE RESHAPED*
rename  totalfruitpct heiy1
rename  wholefruitpct heiy2
rename  totalvegpct heiy3
rename  greenbeanpct heiy4
rename  wholegrainpct heiy5
rename  dairypct heiy6
rename  totprotpct heiy7
rename  seaplantpct heiy8
rename  fattyacidpct heiy9
rename  refinedgrainpct heiy10
rename  sodiumpct heiy11
rename  addsugarpct heiy12
rename  satfapct heiy13

*RESHAPE TO LONG*
reshape long heiy, i(metsyn) j(comp)

*NOW RESHAPE TO WIDE SO THAT THERE ARE 13 LINES - 1 FOR EACH COMPONENT*
reshape wide heiy, i(comp) j(metsyn)

label define comp 1 "Total Fruits" 2 "Whole Fruits" 3 "Total Vegetables" 4 "Greens & Beans" 5 "Whole Grains" 6 "Dairy" 7 "Total Protein Foods" 8 "Seafood & Plant Proteins" 9 "Fatty Acids" 10 "Refined Grains" 11 "Sodium" 12 "Added Sugars" 13 "Saturated Fats"
label values comp comp

*CONVERT COMP INTO A STRING VARIABLE AS REQUIRED BY THE RADAR COMMAND*
decode comp, gen(heicomp)

*RENAME THE VARIABLES FOR DISPLAY IN THE LEGEND*
rename heiy0 MetSx0
rename heiy1 MetSx1
gen Perfect=100


*MAKE THE RADER PLOT IN STATA*
radar heicomp MetSx0 MetSx1 Perfect, aspect(1) labsize(small) rlabel(0 10 20 30 40 50 60 70 80 90 100) legend(off) graphregion(fcolor(white) lcolor(none) ifcolor(white) ilcolor(white)) note(..., span) lc("17 119 51*1.2" "136 204 238*1.2" pink) lp(solid solid dash)
*must go into the graph editor and remove the axis gridlines**




*********************************************************************************
**FIGURE 2 																		*
*Association between Healthy Eating Index-2020 total score and the adjusted 	*
*absolute prevalence of postpartum metabolic syndrome		 					*
*********************************************************************************

*READ IN THE SAVED RUBIN'S RULES COMBINED MARGINS ESTIMATES FROM MODELING DO FILE*
clear all
use "I:\Bodnar Qianhui\HEI MetSx\Files\Risk datasets\adjusted_mimrgns_HEI_total_bp80.dta"
gen heiy1_total_v1 = 38 + (96 - 38) * (_n - 1) / (61 - 1)
keep heiy1_total_v1 estimate min95 max95
rename estimate ir
rename min95 lb_ir
rename max95 ub_ir
sort  heiy1_total_v1

*GRAPHING THE ESTIMATES*		  
		 twoway  (rcap lb_ir ub_ir heiy1_total_v1 if inrange(heiy1_total_v1, 38,96), lcolor(navy) lwidth(thin)) ///
        (scatter ir heiy1_total_v1 if inrange(heiy1_total_v1, 38,96), mcolor(navy) msize(tiny)) ///
        , ytitle("Adjusted prevalence (95% CI)", size(large) ) ///
          xtitle("HEI-2020 total score", size(large) height(6)) ///
          graphregion(fcolor(white) lcolor(white)) ///
		  xlabel(35 (10) 100, nogrid labsize(4)) ///
		  name(adjust_metsx, replace) ///
		  ylabel(0.00 (0.04) 0.28 , nogrid format(%9.2f) labsize(4)) ///  
	      graphregion(fcolor(white) lcolor(white) margin(medium)) ///
	      aspect(1) ///
	      legend(off) ///
          plotregion(style(none) margin(zero))          
		  



**************************************************************************************
*FIGURE S2 																		     *
*Distribution of Healthy Eating Index–2020 total score according to postpartum MetSx * 
**************************************************************************************

*READ IN THE IMPUTED DATASET* 
clear all
use  "HEI and MetS Imputed_20qj_FINAL_mid_wc.dta"

twoway (histogram heiy1_total_v1 if metsx_v5 == 0, width(2) fcolor("136 204 238*1.2") lcolor("136 204 238*1.2")) ///
       (histogram heiy1_total_v1 if metsx_v5 == 1, width(2) fcolor(none) lcolor(pink*1.6)), ///
	   ytitle("Density", size(large) ) ///
	    xtitle("HEI-2020 total score", size(large) height(6)) ///
		xlabel(20 (10) 100, nogrid labsize(4)) ///
		ylabel(0 (0.01) 0.04, nogrid labsize(4) format(%9.2f)) ///
	    legend(off)                 
   


****************************************************************************************
*FIGURE S3																			   *
*Association between Healthy Eating Index-2020 total score and the adjusted prevalence *
*difference for postpartum metabolic syndrome compared with an HEI score of 60	       *
****************************************************************************************

*READ IN THE SAVED RUBIN'S RULES COMBINED PREVALENCE DIFFERENCE ESTIMATES FROM MODELING DO FILE*
clear all		  
use "I:\Bodnar Qianhui\HEI MetSx\Files\Risk Difference datasets\RD_HEI_total_bp80.dta",clear

gen heiy1_total_v1 = 38 + (96 - 38) * (_n - 1) / (61 - 1)
keep heiy1_total_v1 estimate min95 max95
rename estimate rd
rename min95 lb_rd
rename max95 ub_rd
sort  heiy1_total_v1


twoway  (rcap lb_rd ub_rd heiy1_total_v1 if inrange(heiy1_total_v1, 38,96), lcolor(navy) msize(tiny)) ///
        (scatter rd heiy1_total_v1 if inrange(heiy1_total_v1, 38,96), mcolor(navy) msize(tiny)) ///
        , ytitle("Adjusted prevalence difference (95% CI)", size( large)) ///
          xtitle("HEI-2020 total score", size(large) ) ///
          ylabel(-0.2 (0.05) 0.1, nogrid format(%9.2f) labsize(4) ) ///  
          graphregion(fcolor(white) lcolor(white)) ///
          xlabel(35 (10) 100, nogrid labsize(4) ) ///
          name(metsx_rd, replace) ///
          aspect(1) ///
          graphregion(fcolor(white) lcolor(white) margin(medium)) ///
          plotregion(style(none) margin(zero)) ///
          legend(label(1 "95% CI") label(2 "Estimate ")) ///
		  legend(off) ///
		  yline(0, lcolor(gray))

		  

*****************************************************************************************
*FIGURE S4 Sensitivity analysis 														*
*Comparing the association between Healthy Eating Index-2020 							*
*total score and the adjusted absolute prevalence of postpartum metabolic syndrome, 	*
*when metabolic syndrome is defined using the iliac crest waist circumference measure 	*
*(Panel A) versus the midpoint waist circumference measure (Panel B)					*
*****************************************************************************************

*READ IN THE SAVED RUBIN'S RULES COMBINED MARGINS ESTIMATES FROM MODELING DO FILE*

*PANEL A*
clear all
use "I:\Bodnar Qianhui\HEI MetSx\Files\Risk datasets\adjusted_mimrgns_HEI_total_bp80.dta"

gen heiy1_total_v1 = 38 + (96 - 38) * (_n - 1) / (61 - 1)
keep heiy1_total_v1 estimate min95 max95
rename estimate ir
rename min95 lb_ir
rename max95 ub_ir
sort  heiy1_total_v1

*GRAPHING THE ESTIMATES*		  
		 twoway  (rcap lb_ir ub_ir heiy1_total_v1 if inrange(heiy1_total_v1, 38,96), lcolor(navy) lwidth(thin)) ///
        (scatter ir heiy1_total_v1 if inrange(heiy1_total_v1, 38,96), mcolor(navy) msize(tiny)) ///
        , ytitle("Adjusted prevalence (95% CI)", size(large) ) ///
          xtitle("HEI-2020 total score", size(large) height(6)) ///
          graphregion(fcolor(white) lcolor(white)) ///
		  xlabel(35 (10) 100, nogrid labsize(4)) ///
		  name(adjust_metsx, replace) ///
		  ylabel(0.00 (0.04) 0.28 , nogrid format(%9.2f) labsize(4)) ///  
	      graphregion(fcolor(white) lcolor(white) margin(medium)) ///
	      aspect(1) ///
	      legend(off) ///
          plotregion(style(none) margin(zero)) 
		  
*PANEL B
clear all
use "I:\Bodnar Qianhui\HEI MetSx\Files\Risk datasets\adjusted_mimrgns_HEI_total_bp80_midwc.dta"
gen heiy1_total_v1 = 38 + (96 - 38) * (_n - 1) / (61 - 1)
keep heiy1_total_v1 estimate min95 max95
rename estimate ir
rename min95 lb_ir
rename max95 ub_ir
sort  heiy1_total_v1

*GRAPHING THE ESTIMATES*	  
		 twoway  (rcap lb_ir ub_ir heiy1_total_v1 if inrange(heiy1_total_v1, 38,96), lcolor(navy) lwidth(thin)) ///
        (scatter ir heiy1_total_v1 if inrange(heiy1_total_v1, 38,96), mcolor(navy) msize(tiny)) ///
        , ytitle("Adjusted prevalence (95% CI)", size(large) ) ///
          xtitle("HEI-2020 total score", size(large) height(6)) ///
          graphregion(fcolor(white) lcolor(white)) ///
		  xlabel(35 (10) 100, nogrid labsize(4)) ///
		  name(adjust_metsx, replace) ///
		  ylabel(0.00 (0.04) 0.28 , nogrid format(%9.2f) labsize(4)) ///  
	      graphregion(fcolor(white) lcolor(white) margin(medium)) ///
	      aspect(1) ///
	      legend(off) ///
          plotregion(style(none) margin(zero))          
		 		  
		  
		  
		  
