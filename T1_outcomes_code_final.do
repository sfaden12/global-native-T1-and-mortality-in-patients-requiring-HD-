clear all

********************************************************************************
*** SET DIRECTORIES ***

global fig_filepath "Z:/My Documents/CYCLE/Figures"

global data_filepath "Z:/My Documents/CYCLE/Data"

********************************************************************************
*** SAVE NEW CYCLE DATA (WITH DIABETES CODING) AS .DTA FILE ***

import excel using "${data_filepath}/Copy of CYCLE and glasgow patients July 2024" ///
	, clear firstrow
	
keep in 1/156	
	
save "${data_filepath}/data_copy.dta", replace

********************************************************************************
*** LOAD GLASGOW DATA ***

import excel using "${data_filepath}/glasgow", clear firstrow
gen site=2

********************************************************************************
*** TIDY DATA ***

destring new_age, gen(age_new) // string to number

destring StudyID, gen(Study_ID)

rename Dead Dead1
destring Dead1, gen(Dead)

gen DateDeath=date(Datedeath, "MDY") if Study_ID!=214  
replace DateDeath=date(Datedeath, "DMY") if Study_ID==214 // recorded wrong

format DateDeath %td

gen enroldate1=date(enroldate, "DMY") 
format enroldate1 %td

// need date of study exit for those who are censored
gen exit = DateDeath
replace exit = mdy(04,06,2023) if exit==.
format exit %td

gen TxDate1=date(TxDate, "MDY") 
format TxDate1 %td

keep Study_ID Sex Dialysis_Vintage Diabetes HxTransplant ReceivedTx ///
	Dead DateDeath MRI_LVMi_Dias_Baseline MRI_gl_Native_T1_Baseline ///
	MRI_PWV_Baseline MRI_LVEF_Baseline MRI_gl_pss_Long_Baseline TxDate ///
	MRI_gl_pss_Circ_Baseline enroldate1 site age_new SmokingStatus Ethnicity ///
	exit TxDate1

// alternative functional forms for sensitivity analysis
gen age_sq = age_new^2
gen ln_age = ln(age_new)
rcsgen age_new, gen(age_rcs) df(3)

gen Dialysis_Vintage_sq = Dialysis_Vintage^2
gen ln_Dialysis_Vintage = ln(Dialysis_Vintage)
rcsgen Dialysis_Vintage, gen(Dialysis_Vintage_rcs) df(3)

save "${data_filepath}/glasgow.dta", replace

********************************************************************************
*** MERGE CYCLE AND GLASGOW DATASETS ***

use "${data_filepath}/cycle.dta", clear

// append Glasgow data to CYCLE data
append using "${data_filepath}/glasgow.dta"

keep Study_ID Sex Dialysis_Vintage Diabetes HxTransplant ReceivedTx ///
	Dead MRI_LVMi_Dias_Baseline MRI_gl_Native_T1_Baseline exit TxDate1 ///
	MRI_PWV_Baseline MRI_LVEF_Baseline MRI_gl_pss_Long_Baseline ///
	MRI_gl_pss_Circ_Baseline enroldate1 site age_new SmokingStatus Ethnicity ///
	age_sq ln_age age_rcs* Dialysis_Vintage_sq ln_Dialysis_Vintage ///
	Dialysis_Vintage_rcs*
	
format exit %td	

label define sitelabel 1 "Original" 2 "Glasgow"
label values site sitelabel

save "${data_filepath}/cycle_comb.dta", replace

// remerge new correct diabetes values
drop Diabetes
merge 1:1 Study_ID using "${data_filepath}/data_copy", keepusing(AD)

rename AD Diabetes
label define Diabetes_label 1 "Yes" 2 "No"
label values Diabetes Diabetes_label

********************************************************************************
*** COMPARE BASELINE ***

tab site Sex, missing row

tab site Ethnicity, missing row

tab site SmokingStatus, missing row

tab site Diabetes, missing row

bysort site: su Dialysis_Vintage, detail // Glasgow longer on dialysis at entry
hist Dialysis_Vintage, by(site)

bysort site: su MRI_LVEF_Baseline, detail
hist MRI_LVEF_Baseline, by(site)

bysort site: su MRI_LVMi_Dias_Baseline, detail
hist MRI_LVMi_Dias_Baseline, by(site)

bysort site: su MRI_gl_Native_T1_Baseline, detail
hist MRI_gl_Native_T1_Baseline, by(site)

bysort site: su MRI_gl_pss_Circ_Baseline, detail
hist MRI_gl_pss_Circ_Baseline, by(site)

bysort site: su MRI_gl_pss_Long_Baseline, detail
hist MRI_gl_pss_Long_Baseline, by(site)

tab site ReceivedTx, missing row

tab site Dead, missing row

gen tempyear=year(enroldate1)
tab site tempyear, missing

bysort site: su age_new, detail 
hist age_new, by(site)

gen tempyear2=year(exit)
tab site tempyear2, missing

********************************************************************************
*** TIME-TO-EVENT ANALYSIS ***

// time-varying covariate analysis set up
// generate two obs if transplant received
expand 2 if ReceivedTx==1
bysort Study_ID: gen flag=_n

// change origin to transplant date for second obs
gen origin=enroldate1
replace origin=TxDate1 if ReceivedTx==1 & flag==2

// change exit to transplant date for first obs
gen exit2=exit
replace exit2=TxDate1 if ReceivedTx==1 & flag==1

format origin exit2 %td

// change vital status for first obs
gen Dead1=Dead
replace Dead1=0 if ReceivedTx==1 & flag==1

// change transplant status for first obs
gen ReceivedTx1=ReceivedTx
replace ReceivedTx1=0 if ReceivedTx==1 & flag==1

list Study_ID site flag enroldate1 origin TxDate1 ReceivedTx ReceivedTx1 Dead Dead1 exit exit2

list Study_ID site enroldate1 TxDate1 ReceivedTx ReceivedTx1 Dead exit if TxDate1<enroldate1

// set data for survival analysis
// restrict follow-up to five years
stset exit2, fail(Dead1) id(Study_ID) origin(origin) scale(365.25) ///
	exit(time enroldate1+365.25*5) 
	
su _t if site==1
su _t if site==2

tab site flag if _t<(365.25*5) & _st==1 // transplant counts

tab site _d if _t<(365.25*5) & _st==1 // death counts

// univariate survival models	
stcox MRI_gl_Native_T1_Baseline
stcox MRI_LVEF_Baseline
stcox MRI_LVMi_Dias_Baseline

// analysis model
stcox c.MRI_gl_Native_T1_Baseline ReceivedTx1##c.age_new Diabetes c.Dialysis_Vintage Sex

********************************************************************************
*** SENSITIVITY ANALYSIS ***

local MRI_list "MRI_gl_Native_T1_Baseline MRI_LVMi_Dias_Baseline MRI_LVEF_Baseline"

frame create sens_a_frame_g5
frame change sens_a_frame_g5

set obs 13

gen MRI_gl_Native_T1_Baseline=.
gen MRI_gl_Native_T1_Baseline_lci=.
gen MRI_gl_Native_T1_Baseline_uci=.

gen MRI_LVMi_Dias_Baseline=.
gen MRI_LVMi_Dias_Baseline_lci=.
gen MRI_LVMi_Dias_Baseline_uci=.

gen MRI_LVEF_Baseline=.
gen MRI_LVEF_Baseline_lci=.
gen MRI_LVEF_Baseline_uci=.

frame change default

foreach MRI in `MRI_list' {
	
	// 1
	stcox c.`MRI' ReceivedTx1 c.age_new
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 1
		replace `MRI'_lci = r(table)[5,1] in 1
		replace `MRI'_uci = r(table)[6,1] in 1
	}

	// 2
	stcox c.`MRI' ReceivedTx1##c.age_new 
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 2
		replace `MRI'_lci = r(table)[5,1] in 2
		replace `MRI'_uci = r(table)[6,1] in 2
	}
	
	// 3
	stcox c.`MRI' ReceivedTx1 c.age_new c.age_sq
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 3
		replace `MRI'_lci = r(table)[5,1] in 3
		replace `MRI'_uci = r(table)[6,1] in 3
	}
	
	// 4
	stcox c.`MRI' ReceivedTx1 c.ln_age
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 4
		replace `MRI'_lci = r(table)[5,1] in 4
		replace `MRI'_uci = r(table)[6,1] in 4
	}
	
	// 5
	stcox c.`MRI' ReceivedTx1 age_rcs*
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 5
		replace `MRI'_lci = r(table)[5,1] in 5
		replace `MRI'_uci = r(table)[6,1] in 5
	}
	
	// 6
	stcox c.`MRI' ReceivedTx1 c.age_new Sex
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 6
		replace `MRI'_lci = r(table)[5,1] in 6
		replace `MRI'_uci = r(table)[6,1] in 6
	}
	
	// 7
	stcox c.`MRI' ReceivedTx1 c.age_new Sex Diabetes 
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 7
		replace `MRI'_lci = r(table)[5,1] in 7
		replace `MRI'_uci = r(table)[6,1] in 7
	}
	
	// 8
	stcox c.`MRI' ReceivedTx1 c.age_new Sex Diabetes c.Dialysis_Vintage
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 8
		replace `MRI'_lci = r(table)[5,1] in 8
		replace `MRI'_uci = r(table)[6,1] in 8
	}

	// 9
	stcox c.`MRI' ReceivedTx1##c.age_new Diabetes c.Dialysis_Vintage Sex
	
	estimates store `MRI'_6
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 9
		replace `MRI'_lci = r(table)[5,1] in 9
		replace `MRI'_uci = r(table)[6,1] in 9
	}
	
	// 10
	stcox c.`MRI' ReceivedTx1 c.age_new##Diabetes c.Dialysis_Vintage Sex
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 10
		replace `MRI'_lci = r(table)[5,1] in 10
		replace `MRI'_uci = r(table)[6,1] in 10
	}

	// 11
	stcox c.`MRI' ReceivedTx1 c.age_new Sex Diabetes c.Dialysis_Vintage ///
		c.Dialysis_Vintage_sq
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 11
		replace `MRI'_lci = r(table)[5,1] in 11
		replace `MRI'_uci = r(table)[6,1] in 11
	}
	
	// 12
	stcox c.`MRI' ReceivedTx1 c.age_new Sex Diabetes c.ln_Dialysis_Vintage
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 12
		replace `MRI'_lci = r(table)[5,1] in 12
		replace `MRI'_uci = r(table)[6,1] in 12
	}
	
	// 13
	stcox c.`MRI' ReceivedTx1 c.age_new Sex Diabetes c.Dialysis_Vintage_rcs*
	
	frame sens_a_frame_g5 {
		replace `MRI' = r(table)[1,1] in 13
		replace `MRI'_lci = r(table)[5,1] in 13
		replace `MRI'_uci = r(table)[6,1] in 13
	}
	
}

frame sens_a_frame_g5: format MRI* %9.3f

frame sens_a_frame_g5: save "${data_filepath}/sens_analysis_g5.dta", replace

foreach MRI in `MRI_list' {
	estimates restore `MRI'_6
}

********************************************************************************
*** SENSITIVITY ANALYSIS FOREST PLOTS ***

frame sens_a_frame_g5: gen m=_n

frame sens_a_frame_g5 {
	
	twoway (rspike MRI_gl_Native_T1_Baseline_uci MRI_gl_Native_T1_Baseline_lci m ///
				,                                                                ///
				horizontal)                                                      ///
		   (scatter m MRI_gl_Native_T1_Baseline                                  ///
				,                                                                ///
				mcolor(green))                                                   ///
				,                                                                ///
				name(MRI_gl_Native_T1_Baseline_g5, replace)                      ///
				legend(off)                                                      ///
				title(Global native T1 (ms))                                     ///
				yscale(range(1 13))                                              ///
				ylabel(1(1)13)                                                   ///
				ytitle(Model)                                                    ///
				xtitle(Hazard Ratio)                                             ///
				xlabel(0.940 1.008 1.040)                                        ///
				xline(1, lcolor(red) lpattern(solid))
				
	graph export "${fig_filepath}/MRI_gl_Native_T1_Baseline_g5.pdf", replace	
	graph export "${fig_filepath}/MRI_gl_Native_T1_Baseline_g5.png", replace	

	twoway (rspike MRI_LVMi_Dias_Baseline_uci MRI_LVMi_Dias_Baseline_lci m ///
				,                                                          ///
				horizontal)                                                ///
		   (scatter m MRI_LVMi_Dias_Baseline                               ///
				,                                                          ///
				mcolor(green))                                             ///
				,                                                          ///
				name(MRI_LVMi_Dias_Baseline_g5, replace)                   ///
				legend(off)                                                ///
				title(LV mass index (g/m2))                                ///
				yscale(range(1 13))                                        ///
				ylabel(1(1)13)                                             ///
				ytitle(Model)                                              ///
				xtitle(Hazard Ratio)                                       ///
				xlabel(0.940 1.022 1.040)                                  ///
				xline(1, lcolor(red) lpattern(solid))
	
	graph export "${fig_filepath}/MRI_LVMi_Dias_Baseline_g5.pdf", replace
	graph export "${fig_filepath}/MRI_LVMi_Dias_Baseline_g5.png", replace
		
	twoway (rspike MRI_LVEF_Baseline_uci MRI_LVEF_Baseline_lci m ///
				,                                                ///
				horizontal)                                      ///
		   (scatter m MRI_LVEF_Baseline                          ///
				,                                                ///
				mcolor(green))                                   ///
				,                                                ///
				name(MRI_LVEF_Baseline_g5, replace)              ///
				legend(off)                                      ///
				title(LV ejection fraction (%))                  ///
				yscale(range(1 13))                              ///
				ylabel(1(1)13)                                   ///
				ytitle(Model)                                    ///
				xtitle(Hazard Ratio)                             ///
				xlabel(0.940 0.971 1.040)                        ///
				xline(1, lcolor(red) lpattern(solid))
				
	graph export "${fig_filepath}/MRI_LVEF_Baseline_g5.pdf", replace
	graph export "${fig_filepath}/MRI_LVEF_Baseline_g5.png", replace
			
	graph combine MRI_gl_Native_T1_Baseline_g5 MRI_LVMi_Dias_Baseline_g5 ///
				  MRI_LVEF_Baseline_g5                					 ///
					, 													 ///
					row(1) xcommon
		
	graph export "${fig_filepath}/MRI_all_g5.pdf", replace
	graph export "${fig_filepath}/MRI_all_g5.png", replace
	
}		

********************************************************************************
*** KAPLAN MEIER PLOTS ***

// t1
recode MRI_gl_Native_T1_Baseline (min/1185=1) (1185/1245=2) (1245/max=3) ///
	, gen(t1_group)

tab t1_group, missing

sts graph,                                                      ///
	by(t1_group) name(t1_km_g5, replace)                        ///
	legend(order(1 "Below normal range" 2 "Within normal range" ///
				 3 "Above normal range")                        ///
	ring(0) pos(8))
	
graph export "${fig_filepath}/t1_km_g5.pdf", replace
graph export "${fig_filepath}/t1_km_g5.png", replace

gen t1_group2=1 if t1_group==1 | t1_group==3
replace t1_group2=2 if t1_group==2

sts graph,                                                         ///
	by(t1_group2) name(t1_km2_g5, replace)                         ///
	legend(order(2 "Within normal range" 1 "Outside normal range") ///
	ring(0) pos(8))		
	
graph export "${fig_filepath}/t1_km2_g5.pdf", replace
graph export "${fig_filepath}/t1_km2_g5.png", replace

recode MRI_gl_Native_T1_Baseline ///
	(min/1220=1) (1220/1260=2) (1260/1300=3) (1300/1340=4) (1340/max=5) ///
	, gen(t1_group_small)

sts graph,                                                                      ///
	by(t1_group_small) name(t1_km3_g5, replace)                                 ///
	legend(order(1 "<1220" 2 "1220-1260" 3 "1260-1300" 4 "1300-1340" 5 "1340<") ///
	ring(0) pos(8))		
	
graph export "${fig_filepath}/t1_km3_g5.pdf", replace	
graph export "${fig_filepath}/t1_km3_g5.png", replace

// lvmi
recode MRI_LVMi_Dias_Baseline (min/39=1) (39/85=2) (85/max=3), gen(lvmi_m_group)
recode MRI_LVMi_Dias_Baseline (min/30=1) (30/68=2) (68/max=3), gen(lvmi_f_group)

tab lvmi_m_group if Sex==1, missing
tab lvmi_f_group if Sex==2, missing

sts graph if Sex==1,                                            ///
	by(lvmi_m_group) name(lvmi_m_km_g5, replace)                ///
	legend(order(1 "Below normal range" 2 "Within normal range" ///
				 3 "Above normal range")                        ///
	ring(0) pos(8))
	
graph export "${fig_filepath}/lvmi_m_km_g5.pdf", replace
graph export "${fig_filepath}/lvmi_m_km_g5.png", replace
	
sts graph if Sex==2,                                            ///
	by(lvmi_f_group) name(lvmi_f_km_g5, replace)                ///
	legend(order(1 "Below normal range" 2 "Within normal range" ///
				 3 "Above normal range")                        ///
	ring(0) pos(8))

graph export "${fig_filepath}/lvmi_f_km_g5.pdf", replace	
graph export "${fig_filepath}/lvmi_f_km_g5.png", replace
	
gen lvmi_m_group2=1 if lvmi_m_group==1 | lvmi_m_group==3
replace lvmi_m_group2=2 if lvmi_m_group==2

sts graph if Sex==1, 											   ///
	by(lvmi_m_group2) name(lvmi_m_km2_g5, replace) 				   ///
	legend(order(2 "Within normal range" 1 "Outside normal range") ///
	ring(0) pos(8))
	
graph export "${fig_filepath}/lvmi_m_km2_g5.pdf", replace
graph export "${fig_filepath}/lvmi_m_km2_g5.png", replace

gen lvmi_f_group2=1 if lvmi_f_group==1 | lvmi_f_group==3
replace lvmi_f_group2=2 if lvmi_f_group==2

sts graph if Sex==2,                                               ///
	by(lvmi_f_group2) name(lvmi_f_km2_g5, replace)                 ///
	legend(order(2 "Within normal range" 1 "Outside normal range") ///
	ring(0) pos(8))
	
graph export "${fig_filepath}/lvmi_f_km2_g5.pdf", replace
graph export "${fig_filepath}/lvmi_f_km2_g5.png", replace

// lvef
recode MRI_LVEF_Baseline (min/49=1) (49/79=2) (79/max=3), gen(lvef_m_group)
recode MRI_LVEF_Baseline (min/52=1) (52/79=2) (79/max=3), gen(lvef_f_group)

tab lvef_m_group if Sex==1, missing
tab lvef_f_group if Sex==2, missing

sts graph if Sex==1,                                               ///
	by(lvef_m_group) name(lvef_m_km_g5, replace)                   ///
	legend(order(2 "Within normal range" 1 "Outside normal range") ///
	ring(0) pos(8))	
	
graph export "${fig_filepath}/lvef_m_km_g5.pdf", replace
graph export "${fig_filepath}/lvef_m_km_g5.png", replace
	
sts graph if Sex==2, 											   ///
	by(lvef_f_group) name(lvef_f_km_g5, replace) 	     		   ///
	legend(order(2 "Within normal range" 1 "Outside normal range") ///
	ring(0) pos(8))
	
graph export "${fig_filepath}/lvef_f_km_g5.pdf", replace	
graph export "${fig_filepath}/lvef_f_km_g5.png", replace	

********************************************************************************
