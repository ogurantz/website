

*************************************************
/* Specify Folders */
*************************************************

global results_path = "G:\XXX\results"

global data_path = "G:\XXX\Data\transfer240625"


*************************************************
/* Load Data */
*************************************************

* treatment group
cd "$data_path"
import delimited "TestSamp_24JUN2024.csv", clear
g treat=1
tempfile treat
save `treat'

* control group
forv i=1/4 {
	import delimited "CtrlSamp_24JUN2024_Part`i'.csv", clear
	g treat=0
	tempfile d`i'
	save `d`i''
}

* combine
use `treat', clear
forv i=1/4 {
append using `d`i''
}

*************************************************
/* Functions */
*************************************************


capture program drop export_regression_results
program define export_regression_results
    syntax, filename(string) groupn(string) subgroupn(string) modeln(string) 

	preserve
    * Check if regression results are present
    capture confirm matrix r(table)
    if _rc {
        display "No regression results found in r(). Please run a regression first."
        exit 198
    }

	* control mean
	lincom _cons
    scalar b_cm = r(estimate)
    scalar se_cm = r(se)
    scalar ci95_cm = invttail(e(df_r),0.025) * se_cm
    scalar ci_low_cm = b_cm - ci95_cm
    scalar ci_high_cm = b_cm + ci95_cm
    scalar p_cm = r(p)
	scalar n_cm = e(N)
	
	* treatment mean
	lincom _cons + treat
    scalar b_tm = r(estimate)
    scalar se_tm = r(se)
    scalar ci95_tm = invttail(e(df_r),0.025) * se_tm
    scalar ci_low_tm = b_tm  - ci95_tm
    scalar ci_high_tm = b_tm + ci95_tm
    scalar p_tm = r(p)
    scalar n_tm = e(N)
	
	* effect
	lincom treat
    scalar b_tau = r(estimate)
    scalar se_tau = r(se)
    scalar ci95_tau = invttail(e(df_r),0.025) * se_tau
    scalar ci_low_tau = b_tau  - ci95_tau
    scalar ci_high_tau = b_tau + ci95_tau
    scalar p_tau = r(p)
	scalar n_tau = e(N)
	
    * effect over control baseline 
    nlcom (b_nlcom: (_b[treat] / _b[_cons])*100)

    * Store nlcom results and compute confidence intervals
    scalar b_taup = el(r(b),1,1)
    scalar se_taup = sqrt(el(r(V),1,1))
    scalar ci95_taup = invttail(`e(df_r)',0.025) * se_taup
    scalar ci_low_taup = b_taup - ci95_taup
    scalar ci_high_taup = b_taup + ci95_taup
    scalar p_taup = 2 * normal(-abs(b_taup / se_taup))
	scalar n_taup = e(N)
	
	tempname mypost
    tempfile results
    postfile `mypost' str24 group str24 subgroupn str24 modeln str24 estn  pe se ci_lower ci_upper p_value obs   using `results' 
    post `mypost' ("`groupn'") ("`subgroupn'") ("`modeln'") ("co_mean")        (b_cm) (se_cm) (ci_low_cm) (ci_high_cm) (p_cm) (n_cm)
	post `mypost' ("`groupn'") ("`subgroupn'") ("`modeln'") ("tr_mean")        (b_tm) (se_tm) (ci_low_tm) (ci_high_tm) (p_tm) (n_tm)
    post `mypost' ("`groupn'") ("`subgroupn'") ("`modeln'") ("tau")         (b_tau) (se_tau) (ci_low_tau) (ci_high_tau) (p_tau)  (n_tau)
    post `mypost' ("`groupn'") ("`subgroupn'") ("`modeln'") ("tau_percent") (b_taup) (se_taup) (ci_low_taup) (ci_high_taup) (p_taup)  (n_taup)
    postclose `mypost'

    use `results', clear
    capture noisily confirm file "`filename'.dta"
    if _rc { 
        save "`filename'.dta", replace
    }
    else {
        append using "`filename'"
        save "`filename'.dta", replace
    }
	restore
end



*************************************************
/* Preprocess */
*************************************************

* treatment 
label var treat "Treatment Group"

* letter delivered
g       mdeliv = treat
replace mdeliv = 0 if mailing_creation_dt==""|forward=="NO"|(forward=="YES"&return_date_2!="")
label var mdeliv "Letter Delivered"

** outcomes
* n400 submission date
gen double natz_file_dt_d = date(natz_file_dt,"MDY")
format natz_file_dt_d %td

g       n400=natz_file_dt_d!=.
replace n400 = n400*100

* online / mail
gen double n400_online = (ndtmeans_desc == "EFile")*100
gen double n400_mail   = (ndtmeans_desc == "Lockbox")*100

* fee waiver
gen double n400_feew_full       = (ndtpay == "Waived Fee" )*100
gen double n400_feew_partial    = (ndtpay == "Partial Fee"  )*100
gen double n400_feew_any        = (ndtpay == "Partial Fee" | ndtpay == "Waived Fee" )*100
gen double n400_fullfee         = (ndtpay == "Full Fee")*100


** randomization blocks
egen rblock=group(lpr_group coa_group icob_group)

** covariates
gen     xx_coa_group = coa_group
replace xx_coa_group = "Employment" if xx_coa_group=="EMPLOY"
replace xx_coa_group = "Family"     if xx_coa_group=="FAMILY"
replace xx_coa_group = "Parent"     if xx_coa_group=="PARENT"
replace xx_coa_group = "Spouse"     if xx_coa_group=="SPOUSE"
replace xx_coa_group = "Humanitarian"     if xx_coa_group=="HUMANI"

gen     xx_icob_group = icob_group
replace xx_icob_group="Other" if xx_icob_group=="Africa, Asia, Europe, Oceania, Other"

gen     xx_igender = ""
replace xx_igender="Male"   if igender=="M"
replace xx_igender="Female" if igender=="F"

gen     age = 2023 - dob_year

gen     xx_age_cat = ""
replace xx_age_cat = "19-29 years" if age >=19 & age <=29
replace xx_age_cat = "30-39 years" if age >=30 & age <=39
replace xx_age_cat = "40-49 years" if age >=40 & age <=49
replace xx_age_cat = "50-59 years" if age >=50 & age <=59
replace xx_age_cat = "60-69 years" if age >=60 & age <=69
replace xx_age_cat = "70-79 years" if age >=70 & age <=79
replace xx_age_cat = "80+ years"   if age >=80 & age <=150

gen     xx_lpr_group = lpr_group
replace xx_lpr_group = "9 years" if xx_lpr_group=="LPR for 9 years"
replace xx_lpr_group = "8 years" if xx_lpr_group=="LPR for 8 years"
replace xx_lpr_group = "7 years" if xx_lpr_group=="LPR for 7 years"
replace xx_lpr_group = "6 years" if xx_lpr_group=="LPR for 6 years"

gen     xx_ims = ""
replace xx_ims = "Divorced/Widowed" if ims=="Widowed" | ims=="Divorced"
replace xx_ims = "Single"           if ims=="Single" 
replace xx_ims = "Married"          if ims=="Married"


* turn numeric for interactions 
foreach x of varlist xx_* {
	encode `x', gen(`x'_n)
}


*************************************************
/* merge in zip code population estimates */
preserve
import delimited "G:\XXX\Data\ACSDT5Y2023.B01003-Data.csv", clear

drop if inlist(_n,1,2)
rename (v2 v3) (badzip population)
destring population, replace 
g izip=substr(badzip, 7, .)
destring izip, replace
keep izip population
tempfile pop
save `pop'
restore

merge m:1 izip using `pop', gen(merge_pop) keep(1 3)


*************************************************
/* merge in medican zip code income */
merge m:1 izip using "G:\XXX\Data\medinc.dta"
drop if _merge==2


*************************************************
/* Data for concentration of LPRs per zip */
*************************************************

rename population pop
g pop2=pop^2
g pop3=pop^3

g medinc2=medinc^2
g medinc3=medinc^3

bysort izip: egen lpr_in_zip=count(izip)
g concentration_per_zip=lpr_in_zip/pop
egen tag1=tag(izip)
su concentration_per_zip if tag1==1, d

* Topcode at 99th percentile
replace concentration_per_zip=.0351902 if concentration_per_zip>.0351902&concentration_per_zip!=. 

* Create interaction given median zip is .0018031
replace concentration_per_zip=(concentration_per_zip-.01)*100
g treatc=treat*concentration_per_zip


*************************************************
/* Checks */
*************************************************
isid id

* PS constant across blocks at .161
tabstat treat , by(rblock) stats(mean n)


*************************************************
/* Main Regs Fig 2 and 3  */
*************************************************
cd "$results_path"

* Full Sample
capture rm reg_results.dta
reghdfe n400 treat ,  vce(robust) absorb(rblock)
export_regression_results, filename("reg_results") groupn("Full Sample") modeln("ITT") subgroupn("Full Sample")


* Class Of Admission
levelsof xx_coa_group, local(levels)
 foreach l of local levels {
 reghdfe n400 treat ,  vce(robust) absorb(rblock) , if xx_coa_group == "`l'"
 export_regression_results, filename("reg_results") groupn("Class of Admission") modeln("ITT") subgroupn("`l'")
}

* Region of Birth
levelsof xx_icob_group, local(levels)
 foreach l of local levels {
 reghdfe n400 treat ,  vce(robust) absorb(rblock) , if xx_icob_group == "`l'"
 export_regression_results, filename("reg_results") groupn("Region of Birth") modeln("ITT") subgroupn("`l'")
}

* LPR Years
levelsof xx_lpr_group, local(levels)
 foreach l of local levels {
 reghdfe n400 treat ,  vce(robust) absorb(rblock) , if xx_lpr_group == "`l'"
 export_regression_results, filename("reg_results") groupn("Years as LPR") modeln("ITT") subgroupn("`l'")
}

* Gender
levelsof xx_igender, local(levels)
 foreach l of local levels {
 reghdfe n400 treat ,  vce(robust) absorb(rblock) , if xx_igender == "`l'"
 export_regression_results, filename("reg_results") groupn("Gender") modeln("ITT") subgroupn("`l'")
}

* Age
levelsof  xx_age_cat, local(levels)
 foreach l of local levels {
 reghdfe n400 treat ,  vce(robust) absorb(rblock) , if xx_age_cat == "`l'"
 export_regression_results, filename("reg_results") groupn("Age") modeln("ITT") subgroupn("`l'")
}

* Marital Status
levelsof xx_ims, local(levels)
 foreach l of local levels {
 reghdfe n400 treat ,  vce(robust) absorb(rblock) , if xx_ims == "`l'"
 export_regression_results, filename("reg_results") groupn("Marital Status") modeln("ITT") subgroupn("`l'")
}


*************************************************
/* Fee Waiver and Online Regressions Fig 3 */
*************************************************

capture rm reg_results_online.dta
reghdfe n400_online treat ,  vce(robust) absorb(rblock)
export_regression_results, filename("reg_results_online") groupn("Mode of Filing") modeln("ITT") subgroupn("Online")

reghdfe n400_mail treat ,  vce(robust) absorb(rblock)
export_regression_results, filename("reg_results_online") groupn("Mode of Filing") modeln("ITT") subgroupn("Mail")

reghdfe n400_fullfee treat ,  vce(robust) absorb(rblock)
export_regression_results, filename("reg_results_online") groupn("Application Fee") modeln("ITT") subgroupn("Full Fee")

reghdfe n400_feew_any treat ,  vce(robust) absorb(rblock)
export_regression_results, filename("reg_results_online") groupn("Application Fee") modeln("ITT") subgroupn("Any Fee Waiver")

reghdfe n400_feew_partial treat ,  vce(robust) absorb(rblock)
export_regression_results, filename("reg_results_online") groupn("Application Fee") modeln("ITT") subgroupn("Partial Fee Waiver")

reghdfe n400_feew_full treat ,  vce(robust) absorb(rblock)
export_regression_results, filename("reg_results_online") groupn("Application Fee") modeln("ITT") subgroupn("Full Fee Waiver")




*************************************************
/* SM Materials */
*************************************************


*************************************************
/* Summary Stats */
*************************************************

*** label vars
* outcomes
gen       yy_n400 = n400/100
label var yy_n400 "N-400 Submitted"

gen       yy_n400_online = n400_online/100
label var yy_n400_online "N-400 Submitted: Online"

gen       yy_n400_mail = n400_mail/100
label var yy_n400_mail "N-400 Submitted: Mail"

gen       yy_n400_fullfee = n400_fullfee/100
label var yy_n400_fullfee "N-400 Submitted: Full Fee"

gen       yy_n400_feew_full = n400_feew_full/100
label var yy_n400_feew_full "N-400 Submitted: Full Fee Waiver"

gen       yy_n400_feew_partial = n400_feew_partial/100
label var yy_n400_feew_partial "N-400 Submitted: Partial Fee Waiver"

gen       yy_n400_feew_any = n400_feew_any/100
label var yy_n400_feew_any "N-400 Submitted: Any Fee Waiver"

* covariates

* strata variables
tab   xx_lpr_group, gen(xxxs_lpr_group)
label var xxxs_lpr_group1 "Years as LPR: 6"
label var xxxs_lpr_group2 "Years as LPR: 7"
label var xxxs_lpr_group3 "Years as LPR: 8"
label var xxxs_lpr_group4 "Years as LPR: 9"

tab xx_coa_group_n, gen(xxxs_coa_group)
label var xxxs_coa_group1 "Class of Admission: Employment"
label var xxxs_coa_group2 "Class of Admission: Family"
label var xxxs_coa_group3 "Class of Admission: Humanitarian"
label var xxxs_coa_group4 "Class of Admission: Parent"
label var xxxs_coa_group5 "Class of Admission: Spouse"

tab xx_icob_group, gen(xxxs_icob_group)
label var xxxs_icob_group1 "Region of Birth: Caribbean"
label var xxxs_icob_group2 "Region of Birth: Central/South America"
label var xxxs_icob_group3 "Region of Birth: China/India"
label var xxxs_icob_group4 "Region of Birth: Mexico"
label var xxxs_icob_group5 "Region of Birth: Other"

* non strata variables
gen xxx_male = xx_igender=="Male"
label var xxx_male "Male"

gen xxx_age = age
label var xxx_age "Age"

tab xx_ims, gen(xxx_ims)
replace xxx_ims1 = 0 if xx_ims==""
replace xxx_ims2 = 0 if xx_ims==""
replace xxx_ims3 = 0 if xx_ims==""

label var xxx_ims1 "Marital Status: Divorced/Widowed"
label var xxx_ims2 "Marital Status: Married"
label var xxx_ims3 "Marital Status: Single"

gen       xxx_state_CA = istate=="California"
label var xxx_state_CA  "State: California"

gen       xxx_state_NY = istate=="New York"
label var xxx_state_NY  "State: New York "

gen       xxx_state_TX = istate=="Texas"
label var xxx_state_TX  "State: Texas"

gen       xxx_state_FL = istate=="Florida"
label var xxx_state_FL  "State: Florida"

gen       xxx_state_NJ = istate=="New Jersey"
label var xxx_state_NJ  "State: New Jersey"

gen       xxx_propensity = propensity
label var xxx_propensity "Propensity to Refile"

gen       xxx_prefiled = natz_file_dt_d < td(26oct2023)
label var xxx_prefiled "Filed N-400 Prior to Letter Initiation"

gen xxx_medinc = medinc
label var xxx_medinc "Median Income of Zip Code (\\$)"

* table for summary stats
cd "$results_path"
cap texdoc close
cap erase "sum_stats.tex"
 texdoc init "sum_stats.tex", force

tex \begin{table}
tex \centering
tex \caption{Summary Statistics \label{tab:summarystats}}
tex \begin{tabular}{l*{2}{c}}
tex \hline\hline
tex  & \multicolumn{1}{c}{Mean} &  SD  \\
tex \hline

* outcomes
tex & \multicolumn{2}{c}{Outcomes} \\
tex \cmidrule{2-3} 

foreach x of varlist yy_* {

                local label : variable label `x'
                sum `x'
                local o_mean: di %9.3f r(mean)
                local o_sd: di %9.3f r(sd)
                local o_n: di %9.0fc r(N)

tex  `label' & `o_mean' & `o_sd'  \\

}

tex \hline \\
tex & \multicolumn{2}{c}{Treatments} \\
tex \cmidrule{2-3} 

foreach x of varlist treat mdeliv {

                local label : variable label `x'
                sum `x'
                local o_mean: di %9.3f r(mean)
                local o_sd: di %9.3f r(sd)
                local o_n: di %9.0fc r(N)

tex  `label' & `o_mean' & `o_sd'  \\

}
tex \hline \\
tex & \multicolumn{2}{c}{Covariates (Strata)} \\
tex \cmidrule{2-3} 

foreach x of varlist xxxs_* {

                local label : variable label `x'
                sum `x'
                local o_mean: di %9.3f r(mean)
                local o_sd: di %9.3f r(sd)
                local o_n: di %9.0fc r(N)

tex  `label' & `o_mean' & `o_sd'  \\

}
tex \hline \\
tex & \multicolumn{2}{c}{Covariates (Non-Strata)} \\
tex \cmidrule{2-3} 

foreach x of varlist xxx_* {

                local label : variable label `x'
                sum `x'
                local o_mean: di %9.3f r(mean)
                local o_sd: di %9.3f r(sd)
                local o_n: di %9.0fc r(N)

tex  `label' & `o_mean' & `o_sd'  \\

}

tex \hline\hline
tex \multicolumn{3}{p{.4\textwidth}}{\scriptsize Means and standard deviations shown. Strata covariates are the covariates that form the randomization blocks. The Propensity to Refile is a variable calculated internally by USCIS. Due to delays between when N-400 data is uploaded to USCIS servers and the timing of randomization and letter mailings, approximately 1.3\% of the randomization pool had already submitted the Form N-400 before the first letters were mailed on October 26, 2023. These pre-submissions are captured by the Filed N-400 Prior to Letter Initiation covariate. (N=1,861,178).}
tex \end{tabular}
tex \end{table}

*************************************************
/* Balance Checks */
*************************************************

cd "$results_path"
cap texdoc close
cap erase "balance_stats.tex"
 texdoc init "balance_stats.tex", force

tex \begin{table}
tex \centering
tex \caption{Balance Checks \label{tab:balancestats}}
tex \begin{tabular}{l*{6}{c}}
tex \hline\hline

tex & \multicolumn{2}{c}{Means} & \multicolumn{2}{c}{Difference} & \\
tex  &   Treatment &       Control & Raw & Standardized & p-value  \\
tex   & (1) & (2) & (3) & (4) &(5)  \\
tex \hline
tex \\

tex & \multicolumn{5}{c}{Strata Covariates} \\
tex \cmidrule{2-6} \\

foreach x of varlist xxxs_* {

                local label : variable label `x'
                * pooled SD			
			    sum `x'
                local o_sd: di %9.3f r(sd)
			
                * treated
                sum `x' if treat==1
                local t_mean: di %9.3f r(mean)

                * controls
                sum `x' if treat==0
                local c_mean: di %9.3f r(mean)
			
			    * delta raw
				local delta_raw = `t_mean' - `c_mean'
				local delta_raw: di %9.3f `delta_raw'
				
				* delta std
			    local delta_std = (`t_mean' - `c_mean')/ (`o_sd')
				local delta_std: di  %9.3f `delta_std'
			
tex  `label' & `t_mean' & `c_mean' & `delta_raw' & `delta_std' & - \\

}

tex & \multicolumn{5}{c}{Non-Strata Covariates} \\
tex \cmidrule{2-6} \\

foreach x of varlist xxx_* {

                local label : variable label `x'
			
			* pooled sd
			    sum `x'
                local o_sd: di %9.3f r(sd)
			
                * treated
                sum `x' if treat==1
                local t_mean: di %9.3f r(mean)

                * controls
                sum `x' if treat==0
                local c_mean: di %9.3f r(mean)
			
				* delta raw
				local delta_raw = `t_mean' - `c_mean'
				local delta_raw: di %9.3f `delta_raw'
				
				* delta std
			    local delta_std = (`t_mean' - `c_mean')/ (`o_sd')
				local delta_std: di  %9.3f `delta_std'
			
                * p-value
				reghdfe `x' treat ,  vce(robust) absorb(rblock)
                lincom treat
                local p_mean: di %9.3f r(p)

tex  `label' & `t_mean' & `c_mean' & `delta_raw' & `delta_std' & `p_mean' \\

}



tex \hline\hline
tex \multicolumn{6}{p{.7\textwidth}}{\scriptsize The table provides mean values for both the treatment group (n=300,000) and the control group (n=1,561,178). Standardized differences are calculated as raw differences divided by the pooled standard deviation. p-values are derived from regressions that estimate differences between the treatment and control group, including fixed effects for randomization strata.  Strata covariates are the covariates that form the randomization blocks. No p-values are showsn for strata covariates since they are balanced by design. The Propensity to Refile is an internally calculated value by USCIS. Due to delays between the upload of N-400 data to USCIS servers and the timing of randomization and letter mailings, approximately 1.3\% of the randomization pool had already submitted the Form N-400 by October 26, 2023. These pre-submissions are captured by the Filed N-400 Prior to Letter Initiation covariate.}
tex \end{tabular}
tex \end{table}



*************************************************
/* Test for Significant Differences in ITTs */
*************************************************

* class of admission, origin, and LPR years
cap estimates drop model*
reghdfe n400 treat c.treat#ib3.xx_coa_group_n ,  vce(robust) absorb(rblock)
estimates store model1

lincom 1.xx_coa_group_n#c.treat - 2.xx_coa_group_n#c.treat
lincom 1.xx_coa_group_n#c.treat - 4.xx_coa_group_n#c.treat

lincom 5.xx_coa_group_n#c.treat - 2.xx_coa_group_n#c.treat
lincom 5.xx_coa_group_n#c.treat - 4.xx_coa_group_n#c.treat


reghdfe n400 treat c.treat##ib3.xx_icob_group_n ,  vce(robust) absorb(rblock)
estimates store model2

lincom 2.xx_icob_group_n#c.treat - 1.xx_icob_group_n#c.treat
lincom 2.xx_icob_group_n#c.treat - 4.xx_icob_group_n#c.treat

lincom 5.xx_icob_group_n#c.treat - 1.xx_icob_group_n#c.treat
lincom 5.xx_icob_group_n#c.treat - 4.xx_icob_group_n#c.treat


reghdfe n400 treat c.treat#ib3.xx_lpr_group_n ,  vce(robust) absorb(rblock)
estimates store model3


lincom 1.xx_lpr_group_n#c.treat - 4.xx_lpr_group_n#c.treat
lincom 2.xx_lpr_group_n#c.treat - 4.xx_lpr_group_n#c.treat


esttab model? using "hettable1.tex",          ///
scalars("N N") noobs       ///
b(2) se(2) sfmt(0 2)  noomitted  keep(treat *c.treat)                   ///
label interaction(" \$\cdot\$ ")           ///
mtitles("Class of Admission" "Region of Birth" "Years as LPR") ///
mgroups("N-400 Submission", pattern(1 0 0 ) ///
 prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}) )  ///
star(* 0.05) ///
booktabs                       ///
title("Effect Heterogeneity: Class of Admission, Region of Birth, and Years as LPR \label{tab:hetcoa}") replace ///
postfoot("\hline \multicolumn{4}{p{.5\textwidth}}{\scriptsize Regression coefficients with robust standard errors in parentheses. p $<$ 0.05. " ///
	"Reference categories: Humanitarian for Class of admission, China/India for Region of birth, and 8 for Years as LPR. " ///
	"All models include randomization strata fixed effects.}" ///
	"\end{tabular} \end{table}" ///
	)
	
	
* gender, age, marital status
estimates drop model*
reghdfe n400 treat c.treat#ib1.xx_igender_n ib1.xx_igender_n ,  vce(robust) absorb(rblock)
estimates store model1

reghdfe n400 treat c.treat#ib1.xx_age_cat_n ib1.xx_age_cat_n  ,  vce(robust) absorb(rblock)
estimates store model2

reghdfe n400 treat c.treat#ib3.xx_ims_n  ib3.xx_ims_n ,  vce(robust) absorb(rblock)
estimates store model3

esttab model? using "hettable2.tex",          ///
scalars("N N") noobs       ///
b(2) se(2) sfmt(0 2)  noomitted                     ///
label interaction(" \$\cdot\$ ")  kee(treat *c.treat)         ///
mtitles("Gender" "Age" "Marital Status") ///
mgroups("N-400 Submission", pattern(1 0 0 ) ///
 prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}) )  ///
star(* 0.05) ///
booktabs                       ///
title("Effect Heterogeneity: Gender, Age, and Marital Status \label{tab:hetgender}") replace ///
postfoot("\hline \multicolumn{4}{p{.5\textwidth}}{\scriptsize Regression coefficients with robust standard errors in parentheses. p $<$ 0.05. " ///
	"Reference categories: Female for Gender, 19-29 years old for Age, and Single for Marital Status. " ///
	"All models include randomization strata fixed effects and main effects for each moderator category (gender, age, and marital status).}" ///
	"\end{tabular} \end{table}" ///
	)


*************************************************
/* Robustness Checks: ITT */
*************************************************

rename xxx_prefiled prefiled
drop   xxx_medinc

cap estimates drop model*

reghdfe n400 treat ,  vce(robust) absorb(rblock)
estadd local droppre  "No"
estadd local ccon "No"
estadd local lin "No"
estimates store model1

reghdfe n400 treat ,  vce(robust) absorb(rblock) , if prefiled!=1
estadd local droppre  "Yes"
estadd local ccon "No"
estadd local lin "No"
estimates store model2

reghdfe n400 treat xxx_*,  vce(robust) absorb(rblock) 
estadd local droppre  "No"
estadd local ccon "Yes"
estadd local lin "No"
estimates store model3

reghdfe n400 treat xxx_*,  vce(robust) absorb(rblock) , if prefiled!=1
estadd local droppre  "Yes"
estadd local ccon "Yes"
estadd local lin "No"
estimates store model4

* lin estimator
cap drop iii* ppp* 
foreach x of varlist xxx_* {
	egen double meanvar = mean(`x')
	gen double iii`x' = `x' - meanvar
	drop meanvar
	gen double ppp`x' = treat * iii`x'
}

reghdfe n400 treat iii* ppp*,  vce(robust) absorb(rblock) 
estadd local droppre "No"
estadd local ccon "Yes"
estadd local lin "Yes"
estimates store model5

drop iii* ppp* 
foreach x of varlist xxx_* {
	egen double meanvar = mean(`x')   if prefiled!=1 
	gen double iii`x' = `x' - meanvar if prefiled!=1
	drop meanvar
	gen double ppp`x' = treat * iii`x' if prefiled!=1
}

reghdfe n400 treat iii* ppp*,  vce(robust) absorb(rblock) , if prefiled!=1
estadd local droppre  "Yes"
estadd local ccon "Yes"
estadd local lin "Yes"
estimates store model6


esttab model? using "robtable.tex",          ///
scalars("N N" "droppre Drop Pre-Submissions" "ccon Controls" "lin Lin-Estimator")        ///
b(2) se(2) sfmt(0 2)  noomitted  keep(treat)                   ///
star(* 0.05) label ///
mgroups("N-400 Submission", pattern(1 0 0 0 0 0 ) ///
 prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}) )  ///
nomtitles ///
booktabs                       ///
title("Robustness Checks: Effect of Letter on N-400 Submissions \label{tab:robcheck}") replace ///
postfoot("\hline \multicolumn{7}{p{.6\textwidth}}{\scriptsize Regression estimates shown of intention-to-treat effects with robust standard errors in parentheses. * p $<$ 0.05. Column (1) is based on a linear regression of N-400 submission on treatment assignment and randomization strata fixed effects. Columns (2), (4), and (6) exclude all observations of N-400 submissions made between the randomization date and the start of letter mailings on October 26, 2023; this exclusion was balanced between the treatment and control group. Columns (3), (5), and (6) add controls for marital status, age, gender, state, and propensity to refile. Columns (5) and (6) estimate treatment effects using the method from Lin (2013), which adds the control variables mean centered and interacted with treatment assignment.}" ///
	"\end{tabular} \end{table}" ///
	)

	
	
	
*************************************************
/* Robustness Checks: ITT and Standard errors */
*************************************************

cap estimates drop model*

reghdfe n400 treat ,  vce(robust) absorb(rblock)
estadd local droppre  "No"
estadd local ccon "No"
estadd local lin "No"
estadd local sea "None, Robust"
estimates store model1

reghdfe n400 treat ,  cluster(izip) absorb(rblock) 
estadd local droppre  "No"
estadd local ccon "No"
estadd local lin "No"
estadd local sea "Zip Code"
estimates store model2

reghdfe n400 treat xxx_*,  cluster(icbsa) absorb(rblock) 
estadd local droppre  "No"
estadd local ccon "No"
estadd local lin "No"
estadd local sea "CBSA"
estimates store model3

reghdfe n400 treat xxx_*,  cluster(rblock) absorb(rblock) 
estadd local droppre  "No"
estadd local ccon "No"
estadd local lin "No"
estadd local sea "Randomization Block"
estimates store model4

esttab model? using "robtable_se.tex",          ///
scalars("N N" "droppre Drop Pre-Submissions" "ccon Controls" "lin Lin-Estimator" "sea Cluster")        ///
b(2) se(2) sfmt(0 2)  noomitted  keep(treat)                   ///
star(* 0.05) label ///
mgroups("N-400 Submission", pattern(1 0 0 0) ///
 prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}) )  ///
nomtitles ///
booktabs                       ///
title("Robustness Checks: Effect of Letter on N-400 Submissions with Standard Error Adjustments \label{tab:robcheck_se}") replace ///
postfoot("\hline \multicolumn{5}{p{.6\textwidth}}{\scriptsize Regression estimates shown of intention-to-treat effects with robust standard errors in parentheses. * p $<$ 0.05. All columns are based on a linear regression of N-400 submission on treatment assignment and randomization strata fixed effects. Column (1) uses robust standard errors whereas Columns (2), (3), and (4) cluster by zip code, CBSA, and randomization block. .}" ///
	"\end{tabular} \end{table}" ///
)

	

	
	
	
*************************************************
/* Robustness Checks: LATE */
*************************************************


estimates drop model*

* first stage
reghdfe mdeliv treat , vce(robust) absorb(rblock) 
test treat==0
local  fstat: di %8.0g r(F)
estadd local fstat = `fstat'
estadd local controls = "No"
estadd local droppre    "No"
estimates store model1

* first stage: drop prefile
reghdfe mdeliv treat , vce(robust) absorb(rblock) , if prefiled!=1
test treat==0
local  fstat: di %8.0g r(F)
estadd local fstat = `fstat'
estadd local controls = "No"
estadd local droppre    "Yes"
estimates store model2


* first stage: add controls
reghdfe mdeliv treat xxx_*, vce(robust) absorb(rblock) 
test treat==0
local  fstat: di %8.0g r(F)
estadd local fstat = `fstat'
estadd local controls = "Yes"
estadd local droppre    "No"
estimates store model3

* first stage: both
reghdfe mdeliv treat xxx_*, vce(robust) absorb(rblock) , if prefiled!=1
test treat==0
local  fstat: di %8.0g r(F)
estadd local fstat = `fstat'
estadd local controls = "Yes"
estadd local droppre    "Yes"
estimates store model4

* late
eststo s5: ivregress 2sls n400 (mdeliv=treat) i.rblock , robust 
estadd local controls = "No"
estadd local droppre    "No"
estimates store model5


* late: drop prefile
ivregress 2sls n400 (mdeliv=treat) i.rblock , robust , if prefiled!=1
estadd local controls = "No"
estadd local droppre    "Yes"
estimates store model6


* late: add controls
ivregress 2sls n400 (mdeliv=treat) i.rblock xxx_*, robust 
estadd local controls = "Yes"
estadd local droppre    "No"
estimates store model7

* late: both
ivregress 2sls n400 (mdeliv=treat) i.rblock xxx_*, robust , if prefiled!=1
estadd local controls = "Yes"
estadd local droppre    "Yes"
estimates store model8

esttab model? using "IVtable.tex",          ///
scalars("N N" "fstat Partial F-Statistic" "droppre Drop Pre-Submissions" "controls Controls" )        ///
b(2) se(2) sfmt(0 2)  noomitted  keep(treat mdeliv)                   ///
star(* 0.05) label ///
nomtitles ///
mgroups("Letter Delivered" "N-400 Submission", pattern(1 0 0 0 1 0 0 0) ///
 prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}) )  ///
booktabs                       ///
title("Robustness Checks: Effect of Letter Delivered on N-400 Submissions \label{tab:IVcheck}") replace ///
postfoot("\hline \multicolumn{9}{p{.8\textwidth}}{\scriptsize Regression estimates shown with robust standard errors in parentheses. * p $<$ 0.05. Columns (1) to (4) show first stage estimates based on a linear regression of letter delivery on treatment assignment and randomization strata fixed effects. Columns (5) to (8) show local average treatment effect (LATE) estimates based on instrumental variable regressions of N-400 submissions on letter delivery, instrumented by treatment assignment, and randomization strata fixed effects. Columns (2), (4), (6) and (8) exclude all observations of N-400 submissions made between the randomization date and the start of letter mailings on October 26, 2023; this exclusion was balanced between the treatment and control group. Columns (3), (4), (7) and (8) add controls for marital status, age, gender, state, and propensity to refile.}" ///
	"\end{tabular} \end{table}" ///
	)



*************************************************
/* Robustness Checks: Concentration */
*************************************************

cap estimates drop model*

label var treatc "Treatment Group * Concentration of LPRs"

* By concentration in the zip 
reghdfe n400 treat if concentration_per_zip!=.&pop!=.&medinc!=.,  vce(robust) absorb(rblock)
estadd local basic  "No"
estadd local linear "--"
estimates store model1

reghdfe n400 treat concentration_per_zip if concentration_per_zip!=.&pop!=.&medinc!=.,  vce(robust) absorb(rblock)
estadd local basic  "Yes"
estadd local linear "--"
estimates store model2

reghdfe n400 treat treatc concentration_per_zip if concentration_per_zip!=.&pop!=.&medinc!=.,  vce(robust) absorb(rblock)
estadd local basic  "Yes"
estadd local linear "--"
estimates store model3

reghdfe n400 treat treatc concentration_per_zip pop medinc if concentration_per_zip!=.&pop!=.&medinc!=. ,  vce(robust) absorb(rblock)
estadd local basic  "Yes"
estadd local linear "Linear"
estimates store model4

reghdfe n400 treat treatc concentration_per_zip pop medinc pop2 pop3 medinc2 medinc3 if concentration_per_zip!=.&pop!=.&medinc!=.,  vce(robust) absorb(rblock)
estadd local basic  "Yes"
estadd local linear "Cubic"
estimates store model5


esttab model? using "concentration.tex",          ///
scalars("N N" "basic Concentration control" "linear Income and Population controls")        ///
b(2) se(2) sfmt(0 2)  noomitted  keep(treat treatc)                   ///
star(* 0.05) label ///
mgroups("N-400 Submission", pattern(1 0 0 0 0 ) ///
 prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}) )  ///
nomtitles ///
booktabs                       ///
title("Heterogeneity: Effect of Letter on N-400 Submissions by Concentration of LPRs per Zip Code") replace ///
postfoot("\hline \multicolumn{6}{p{.6\textwidth}}{\scriptsize Regression estimates shown of intention-to-treat effects with robust standard errors in parentheses. * p $<$ 0.05. Column (1) is based on a linear regression of N-400 submission on treatment assignment and randomization strata fixed effects after removing individuals with missing zip code data. Columns (2) through (5) include a control for the concentration of LPRs, defined as the number of LPRs in our analytic sample divided by the 2023 ACS 5-Year Zip Code population estimates. Column (3) interacts concentration with our treatment indicator, Column (4) adds linear controls for zip code population size and median income, and column (5) includes cubic controls for zip code population size and median income.}" ///
	"\end{tabular} \end{table}" ///
	)

* Predictions
su concentration_per_zip if tag1==1, d
