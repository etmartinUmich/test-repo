****************************************************************************************************************
*This was used for the case crossover comparison in the Natural History manuscript
*Comparison windows were changed to 7d for the final manuscript in response to 
*reviewer comments on a previous submission.
****************************************************************************************************************

/*These numbers changed after I realized that patients with all 
missing were reported as having a sum of 0 in SalivaAnalysis code. This was fixed
by using total() instead of sum()*/


clear
insheet using "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\CC6Analysis_7d_30SEP.raw", clear
describe
reshape long case6 cc6pcpill cc6fever cc6cough cc6runnynose cc6vomit cc6diarrhea cc6generalrash cc6localrash cc6diaperrash cc6fussy cc6seizure cc6roseolla cc6sxbi_all cc6acute cc6newcough cc6newrunnynose cc6resp cc6gi, i(patient primary6date) j(tag)

foreach x in primary6date {
gen `x'_f = date(`x', "DMY")
drop `x'
gen `x' = `x'_f
format `x' %d
drop `x'_f
describe `x'
}

*regular cough and runny nose excluded below. Only newcough and newrunnynose used throughout
*Data presented by case=1 (table column 1) or case=0 (table column 2)
sort case6
by case6: tab1 cc6pcpill cc6fever cc6vomit cc6diarrhea cc6generalrash cc6localrash /// use this to check missing data elements
cc6fussy cc6sxbi_all cc6acute cc6resp cc6gi cc6newcough cc6newrunnynose, miss
by case6: tab1 cc6pcpill cc6fever cc6vomit cc6diarrhea cc6generalrash cc6localrash /// use this to get valid percent
cc6fussy cc6sxbi_all cc6acute cc6resp cc6gi cc6newcough cc6newrunnynose		 ///**********cc6resp is "New Respiratory Symptoms" on Table 3**********

*MANUSCRIPT REPORTS CLOGIT RESULTS FOR ALL BELOW
foreach x in cc6pcpill cc6fever { //regular cough and runnynose not in manuscript, limited to new onset below
clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}


foreach x in cc6vomit cc6diarrhea cc6generalrash cc6localrash {
clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}

foreach x in cc6fussy cc6sxbi_all cc6acute cc6resp cc6gi {
clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}

foreach x in cc6newcough cc6newrunnynose {
clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}


*Do a sensitivity analysis for season
gen primarymonth= month(primary6date)
tab primarymonth

gen respseas=1
  replace respseas=0 if primarymonth>=3 & primarymonth<=10
  
  tab respseas primarymonth, miss
  
foreach x in cc6pcpill cc6fever { //regular cough and runnynose not in manuscript, limited to new onset below
bysort respseas: clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}

gen interact= respseas*case6
clogit cc6pcpill case6 respseas interact, group(patient) or


foreach x in cc6diarrhea cc6generalrash {
bysort respseas: clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}

clogit cc6vomit case6 if respseas==1, group(patient) or
*clogit cc6vomit case6 if respseas==0, group(patient) or sample size too small
*clogit cc6localrash case6 if respseas==1, group(patient) or sample size too small
*clogit cc6localrash case6 if respseas==0, group(patient) or sample size too small


foreach x in cc6fussy cc6sxbi_all cc6acute cc6resp cc6gi {
bysort respseas: clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}

foreach x in cc6newcough cc6newrunnynose {
bysort respseas: clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}


*Do a sensitivity analysis about quantity
clear
insheet using "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\CC6AnalysisQuant_7d_30SEP.raw", clear
describe
reshape long case6 cc6pcpill cc6fever cc6cough cc6runnynose cc6vomit cc6diarrhea cc6generalrash cc6localrash cc6diaperrash cc6fussy cc6seizure cc6roseolla cc6sxbi_all cc6acute cc6newcough cc6newrunnynose cc6resp cc6gi, i(patient maxlog) j(tag)

describe
tab case6
summ maxlog if case6==1, detail
*66 Cases

*Select and control periods for sensitivity analysis, based on having a log viral load above 4
gen taghigh=0 if case6==1
  replace taghigh=1 if case6==1 & maxlog>=6 & maxlog!=.

bysort patient: egen highbov=max(taghigh)
tab highbov case6

*36 cases and 72 control periods used

sort case6 highbov
by case6: tab1 cc6pcpill cc6fever cc6vomit cc6diarrhea cc6generalrash cc6localrash /// use this to check missing data elements
cc6fussy cc6sxbi_all cc6acute cc6resp cc6gi cc6newcough cc6newrunnynose if highbov==1, miss
by case6: tab1 cc6pcpill cc6fever cc6vomit cc6diarrhea cc6generalrash cc6localrash /// use this to get valid percent
cc6fussy cc6sxbi_all cc6acute cc6resp cc6gi cc6newcough cc6newrunnynose if highbov==1

foreach x in cc6pcpill cc6fever {
tab `x' case6 if highbov==1, col
clogit `x' case6  if highbov==1, group(patient) or
*xtgee `x' case6 if highbov==1, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}


foreach x in cc6vomit cc6generalrash cc6diarrhea cc6localrash {
tab `x' case6 if highbov==1, col
clogit `x' case6 if highbov==1, group(patient) or
*xtgee `x' case6 if highbov==1, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}

foreach x in cc6fussy cc6sxbi_all cc6acute cc6resp cc6gi {
tab `x' case6 if highbov==1, col
clogit `x' case6 if highbov==1, group(patient) or
*xtgee `x' case6 if highbov==1, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}

foreach x in cc6newcough cc6newrunnynose cc6cough {
tab `x' case6 if highbov==1, col
clogit `x' case6 if highbov==1, group(patient) or
*xtgee `x' case6 if highbov==1, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}

*Extra - not in text at this point
foreach x in cc6cough cc6runnynose {
tab `x' case6 if highbov==1, col
clogit `x' case6 if highbov==1, group(patient) or
xtgee `x' case6 if highbov==1, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}




**************
*Do a recurrent case crossover analysis
****************
clear
insheet using "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\CC6AnalysiswRecur_7d_30SEP.raw"
reshape long case6 cc6pcpill cc6fever cc6cough cc6runnynose cc6vomit cc6diarrhea cc6generalrash cc6localrash cc6diaperrash cc6fussy cc6seizure cc6roseolla cc6sxbi_all cc6acute cc6newcough cc6newrunnynose cc6resp cc6gi rr6pcpill rr6fever rr6cough rr6runnynose rr6vomit rr6diarrhea rr6generalrash rr6localrash rr6diaperrash rr6fussy rr6seizure rr6roseolla rr6sxbi_all rr6acute rr6newcough rr6newrunnynose rr6resp rr6gi, i(patient) j(tag)

*New variable to define interim between end of primary and start of next event
gen interim=recur6date-eventenddt
summ interim, detail

drop if recur6date==.


sort case6
by case6: tab1 rr6pcpill rr6fever rr6vomit rr6diarrhea rr6generalrash rr6localrash /// use this to check missing data elements
rr6fussy rr6sxbi_all rr6acute rr6resp rr6gi rr6newcough rr6newrunnynose, miss
by case6: tab1 rr6pcpill rr6fever rr6vomit rr6diarrhea rr6generalrash rr6localrash /// use this to get valid percent
rr6fussy rr6sxbi_all rr6acute rr6resp rr6gi rr6newcough rr6newrunnynose

foreach x in rr6pcpill rr6fever  { //regular cough and runnynose not in manuscript, limited to new onset below
clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}


foreach x in rr6vomit rr6diarrhea rr6generalrash rr6localrash {
clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
} 

foreach x in rr6fussy rr6sxbi_all rr6acute rr6resp rr6gi {
clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}

foreach x in rr6newcough rr6newrunnynose {
clogit `x' case6, group(patient) or
*xtgee `x' case6, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}
