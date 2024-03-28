**************************************************************************************
*This code contains the main analysis for the HBOV Natural History Manuscript
*Comparison windows were changed to 7d for the final manuscript in response to 
*reviewer comments on a previous submission.
************************************************************************************

/* This modification is to redefine symptom windows to only capture symptoms 
occuring 7 days before the first infection date (rather than 14 days)
*/

clear
clear matrix

use "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\GrantRes1_SampleSx_071313.dta", clear


*Mark earliest entry with oneper==1
gen oneper=0
replace oneper=1 if start==date //start is the first available bovpos test date.  Using oneper will select kids that had bov testing performed
count if oneper==1

/*Optional code for if I want to select out all survey data with no bov testing performed)*/
***********************************************
*Tag kids with survey data but no bov testing
**************************************************
sort patient date
by patient: egen include=max(oneper)


***********************************************
*Eliminate 5 children with spotty testing:
***********************************************

foreach x in 10 20 23 26 32 {
drop if patient==`x'
}
count if oneper==1

*How many children had BoV detected at any point
tab tally if oneper

*Show what the range of followup time is in this group
summ totaltime if oneper==1, detail

*****************************************************************
*Recode event to require 2 positives to be an event
*****************************************************************
gen event6=event
  replace event6=0 if eventdur<6 & eventdur!=.
  
gen event6start=eventstart if event6==1
gen event6startdt=eventstartdt if event6==1

sort patient event6start date
by patient event6start: egen event6seq=seq(), from(1)
replace event6seq=. if event6start!=1 
sort patient date
*list patient date event event6 eventstart event6start eventstartdt event6startdt eventseq event6seq if patient==334
  
sort patient eventseq date
replace event6seq = event6seq[_n-1] if eventseq==eventseq[_n-1] & event6==1

sort patient date
by patient: egen tally6=total(event6start)
gen anytally6=0
  replace anytally6=1 if tally6>=1 & tally6!=.

tab eventdur if eventstart==1 & eventseq==1
****REMEMBER - the primary Event6 group is not simply a subset of the primary event group - there were kids who had an eligible
*event6 primary even after the short events are removed
*Using new longer definition, 66 kids had primary events
tab eventdur if event6start==1 & event6seq==1

**************************************************
*DEMOGRAPHICS AND STUDY DESCRIPTIVES
*************************************************
*Number of kids
tab oneper

*Show what the range of followup time is in this group
summ totaltime if oneper==1, detail


*How many children had BoV detected at any point
tab tally if oneper

tab gender anytally6 if oneper, col
tab race race_text
tab race_text anytally6 if oneper
tab race anytally6 if oneper, col
tab breastfed anytally6 if oneper, col
    gen ageatbfstop=untilwhen-infantdob
	summ ageatbfstop if oneper, detail
    bysort anytally6: summ ageatbfstop if oneper, detail

tab daycare anytally6 if oneper, col
    gen ageatdcstart=daycarestart-infantdob
	summ ageatdcstart if oneper, detail
	tab daycareendc daycare if oneper, miss
	gen ageatdcend=daycareend-infantdob if daycareendc==1
	summ ageatdcend if oneper, detail
	
	gen dcfulltime=0 if daysweek!=""
	  replace dcfulltime=1 if daysweek=="1"|daysweek=="1-2"|daysweek=="2"|daysweek=="2-3"|daysweek=="2.5"
	  replace dcfulltime=2 if daysweek=="3"|daysweek=="3-4"|daysweek=="3-5"|daysweek=="4"|daysweek=="4-5"
	  tab daysweek dcfulltime if oneper, miss
	tab dcfulltime if oneper

tab playgroup anytally6 if oneper & playgroup!="unknown", col
    tab playgroup daycare if oneper, miss
	gen ageatpgstart=playstart-infantdob
	summ ageatpgstart if oneper, detail
	tab playendc playgroup if oneper, miss
	gen ageatpgend=playend-infantdob if playendc==1
	summ ageatpgend if oneper, detail
	
gen twin=1 if comment=="twin 223"|comment=="twin 224"
*there is only one twin pair in the selected set

tab numberoffamily18yrs
    gen anysibs=0 if numberoffamily18yrs!=.
	replace anysibs=1 if numberoffamily18yrs==1|numberoffamily18yrs==2
	tab numberoffamily18yrs anysibs
	tab numberoffamily18yrs if oneper
	tab anysibs anytally6 if oneper, col

gen maternalage=(infantdob-maternaldob)/365.25
summ maternalage if oneper, detail
bysort anytally6: summ maternalage if oneper, detail

gen matagecat=0
  replace matagecat=1 if maternalage>30 & maternalage<40
  replace matagecat=2 if maternalage>40 & maternalage!=.
  
gen yearenroll=year(infantdob)
tab yearenroll

****************************************************************
*Graph and describe first detections by age
*Use months in study for proxy right now
****************************************************************
gen smo_int=int(smo)

foreach x of numlist 0/26 {
  gen primary_smo_`x'=1 if smo_int==`x' & eventstart==1 & eventseq==1
  gen primary6_smo_`x'=1 if smo_int==`x' & event6start==1 & event6seq==1
}

*Descriptive statistics for age at primary (approximated by study dates)
summ smo if eventstart==1 & eventseq==1, detail
summ smo if event6start==1 & event6seq==1, detail
gen agedays=date-infantdob
summ agedays if eventstart==1 & eventseq==1, detail

*In this pilot sample, we only have full data for all 50 up to 18 months.  Only graph that far.
*graph bar (sum) primary_smo_0 - primary_smo_26, showyvars yvaroptions(relabel(1 "0" 2 "1" 3 "2" 4 "3" 5 "4" 6 "5" 7 "6" 8 "7" 9 "8" 10 "9" 11 "10" 12 "11" 13 "12" 14 "13" 15 "14" 16 "15" 17 "16" 18 "17" 19 "18" 20 "19" 21 "20" 22 "21" 23 "22" 24 "23" 25 "24" 26 "25" 27 "26") gap(5) label(angle(vertical))) legend(off)
  *graph save "O:\Analyses\Zerr Saliva\StataGraphs\AgeAtPrimaryStart.gph", replace

*graph bar (sum) primary6_smo_*, showyvars yvaroptions(relabel(1 "0" 2 "1" 3 "2" 4 "3" 5 "4" 6 "5" 7 "6" 8 "7" 9 "8" 10 "9" 11 "10" 12 "11" 13 "12" 14 "13" 15 "14" 16 "15" 17 "16" 18 "17" 19 "18" 20 "19" 21 "20" 22 "21" 23 "22" 24 "23" 25 "24" 26 "25" 27 "26") gap(5) label(angle(vertical))) legend(off)
*****************************************************
*AIM 1: Symptoms
****************************************************

sort patient date
by patient: gen coughlag1=cough[_n-1]
by patient: gen coughlag2=cough[_n-2]
by patient: gen coughlag3=cough[_n-3]
by patient: gen coughlag4=cough[_n-4]
by patient: gen coughlag5=cough[_n-5]
by patient: gen coughlag6=cough[_n-6]
by patient: gen coughlag7=cough[_n-7]

egen prevcough = rowtotal(coughlag1 coughlag2 coughlag3 coughlag4 coughlag5 coughlag6 coughlag7), missing
tab prevcough cough, miss
gen newcough=0 if cough!=.
  replace newcough=1 if prevcough==0 & cough==1
  
sort patient date
by patient: gen runnynoselag1=runnynose[_n-1]
by patient: gen runnynoselag2=runnynose[_n-2]
by patient: gen runnynoselag3=runnynose[_n-3]
by patient: gen runnynoselag4=runnynose[_n-4]
by patient: gen runnynoselag5=runnynose[_n-5]
by patient: gen runnynoselag6=runnynose[_n-6]
by patient: gen runnynoselag7=runnynose[_n-7]

egen prevrunnynose = rowtotal(runnynoselag1 runnynoselag2 runnynoselag3 runnynoselag4 runnynoselag5 runnynoselag6 runnynoselag7), missing
gen newrunnynose=0 if runnynose!=.
  replace newrunnynose=1 if prevrunnynose==0 & runnynose==1


*********Create an aggregate symptom variable
egen sx_all=rowtotal(fever cough runnynose vomit diarrhea generalrash localrash diaperrash fussy seizure roseolla pcpill), miss
*Create a binary variable for any symptoms
gen sxbi_all=0 if sx_all==0
  replace sxbi_all=1 if (sx_all==1 | sx_all>1) & sx_all!=.
*Create a binary variable for acute illness
egen acute=rowtotal(fever newcough vomit diarrhea generalrash seizure pcpill), missing
  replace acute=1 if acute>1 & acute !=.
*Create an aggregate resp variable
egen resp=rowtotal(newcough newrunnynose), missing
  replace resp=1 if resp==2
egen gi=rowtotal(vomit diarrhea), missing
  replace gi=1 if gi==2
 

 

********************************************************
*Symptoms - Define symptom windows
********************************************************
*Recode pcpill to be non-missing if other symptom data is present
replace pcpill=0 if fever!=. & pcpill==.

gen primary6date = date if bovpos==1 & event6start==1 & event6seq==1
bysort patient: egen fillp6date=max(primary6date)
replace primary6date=fillp6date
drop fillp6date
tab primary6date if oneper, miss

tab fever, miss

foreach x in pcpill fever cough runnynose vomit diarrhea generalrash localrash diaperrash fussy seizure roseolla sxbi_all acute resp gi newcough newrunnynose {
  gen count1`x'=`x' if date>=(primary6date-7) & date<=(primary6date+7)
  gen primary`x'dt=date if count1`x'==1 
  bysort patient: egen sum1`x'=total(count1`x'), missing
  gen primary6`x'=sum1`x' //set primary variable to equal tally of symptom days in primary window
  replace primary6`x'=1 if sum1`x'>0 & sum1`x'!=. //set primary variable to a yes/no for presence of symptom
  replace primary6`x'=. if primary6date==.

  gen cc6`x'1=primary6`x'
 
  gen count2`x'=`x' if date>=(primary6date+28) & date<=(primary6date+42)
  gen cc6`x'dt2=date if count2`x'==1
  bysort patient: egen sum2`x'=total(count2`x'), missing
  gen cc6`x'2=sum2`x'
  replace cc6`x'2=1 if sum2`x'>0 & sum2`x'!=.
  replace cc6`x'2=. if primary6date==.

  gen count3`x'=`x' if date>=(primary6date-42) & date<=(primary6date-28)
  gen cc6`x'dt3=date if count3`x'==1
  bysort patient: egen sum3`x'=total(count3`x'), miss
  gen cc6`x'3=sum3`x'
  replace cc6`x'3=1 if sum3`x'>0 & sum3`x'!=.
  replace cc6`x'3=. if primary6date==.
}

*G variables are used to create extended graph in "SymptomWindowstest_06JAN14.do" code
foreach x in pcpill fever cough runnynose vomit diarrhea generalrash localrash diaperrash fussy seizure roseolla sxbi_all acute resp gi newcough newrunnynose {
  gen gcount1`x'=`x' if date>=(primary6date-28) & date<=(primary6date+21)
  gen gprimary`x'dt=date if gcount1`x'==1 
  bysort patient: egen gsum1`x'=total(gcount1`x'), missing
  gen gprimary6`x'=gsum1`x' //set primary variable to equal tally of symptom days in primary window
  replace gprimary6`x'=1 if gsum1`x'>0 & gsum1`x'!=. //set primary variable to a yes/no for presence of symptom
  replace gprimary6`x'=. if primary6date==.
}

tab reasonfordrvisit if count1pcpill==1 //describe reasons for PCP visit

sort patient date
format date %d
format primary6date %d

list patient date pcpill fever bovpos primary6date count1pcpill count1fever count2pcpill count2fever count3pcpill count3fever if patient==155
 
gen case61=1
gen case62=0
gen case63=0
*OFFICIAL DATA FOR SYMPTOM TABLE, COLUMNS 1 through 3:
outsheet patient primary6date case61 case62 case63 cc6pcpill1-cc6newrunnynose3 using "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\CC6Analysis_7d_30SEP.raw" if event6start==1 & event6seq==1, replace


************************************
*Symptom Table 
*************************************
*These should match n of primary events (column one of table)
*cough and runnynose removed - only want new onset
tab primary6pcpill if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6pcpill if oneper==1 & tally>0, miss

tab primary6fever if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6fever if oneper==1 & tally>0, miss

tab primary6newcough if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6newcough if oneper==1 & tally>0, miss

tab primary6newrunnynose if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6newrunnynose if oneper==1 & tally>0, miss

tab primary6vomit if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6vomit if oneper==1 & tally>0, miss

tab primary6diarrhea if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6diarrhea if oneper==1 & tally>0, miss

tab primary6generalrash if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6generalrash if oneper==1 & tally>0, miss

tab primary6localrash if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6localrash if oneper==1 & tally>0, miss

tab primary6diaperrash if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6diaperrash if oneper==1 & tally>0, miss
	
tab primary6fussy if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6fussy if oneper==1 & tally>0, miss

tab primary6seizure if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6seizure if oneper==1 & tally>0, miss

tab primary6roseolla if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6roseolla if oneper==1 & tally>0, miss

*This is NOT in table
tab primary6sxbi_all if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6sxbi_all if oneper==1 & tally>0, miss

*Paper uses primary6acute variable for "any symptom", not in table
tab primary6acute if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6acute if oneper==1 & tally>0, miss

*This is in table
tab primary6resp if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6resp if oneper==1 & tally>0, miss

tab primary6gi if bovpos==1 & event6start==1 & event6seq==1, miss
tab1 primary6gi if oneper==1 & tally>0, miss

********************************************************
*Recurrent Symptoms - Define symptom windows in first recurrence
********************************************************

gen recur6date = date if bovpos==1 & event6start==1 & event6seq==2
bysort patient: egen fillr6date=max(recur6date)
replace recur6date=fillr6date
drop fillr6date
tab recur6date if oneper, miss

*Get recurrence definition
tab event6seq if event6start==1
bysort event6seq: summ eventdur if event6start==1  
gen recurrent6=0 if event6seq==1  //variable used to compare primary vs all recurrent
  replace recurrent6=1 if event6seq>1 & event6seq!=.

bysort recurrent6: summ eventdur if event6start==1, detail  
ranksum eventdur if event6start==1, by(recurrent6)

foreach x in pcpill fever cough runnynose vomit diarrhea generalrash localrash diaperrash fussy seizure roseolla sxbi_all acute resp gi newcough newrunnynose {
  gen countr1`x'=`x' if date>=(recur6date-7) & date<=(recur6date+7)
  gen recur`x'dt=date if countr1`x'==1 
  bysort patient: egen sumr1`x'=total(count1`x'), missing
  gen recur6`x'=sumr1`x' //set primary variable to equal tally of symptom days in primary window
  replace recur6`x'=1 if sumr1`x'>0 & sumr1`x'!=. //set primary variable to a yes/no for presence of symptom
  replace recur6`x'=. if recur6date==.

  gen rr6`x'1=recur6`x'
 
  gen countr2`x'=`x' if date>=(recur6date+28) & date<=(recur6date+42)
  gen rr6`x'dt2=date if countr2`x'==1
  bysort patient: egen sumr2`x'=total(countr2`x'), missing
  gen rr6`x'2=sumr2`x'
  replace rr6`x'2=1 if sumr2`x'>0 & sumr2`x'!=.
  replace rr6`x'2=. if recur6date==.

  gen countr3`x'=`x' if date>=(recur6date-42) & date<=(recur6date-28)
  gen rr6`x'dt3=date if countr3`x'==1
  bysort patient: egen sumr3`x'=total(countr3`x'), miss
  gen rr6`x'3=sumr3`x'
  replace rr6`x'3=1 if sumr3`x'>0 & sumr3`x'!=.
  replace rr6`x'3=. if recur6date==.
}

 
outsheet patient case61 case62 case63 cc6pcpill1-cc6gi3 rr6pcpill1-rr6newrunnynose3 recur6date event6start eventenddt using "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\CC6AnalysiswRecur_7d_30SEP.raw" if event6start==1 & recurrent6==0, replace

*******************************
*SYMPTOM TABLE COMPARE RECURRENT EVENT TO PRIMARY EVENT
******************************
gen recurcompare=0 if event6start==1 & event6seq==1
 replace recurcompare=1 if event6start==1 & event6seq==2

foreach x in pcpill fever newcough newrunnynose vomit diarrhea generalrash localrash diaperrash fussy seizure roseolla sxbi_all acute resp gi {
gen `x'compare=cc6`x'1 if recurcompare==0
  replace `x'compare=rr6`x'1 if recurcompare==1 
}

foreach x in pcpill fever newcough newrunnynose vomit diarrhea generalrash localrash diaperrash fussy resp gi {
tab `x'compare recurcompare, miss
tab `x'compare recurcompare
}
/*THIS ANALYSIS IS PROBLEMATIC clogit doesn't work and gee doesn't converse for almost anything
*

foreach x in pcpill newcough vomit generalrash fussy{
tab `x'compare recurcompare, miss
tab `x'compare recurcompare
xtgee `x'compare recurcompare, i(patient) corr(exchangeable) family(binomial) link(logit) robust eform
}
*/

************************************************************
*Code "asymptomatic"
************************************************************
*sxbi_all defined above as an aggregate symptom variable

sort patient date
gen bovsx=0 if bovpos!=.
by patient: replace bovsx=1 if sxbi_all[_n-1]==1 & bovpos!=.
by patient: replace bovsx=1 if sxbi_all[_n+1]==1 & bovpos!=.
by patient: replace bovsx=. if sxbi_all[_n-1]==. & sxbi_all[_n+1]==. & bovpos!=.
replace bovsx=. if bovpos==.

list patient date bovpos bovsx sxbi_all if include==1 & patient==345

*Use this number for the % of bocavirus detections with no symptom recorded, manuscript text
tab bovpos bovsx if include==1, miss



*********
*Make a dataset of symptom dates in the primary windoww
save "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\BoVPrimaryEvents_100714.dta", replace

********************************************************
*AIM 2: Duration of Shedding and patterns of virus quantity
**********************************************************

********************************************************************
*Graph distribution of event length
********************************************************************
label var eventdur "Duration of Continuous Bocavirus Detection (in days)"
tab eventdur if event6start==1
histogram eventdur if event6start==1, freq bin(16)
summ eventdur if event6start==1, detail


gen eventdurcat=.
  replace eventdurcat=0 if eventdur==0
  replace eventdurcat=1 if eventdur>0 & eventdur<=30
  replace eventdurcat=2 if eventdur>30 & eventdur<=60
  replace eventdurcat=3 if eventdur>60 & eventdur<=90
  replace eventdurcat=4 if eventdur>90 & eventdur<=180
  replace eventdurcat=5 if eventdur>180 & eventdur!=.
label define labdurcat 0 "Single Pos" 1 "7 to 30 days" 2 "31 to 60 days" 3 "61 to 90 days" 4 "91 to 180 days" 5 "181+ days"
label values eventdurcat labdurcat
label var eventdur "Duration of Continuous Bocavirus Detection (in days)"
/*
*Generate dummy variables for bar graph
foreach x of numlist 0/5 {
gen eventdurcat`x'=1 if eventdurcat==`x'
}
label var eventdurcat0 "Single Pos"
label var eventdurcat1 "7 to 30 days"
label var eventdurcat2 "31 to 60 days"
label var eventdurcat3 "61 to 90 days"
label var eventdurcat4 "91 to 180 days"
label var eventdurcat5 "180+ days"

graph bar (sum) eventdurcat0 eventdurcat1 eventdurcat2 eventdurcat3 eventdurcat4 eventdurcat5 if eventstart==1, showyvars yvaroptions(relabel(1 "Single Pos" 2 "7 to 30 days" 3 "31 to 60 days" 4 "61 to 90 days" 5 "91 to 180 days" 6 "181+ days") gap(5) label(angle(vertical))) blabel(total) legend(off)
  *graph save "O:\Analyses\Zerr Saliva\StataGraphs\EventDuration_All.gph", replace

*How many events per kid in the above graph?
bysort patient: egen totalevents=sum(eventstart)
tab totalevents if oneper==1

tab eventdur if eventstart==1
histogram eventdur if eventstart==1, freq bin(16)
summ eventdur if eventstart==1, detail

graph bar (sum) eventdurcat0 eventdurcat1 eventdurcat2 eventdurcat3 eventdurcat4 eventdurcat5 if eventstart==1, showyvars yvaroptions(relabel(1 "Single Pos" 2 "7 to 30 days" 3 "31 to 60 days" 4 "61 to 90 days" 5 "91 to 180 days" 6 "181+ days") gap(5) label(angle(vertical))) blabel(total) legend(off)
  *graph save "O:\Analyses\Zerr Saliva\StataGraphs\EventDuration_All.gph", replace

*How many events total

*How many events per kid in the above graph?
tab totalevents if oneper==1

*Figure eventdur by month for heterogeneity grant
gen mdur=round(eventdur/30, 1)
*/
*Duration of shedding for primary only
tab eventdur if  event6start==1 & event6seq==1
*histogram eventdur if eventseq==1 & eventstart==1, freq bin(16)
summ eventdur if event6start==1 & event6seq==1, detail
tab eventdurcat if event6start==1 & event6seq==1

gen logq=log10(q_per_ml)

sort id event6seq
egen maxlog=max(logq), by(id event6seq)
sum logq, detail
sum maxlog if event6start==1, detail
bysort recurrent6: sum maxlog if event6start==1, detail
ranksum maxlog if event6start==1, by (recurrent6)

*Code virus quant categories:
gen maxcat=0 if maxlog>0 & maxlog<=4 & maxlog!=.
replace maxcat=1 if maxlog>4 & maxlog<=6 & maxlog!=.
replace maxcat=2 if maxlog>6 & maxlog<=8 & maxlog!=.
replace maxcat=3 if maxlog>8 & maxlog!=.

gen samdt=date-infantdob
glm logq samdt id if bovpos==1
*This looks at peak overall by age
glm maxlog samdt id if bovpos==1 & event6start==1

*Look at peak shedding by age, for primary events only

*Make age categorical 
gen agecat=0 if smo<=6 & smo!=.
replace agecat=1 if smo>6 & smo<=12
replace agecat=2 if smo>12 & smo<=18
replace agecat=3 if smo>18 & smo!=.

twoway (scatter maxlog smo) if event6start==1 & recurrent6==0

graph box maxlog if event6start==1 & recurrent6==0, over(agecat) 
reg maxlog smo if event6start==1 & recurrent6==0

tab agecat if event6start==1 & recurrent6==0

xi: reg maxlog i.agecat if event6start==1 & recurrent6==0
anova maxlog agecat if event6start==1 & recurrent6==0
bysort agecat: summ maxlog if event6start==1 & recurrent6==0, detail

summ maxlog if event6start==1 & recurrent6==0

*Repeat for year of birth
gen yearn=0 if yearenroll==1998
  replace yearn=1 if yearenroll==1999
  replace yearn=2 if yearenroll==2000
  replace yearn=3 if yearenroll==2001

twoway (scatter maxlog yearn) if event6start==1 & recurrent6==0
xi: reg maxlog i.yearenroll if event6start==1 & recurrent6==0
anova maxlog yearenroll if event6start==1 & recurrent6==0
bysort yearenroll: summ maxlog if event6start==1 & recurrent6==0, detail



*Get viral load for singletons
summ maxlog if eventdur==0, detail
 
***************************
*SYMPTOM TABLE COLUMNS 4-6 - Quant sensitivity analysis
***************************
outsheet patient maxlog case61 case62 case63 cc6pcpill1-cc6newrunnynose3 using "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\CC6AnalysisQuant_7d_30SEP.raw" if event6start==1 & recurrent6==0, replace




********************************************************
*Allele Switch Symptoms - Define symptom windows in first allele switch
********************************************************

gen alleledate = date if bovpos==1 & switchstart==1 & alleleorder==1
count if alleledate !=.
bysort patient: egen fillalleledate=max(alleledate)
replace alleledate=fillalleledate
drop fillalleledate
tab alleledate if oneper, miss

list patient alleleorder if bovpos==1 & alleleorder!=. 

  gen newid="Child A" if patient==91
replace newid="Child B" if patient==98
replace newid="Child C" if patient==107
replace newid="Child D" if patient==151
replace newid="Child E" if patient==189
replace newid="Child F" if patient==202
replace newid="Child G" if patient==221
replace newid="Child H" if patient==227
replace newid="Child I" if patient==243
replace newid="Child J" if patient==284
replace newid="Child K" if patient==317
replace newid="Child L" if patient==354

sort newid date
list newid date smo allele* if bovpos==1 & newid!="" & alleleorder!=.
outsheet newid date smo allele2842 allele3923 allele4181 allele4229 allele4742 allele4799 allele4862 allele4874 alleleorder alleledate using "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Table4.txt" if bovpos==1 & newid!="" & alleleorder!=., replace

*alleleorder=0 is original type. alleleorder==1 is first switch. alleleorder==2 is second switch
twoway (scatter graphlogbov smo if bovpos==0, yaxis(1) sort mcolor(gray) msize(medium) msymbol(circle_hollow) mcolor(gray)) ///
    (scatter graphlogbov smo if bovpos==1 & alleleorder==0, yaxis(1) sort mcolor(red) msize(medium) msymbol(circle)) ///
	(scatter graphlogbov smo if bovpos==1 & alleleorder==1, yaxis(1) sort mcolor(blue) msize(small) msymbol(square)) ///
	(scatter graphlogbov smo if bovpos==1 & alleleorder==2, yaxis(1) sort mcolor(green) msize(medium) msymbol(triangle)), ///
   by(newid, cols(2) iscale(0.6) imargin(small)) ///
  xscale(r(0 25)) ///
  yscale(r(1 9.5) axis(1)) ///
  ylabel(0(4)9.5, angle(horizontal) axis(1) labsize(small)) ///
  xtitle("Months", size(small)) ///
  ytitle("Log HBoV-1 Viral Load / mL Saliva", size(small)) ///
  xlabel(0(2)25) ///
  ysize(6) ///
  xsize(6) ///
  legend(order(1 "HBoV1(-)" 2 "1st HBoV1 Allele" 3 "2nd HBoV1 Allele" 4 "3rd HBoV1 Allele") colgap(*.5) rows(1) size(vsmall))
  drop newid

foreach x in pcpill fever cough runnynose vomit diarrhea generalrash localrash diaperrash fussy seizure roseolla sxbi_all acute resp gi newcough newrunnynose {
  gen count`x'=`x' if date>=(alleledate-14) & date<=(alleledate+7)
  gen allele`x'dt=date if count`x'==1 
  bysort patient: egen sum`x'=total(count`x'), missing
  gen allele`x'=sum`x' //set recur variable to equal tally of symptom days in recur window
  replace allele`x'=1 if sum`x'>0 & sum`x'!=. //set recur variable to a yes/no for presence of symptom
  replace allele`x'=. if alleledate==.
  drop count`x'
  drop sum`x'
}

*Generate symptom text. Try two methods to choose one per event in case marker variables have developed problems :)
tab1 allelepcpill allelefever allelecough allelenewcough allelerunnynose if oneper
tab1 allelenewrunnynose allelevomit allelediarrhea allelegeneralrash allelelocalrash if oneper
tab1 allelediaperrash allelefussy alleleseizure alleleroseolla if oneper

*Paper uses alleleacute for "any symptom"
tab1 alleleacute alleleresp allelegi if oneper

*Get viral load of allelswitches
sort sample
list sample alleleorder logq if alleleorder!=. //confirm that there is only one quant per pair listed
list sample alleleorder logq patient date in 1/6000 //confirm that this is true for no-allele samples too

gen allelestestcomplete=0
  replace allelestestcomplete=1 if alleleorder!=. & logq!=.

bysort allelestestcomplete: sum logq, detail
bysort allelestestcomplete: sum q_per_ml, detail

gen logq4=1 if logq>=4 & logq!=.
tab logq4 allelestestcomplete, miss


*IDSA 2010 shedding

gen newid="Subject 1" if patient==208
replace newid="Subject 2" if patient==221
replace newid="Subject 3" if patient==229
replace newid="Subject 4" if patient==233


twoway (scatter graphlogbov smo if bovpos==0, yaxis(1) sort mcolor(gray) msize(medium) msymbol(circle_hollow) mcolor(gray)) ///
    (scatter graphlogbov smo if bovpos==1, yaxis(1) sort mcolor(red) msize(medium) msymbol(circle)), ///
   by(newid, cols(2) iscale(0.6) imargin(small)) ///
  xscale(r(0 25)) ///
  yscale(r(1 9.5) axis(1)) ///
  ylabel(0(2)9.5, angle(horizontal) axis(1)) ///
  xtitle("Months", size(small)) ///
  ytitle("Log HBoV-1 Viral Load / mL Saliva", size(small)) ///
  xlabel(0(2)25) ///
  ysize(4) ///
  xsize(6) ///
  legend(order(1 "HBoV-1(-)" 2 "HBoV-1(+)") colgap(*.5) rows(1) size(small))
drop newid
  
  *IDSA recurrence
  gen newid="Subject 5" if patient==91
replace newid="Subject 6" if patient==202
replace newid="Subject 7" if patient==243
replace newid="Subject 8" if patient==284


twoway (scatter graphlogbov smo if bovpos==0, yaxis(1) sort mcolor(gray) msize(medium) msymbol(circle_hollow) mcolor(gray)) ///
    (scatter graphlogbov smo if bovpos==1, yaxis(1) sort mcolor(red) msize(medium) msymbol(circle)), ///
   by(newid, cols(2) iscale(0.6) imargin(small)) ///
  xscale(r(0 25)) ///
  yscale(r(1 9.5) axis(1)) ///
  ylabel(0(2)9.5, angle(horizontal) axis(1)) ///
  xtitle("Months", size(small)) ///
  ytitle("Log HBoV-1 Viral Load / mL Saliva", size(small)) ///
  xlabel(0(2)25) ///
  ysize(4) ///
  xsize(6) ///
  legend(order(1 "HBoV-1(-)" 2 "HBoV-1(+)") colgap(*.5) rows(1) size(small))
  drop newid
*/


*********************************************************************
*AIM 3: Incidence of HBoV and risk factors for initial acquisition
****************************************************************
*Graph month of first detections
format date %d

*Interpret this graph with caution!  Multiple age cohorts are beginning throughout time!
*histogram date if eventseq==1 & eventstart==1, bin(48) frequency

*Graph and describe first detections by age
*Use months in study for proxy right now
****************************************************************
*Descriptive statistics for age at primary (approximated by study dates) are included in the survival code
*gen hazardend = date if eventstart==1 & eventseq==1
*replace hazardend = end if tally==0 & oneper==1

*outsheet using "F:\SCH-Zerr Saliva\Source Data\CoxAnalysis_120810.raw" if hazardend!=., replace


*****************************************************************
*Recode survival event to require 2 positives to be an event
*****************************************************************

gen hazardend6 = date if event6start==1 & event6seq==1
  replace hazardend6=end if tally6==0 & oneper==1
  
outsheet using "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\CoxAnalysisEvent6_120114.raw" if hazardend6!=., replace

***************************************************

save "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\ArchiveAnalysisData_7d_30SEP.dta", replace
