********************************************************************
*This file combines the data for the first set of pilot testing and 
*does some basic descriptive analyses and graphs.

*Analysis was originally done using the pilot_50.txt dataset generated
*from Virology lab file "BoV PCR on Ext SS 052809.xls".

*Analysis now uses pilot_50_v2.txt file, generated from Virology lab
*file "HHV6 samples for BoV and AdV 062409.xls"
*Emily Martin 2009
********************************************************************


clear

set memory 700m

/*
*****************************
*Below is a record of how pilot data was cleaned.
********************************
insheet using "O:\Analyses\Zerr Saliva\SourceData\pilot_50_v2.txt", clear
describe

****Format Date field and code
foreach x in date{
gen `x'_f = date(`x', "MDY")
drop `x'
gen `x' = `x'_f
format `x' %d
drop `x'_f
describe `x'
}


****Format virology results
*Drop the adv results.  We aren't working on that in this analysis.
drop adv

*"None detected" results were marked at -9 by me in the excel file.
*Change those to 0 for the purposes of this analysis.
*A few of the BoV results were randomly redone.  Ignore the redo's (bov_run2)
*unless there was none detected during run1 and virus detected in run2.

replace bov_run1=0 if bov_run1==-9
replace bov_run2=0 if bov_run2==-9

gen q_per_pcr=bov_run1
replace q_per_pcr=bov_run2 if (bov_run1==.|bov_run1==0)&(bov_run2!=. & bov_run2!=0)

*Virology lab reports quantitation in copies per pcr.  To get copies per ml
*multiple everything by 167.

gen q_per_ml=q_per_pcr*167

gen bovpos=0
  replace bovpos=1 if q_per_ml>0 & q_per_ml!=.
tab bovpos, miss
*Pilot Set: 261 BOV+; 916 BOV-

gen pilotset=1

describe
*1177 total observations
sort patient date
save "O:\Analyses\Zerr Saliva\SourceData\pilot_50_formatted.dta", replace
***************************************************************************
*/

/*
Sort and prep pilot data file.
*/

use "C:\SCH-Zerr Saliva\Source Data\pilot_50_formatted.dta", clear
sort patient date
save "C:\SCH-Zerr Saliva\Source Data\pilot_50_formatted.dta", replace

*************************************************************************
*Add in full grant results from Jane
*Previously extracted samples:
****Uses SourceData\BoVResults_Oldextracted_5-12-10.txt from BoV PCR Ext SS 041310.xls from Jane

*Newly extracted samples:
****Uses SourceData\BoVResults_Newextracted_5-13-10.txt from BoV PCR unext SS 4-16 to 5-13-052110.xlsx from Jane
****Uses SourceData\BoVResults_Newextracted_5-14to6-14.txt from BoV PCR unext SS 5-14 to 6-14 061410.xls from Jane
****Uses SourceData\BoVResults_Newextracted_9-20.txt from "Results SS for BoV July2010Request 092010.xls" from Jane

*Update run performed by JP:
****Uses SourceData\Updatebov_ptdate.txt from JP's notes. All tests were negative
*************************************************************************
insheet using "C:\SCH-Zerr Saliva\Source Data\BoVResults_Newextracted_5-13-10.txt", clear
gen ext="batch1"
describe
tab bovct, miss
sort patient date
save data1, replace

insheet using "C:\SCH-Zerr Saliva\Source Data\BoVResults_Newextracted_5-14to6-14.txt", clear
gen ext="batch2"
describe
tab bovct, miss
drop if bovct==""
tab bovct, miss
sort patient date
save data2, replace

insheet using "C:\SCH-Zerr Saliva\Source Data\BoVResults_Newextracted_9-20.txt", clear
gen ext="batch3"
describe
tab bovct, miss
drop if bovct=="miss_samp"
tab bovct, miss
sort patient date
save data3, replace

*format bov fill in run from Oct 2012

insheet using "C:\SCH-Zerr Saliva\Source Data\updatebov_ptdate.txt", clear
gen ext="oct2012"
describe
sort patient date
save data4, replace

insheet using "C:\SCH-Zerr Saliva\Source Data\BoVResults_Oldextracted_5-12-10.txt", clear
gen ext="old"
describe
tab bovct, miss
sort patient date

merge m:m patient date using data1
tab _merge
drop _merge
sort patient date

merge m:m patient date using data2
tab _merge
drop _merge

merge m:m patient date using data3
tab _merge
drop _merge

*add updated bov runs from Oct 2012
merge m:m patient date using data4
tab _merge

describe
tab bovct, miss

gen bovpos=0 if bovct=="neg"
  replace bovpos=1 if q_per_pcr>0 & q_per_pcr!=.
  replace bovpos=0 if updatebovpos==0
tab bovpos, miss

tab bovpos updatebovpos if _merge==3
*******Numbers of observations, BoV+ and BoV- are all fine after the merges:
*******3248 total tested samples.  2646 negative. 602 positive.
drop _merge

foreach x in date{
gen `x'_f = date(`x', "MDY")
drop `x'
gen `x' = `x'_f
format `x' %d
drop `x'_f
describe `x'
}

gen q_per_ml=q_per_pcr*167

drop labid

gen pilotset=0

describe
*3114 total observations
sort patient date

merge m:m patient date using "C:\SCH-Zerr Saliva\Source Data\pilot_50_formatted.dta"
tab _merge
*Both datasets are unique
tab pilotset _merge
drop _merge
describe
*Note bov_run1 & bov_run2 are only in the pilot set; bovct is only in the grant set
*Adds 1177 observations from the pilot data.  Now is 4291 total observations

replace ext="pilot" if ext==""

*Mark the earliest bov positive
sort patient date
sort patient bovpos, stable
by patient bovpos: gen order=_n
gen bovstartdt=date if order==1 & bovpos==1
sort patient date
list patient date bovpos bovstartdt order if patient==208
drop order

bysort patient: egen fillbovstartdt=min(bovstartdt)
replace bovstartdt=fillbovstartdt
drop fillbovstartdt
gen bovstart=1 if bovstartdt==date
sort patient date
list patient date bovpos bovstartdt bovstart if patient==208

*Identify continuous shedding patterns - mark events, defined as continuous detections with no more than 1 interim negative,
*(i.e. - 2 negatives=the end of the event)
sort patient date
by patient: gen bovlag1 = bovpos[_n-1]
by patient: gen bovlag2 = bovpos[_n-2]
by patient: gen datelag1 = date[_n-1]
by patient: gen datelag2 = date[_n-2]
by patient: gen bovlead1 = bovpos[_n+1]

gen event=1 if bovpos==1 
replace event=1 if bovpos==0 & (bovlead1==1) & (bovlag1==1 | datelag1==bovstartdt | date==bovstartdt)

sort patient date
list patient date bovstartdt bovstart bovpos bovlag1 bovlag2 datelag1 datelag2 event if patient==208|334|221|285|233|271

*Mark the start of each event
sort patient date
gen eventstart=1 if (event==1 & event[_n-1]==.) | date==bovstartdt
gen eventstartdt=date if eventstart==1
list patient date bovstartdt bovstart bovpos event eventstart eventstartdt if patient==334

*Sequentially tag the events (Event 1, Event 2, etc.)
sort patient eventstart date
by patient eventstart: egen eventseq=seq(), from(1)
replace eventseq=. if eventstart!=1 
list patient date eventstart eventstartdt eventseq if patient==334

*Fill in the event sequence into the whole event
sort patient date
replace eventseq = eventseq[_n-1] if eventseq >= . & event==1
list patient date eventstart eventstartdt eventseq if patient==334

*Determine length of events and fill in event start date
sort patient eventseq date
by patient eventseq: egen neweventstartdt=min(eventstartdt)
replace eventstartdt=neweventstartdt if event!=1
drop neweventstartdt
by patient eventseq: egen eventenddt=max(date)
replace eventenddt=. if event!=1
gen eventdur=eventenddt-eventstartdt if event==1
sort patient date
replace eventdur = eventdur[_n-1] if eventdur >= . & event==1

sort patient date
list patient date event eventstart eventstartdt eventenddt  eventdur eventseq  if patient==334

*********************************************
*Variable coding - BoV test data
*********************************************
sort patient date
by patient: egen start=min(date)
by patient: egen end=max(date)
gen totaltime=end-start

gen smo=(date-start)/30.25

*Tally total number of positive detections in each patient
sort patient date
by patient: egen tally=total(bovpos)
tab tally if date==start

*Mark patients tested as part of the pilot phase 
sort patient date
by patient: egen newpilot=max(pilotset)
replace pilotset=newpilot
drop newpilot
tab pilotset, miss

*Mark patients tested for BoV at any time
sort patient date
by patient: egen anytest=max(bovpos)
replace anytest=1 if anytest>=0
tab anytest if date==start, miss
tab tally anytest if date==start, miss

*Create a log bov variable with the 0 point moved slightly higher
*to help with graph visualization
gen graphlogbov=log10(q_per_ml)
  replace graphlogbov=0.5 if bovpos==0

sort patient
save "C:\SCH-Zerr Saliva\Source Data\SampleData_061113.dta", replace

************************
*Format demographics dataset
****************************
insheet using "C:\SCH-Zerr Saliva\Source Data\demo.txt", clear
gen patient=id
foreach x in maternaldob infantdob stopdate untilwhen daycarestart daycareend playstart playend lastfuvisit {
gen `x'_f = date(`x', "MDY")
drop `x'
gen `x' = `x'_f
format `x' %d
drop `x'_f
describe `x'
}

*Merge with BoV testing results
sort patient
merge 1:m patient using "C:\SCH-Zerr Saliva\Source Data\SampleData_061113.dta"
tab _merge
drop if _merge==1
save "C:\SCH-Zerr Saliva\Source Data\SampleDemsData_061113.dta", replace
*4291 BoV tests, with demographic data
*************************************
*Add in symptom data
*************************************
/*
Survey data already cleaned as follows:
insheet using "C:\Users\Emily Martin\Documents\Analyses\SCH-Zerr Saliva\Source Data\survey1.txt", clear

foreach x in date{
gen `x'_f = date(`x', "MDY")
drop `x'
gen `x' = `x'_f
format `x' %d
drop `x'_f
describe `x'
}
*drop pcr variable - it refers to HHV6
drop pcr hhv6 hhv6copies

save "C:\Users\Emily Martin\Documents\Analyses\SCH-Zerr Saliva\Source Data\survey1.dta", replace 

insheet using "C:\Users\Emily Martin\Documents\Analyses\SCH-Zerr Saliva\Source Data\survey2.txt", clear

foreach x in date{
gen `x'_f = date(`x', "MDY")
drop `x'
gen `x' = `x'_f
format `x' %d
drop `x'_f
describe `x'
}
*drop pcr variable - it refers to HHV6
drop pcr hhv6 hhv6copies

save "C:\Users\Emily Martin\Documents\Analyses\SCH-Zerr Saliva\Source Data\survey2.dta", replace 
*/
***********************
*Append survey data to the end of BoV result data
*Remember that the demographic data is attached to the BoV result data
************************
use "C:\SCH-Zerr Saliva\Source Data\SampleDemsData_061113.dta", clear
append using "C:\SCH-Zerr Saliva\Source Data\survey1.dta"
append using "C:\SCH-Zerr Saliva\Source Data\survey2.dta"


*Drop kids followed for less than 18 months for consistency sake
drop if totaltime<547

********************************************
*Variable Coding
**********************************************
*HCP data coding for dr visit type:
*-1=reason unknown; 0=checkup; 1=acute illness; 2=trauma
gen pcpill=0 if drvisittype!=.
replace pcpill=1 if drvisittype==1
tab pcpill, miss


save "F:\SCH-Zerr Saliva\Source Data\GrantRes1_SampleSx_061113.dta", replace

clear 
clear matrix
set memory 700m

use "F:\SCH-Zerr Saliva\Source Data\GrantRes1_SampleSx_061113.dta", clear
drop _merge
sort patient date
merge n:n patient date using "U:\My Documents\Bocavirus\SourceData\alleleswitches.dta"
save "F:\SCH-Zerr Saliva\Source Data\GrantRes1_SampleSx_071313.dta", replace

