**********************************************************************************************
*This code created the figure that graphs symptom prevalence in the Natural history manuscript
*************************************************************************************************

clear
set memory 500m

use "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\BoVPrimaryEvents_100714.dta", clear
sort patient date
describe
*keep if date<=(primary6date+21) & date>=(primary6date-28)
describe
egen ptag=tag(patient)

gen stddate=date-primary6date
tab stddate


/*
egen pcpstart=min(gprimarypcpilldt), by(patient)
gen pcpdiff=pcpstart-primary6date
tab pcpdiff if ptag==1, miss
summ pcpdiff if ptag==1, detail

*Make a variable of date of symptom, standardized by primary onset date
gen stdpcp=gprimarypcpilldt-primary6date
tab stdpcp pcpill
histogram stdpcp, bin(40) freq

*<symptom>after variables are based on the symptom start date that was filled into all observations by patient
*For this reason, we can select any observation for a patient and get a pcpafter result
gen pcpafter=0 if pcpdiff!=.
  replace pcpafter=1 if pcpdiff>=0 & pcpdiff!=.
tab pcpafter if ptag==1, miss


egen feverstart=min(gprimaryfeverdt), by(patient) 
gen feverdiff=feverstart-primary6date
tab feverdiff if ptag==1, miss
summ feverdiff if ptag==1, detail

gen feverafter=0 if feverdiff!=.
replace feverafter=1 if feverdiff>0 & feverdiff!=.
tab feverafter if ptag==1, miss
*/
*Make a variable of date of symptom, standardized by primary onset date

bysort stddate: egen pcpilldayct=total(pcpill)
bysort stddate: egen feverdayct=total(fever)
bysort stddate: egen newcoughdayct=total(newcough)

egen oneperdate=tag(stddate)


set scheme s1mono
twoway (area pcpilldayct stddate, yaxis(1))(mspline pcpilldayct stddate) if oneperdate==1 & stddate<=21 & stddate>=-28, legend(off) ytitle("PCP visit (n)") xtitle("") xlabel(none) yscale(range(0 5)) ysize(5) xsize(16) xline(-7, lwidth(1)) xline(7, lwidth(1))
graph save pcpill, replace
twoway (area feverdayct stddate, yaxis(1))(mspline feverdayct stddate) if oneperdate==1 & stddate<=21 & stddate>=-28, legend(off) ytitle("Fever (n)") xtitle("") xlabel(none) yscale (range(0 5)) ysize(5) xsize(16) xline(-7, lwidth(1)) xline(7, lwidth(1))
graph save fever, replace
twoway (area newcoughdayct stddate, yaxis(1))(mspline newcoughdayct stddate) if oneperdate==1 & stddate<=21 & stddate>=-28, legend(off) ytitle("New Cough (n)") xtitle("Days to Start of HBoV-1 Primary Event") yscale(range(0 5)) ysize(5) xsize(16) xline(-7, lwidth(1)) xline(7, lwidth(1))
graph save newcough, replace

graph combine pcpill.gph fever.gph newcough.gph, cols(1) xcommon ycommon imargin(0 0 0 0) scheme(s1mono) ///THIS IS FIGURE 1!!!

*NEW FIGURE 1 is here:
set scheme s1mono
twoway (bar pcpilldayct stddate, yaxis(1))(mspline pcpilldayct stddate) if oneperdate==1 & stddate<=21 & stddate>=-28, legend(off) ytitle("PCP visit (n)") xtitle("") xlabel(none) yscale(range(0 5)) ysize(5) xsize(16) xline(-7, lwidth(1)) xline(7, lwidth(1))
graph save pcpillbar, replace
twoway (bar feverdayct stddate, yaxis(1))(mspline feverdayct stddate) if oneperdate==1 & stddate<=21 & stddate>=-28, legend(off) ytitle("Fever (n)") xtitle("") xlabel(none) yscale (range(0 5)) ysize(5) xsize(16) xline(-7, lwidth(1)) xline(7, lwidth(1))
graph save feverbar, replace
twoway (bar newcoughdayct stddate, yaxis(1))(mspline newcoughdayct stddate) if oneperdate==1 & stddate<=21 & stddate>=-28, legend(off) ytitle("New Cough (n)") xtitle("Days to Start of HBoV-1 Primary Event") yscale(range(0 5)) ysize(5) xsize(16) xline(-7, lwidth(1)) xline(7, lwidth(1))
graph save newcoughbar, replace

graph combine pcpillbar.gph feverbar.gph newcoughbar.gph, cols(1) xcommon ycommon imargin(0 0 0 0) scheme(s1mono) ///THIS IS FIGURE 1!!!






twoway (area pcpilldayct stddate, yaxis(1))(area feverdayct stddate, yaxis(1))(area newcoughdayct stddate, yaxis(1)) if oneperdate==1, ysize(4) xsize(16)

twoway (mspline pcpilldayct stddate)(mspline feverdayct stddate)(mspline newcoughdayct stddate) if oneperdate==1, ysize(4) xsize(8)


gen stdpcpill=gprimarypcpilldt-primary6date
tab stdpcpill pcpill
sort stdpcpill
by stdpcpill: egen pcpilldayct=total(pcpill)
egen onepcpill=tag(stdpcpill)
twoway (mspline pcpilldayct stdpcpill)(line pcpilldayct stdpcpill) if onepcpill==1, 

gen stdpcpill=gprimarypcpilldt-primary6date
tab stdpcpill pcpill
sort stdpcpill
by stdpcpill: egen pcpilldayct=total(pcpill)
egen onepcpill=tag(stdpcpill)
twoway (mspline pcpilldayct stdpcpill if onepcpill!=1)(line pcpilldayct stdpcpill if onepcpill!=1)

gen stdnewcough=gprimarynewcoughdt-primary6date
tab stdnewcough newcough
sort stdnewcough
by stdnewcough: egen newcoughdayct=total(newcough)
egen onenewcough=tag(stdnewcough)
twoway (mspline newcoughdayct stdnewcough if onenewcough!=1)(line newcoughdayct stdnewcough if onenewcough!=1)  


gen stdfever=gprimaryfeverdt-primary6date
tab stdfever fever
*histogram stdfever, bin(40) freq 

sort stdfever 
by stdfever: egen feverdayct=total(fever)
egen onefever=tag(stdfever)
replace feverdayct=. if onefever!=1
replace stdfever=. if onefever!=1
twoway (mspline feverdayct stdfever)(line feverdayct stdfever)



egen newcoughstart=min(gprimarynewcoughdt), by(patient)
gen newcoughdiff=newcoughstart-primary6date
tab newcoughdiff if ptag==1, miss
summ newcoughdiff if ptag==1, detail

gen newcoughafter=0 if newcoughdiff!=.
replace newcoughafter=1 if newcoughdiff>0 & newcoughdiff!=.
tab newcoughafter if ptag==1, miss

egen newrunnynosestart=min(gprimarynewrunnynosedt), by(patient) 
gen newrunnynosediff=newrunnynosestart-primary6date
tab newrunnynosediff if ptag==1, miss
summ newrunnynosediff if ptag==1, detail

gen newrunnynoseafter=0 if newrunnynosediff!=.
replace newrunnynoseafter=1 if newrunnynosediff>0 & newrunnynosediff!=.
tab newrunnynoseafter if ptag==1, miss

egen diarrheastart=min(gprimarydiarrheadt), by(patient) 
gen diarrheadiff=diarrheastart-primary6date
tab diarrheadiff if ptag==1, miss
summ diarrheadiff if ptag==1, detail

gen diarrheaafter=0 if diarrheadiff!=.
replace diarrheaafter=1 if diarrheadiff>0 & diarrheadiff!=.
tab diarrheaafter if ptag==1, miss

egen fussystart=min(gprimaryfussydt), by(patient) 
gen fussydiff=fussystart-primary6date
tab fussydiff if ptag==1, miss
summ fussydiff if ptag==1, detail

gen fussyafter=0 if fussydiff!=.
replace fussyafter=1 if fussydiff>0 & fussydiff!=.
tab fussyafter if ptag==1, miss

egen acutestart=min(gprimaryacutedt), by(patient) 
gen acutediff=acutestart-primary6date
tab acutediff if ptag==1, miss
summ acutediff if ptag==1, detail

gen acuteafter=0 if acutediff!=.
replace acuteafter=1 if acutediff>0 & acutediff!=.
tab acuteafter if ptag==1, miss

egen respstart=min(gprimaryrespdt), by(patient) 
gen respdiff=respstart-primary6date
tab respdiff if ptag==1, miss
summ respdiff if ptag==1, detail

gen respafter=0 if respdiff!=.
replace respafter=1 if respdiff>0 & respdiff!=.
tab respafter if ptag==1, miss

egen gistart=min(gprimarygidt), by(patient) 
gen gidiff=gistart-primary6date
tab gidiff if ptag==1, miss
summ gidiff if ptag==1, detail

gen giafter=0 if gidiff!=.
replace giafter=1 if gidiff>0 & gidiff!=.
tab giafter if ptag==1, miss

*twoway (kdensity feverdiff) (kdensity newcoughdiff) (kdensity pcpdiff) (kdensity newrunnynosediff) (kdensity diarrheadiff) (kdensity fussydiff)

*twoway (kdensity acutediff) (kdensity respdiff) (kdensity gidiff)

twoway (kdensity feverdiff) (kdensity newcoughdiff) (kdensity pcpdiff)

 
*twoway (kdensity respdiff)
