********************************************************************************************************
*This code looks at time to HBoV acquisition. Used in Natural History Manuscript
********************************************************************************************************


/* This code requires at least two positive samples to establish a primary event*/
/*Everyone with no event6 is right censored at study end*/
insheet using "U:\Secure\MartinEpi\Analysis\SCH-Zerr Saliva\Source Data\CoxAnalysisEvent6_120114.raw", clear

gen time=hazardend6-start
summ time, detail //this is the same as age at acquisition
replace event6=0 if tally6==0
tab event6, miss

stset time, id(patient) failure(event6)
*This is only time to single event.
*So, incidence rate is incidence of primary infection
*Confirm that only 66 events are included in all analyses (first events)
strate, per(36500)
sts graph
sts list

*Try changing to child-months
strate, per(3000)
sts list
sts graph, 

*Median age at acquisition is now 13m, because we are requiring two positives!!

foreach x in daycarestart maternaldob infantdob {
gen `x'_f = date(`x', "DMY")
drop `x'
gen `x' = `x'_f
format `x' %d
drop `x'_f
describe `x'
}

gen male=1 if gender=="M"
replace male=0 if gender=="F"

gen bfed=0 if breastfed=="N"
replace bfed=1 if breastfed=="Y"

gen dc=0 if daycare=="N"
replace dc=1 if daycare=="Y"

gen pg=0 if playgroup=="N"
replace pg=1 if playgroup=="Y"

gen anygroup=1 if pg==1|dc==1
replace anygroup=0 if pg==0 & dc==0
tab anygroup
sts test anygroup, logrank
strate anygroup, per(36500)
tab anygroup event

gen birthmonth=month(start)
gen birthseas=.
replace birthseas=0 if birthmonth==6|birthmonth==7|birthmonth==8
replace birthseas=1 if birthmonth==9|birthmonth==10|birthmonth==11
replace birthseas=2 if birthmonth==12|birthmonth==1|birthmonth==2
replace birthseas=3 if birthmonth==3|birthmonth==4|birthmonth==5
tab birthseas

*sts graph, by (birthseas) failure
sts test birthseas, logrank
xi: stcox i.birthseas
sts graph, by(birthseas) failure

gen fallbaby=0
replace fallbaby=1 if birthseas==1
sts graph, by (fallbaby) failure
sts test fallbaby, logrank

sts graph, by(birthseas) failure
sts test birthseas, logrank

gen race4=0
replace race4=1 if race==4
sts graph, by (race4) failure
sts test race4, logrank


stcox maternalage
sts graph, by (matagecat) failure
sts test matagecat, logrank

stcox ageatbfstop

stcox fallbaby
stcox male
xi: stcox i.race
stcox bfed
stcox dc
xi: stcox i.dcfulltime
stcox pg
stcox anysibs
xi: stcox i.numberoffamily18yrs
stcox i.matagecat
stcox race4
*Look at twins:
list patient event time if twin==1

stcox yearn

save "F:\SCH-Zerr Saliva\Source Data\CoxAnalysisUnsplit14JUN11.dta", replace
use "F:\SCH-Zerr Saliva\Source Data\CoxAnalysisUnsplit14JUN11.dta", clear

*Try daycare as a time-varying variable
gen difftest=time-ageatdcstart
tab difftest, miss
count if difftest>0 & difftest!=.
gen dcstartsplit=1 if difftest>0 & difftest!=.
*17 kids started daycare AFTER primary infection

stsplit dcsplit, after(time=ageatdcstart) at(1)
gen newdaycare=0 if dc==0
  replace newdaycare=dcsplit if (ageatdcstart<time) & ageatdcstart!=.
  replace newdaycare=0 if ageatdcstart>time
*list patient time event dc difftest dcsplit newdaycare
  
gen enddiff=time-ageatdcend
tab enddiff
*no one ended before acquisition.  No need to do this split
tab event dc, miss
stcox dc
stcox newdaycare

****Patient 179 was bov positive at first sample.  Is not included in cox models

  *Try breastfeeding as a time-varying variable
    *BE CAREFUL - weaning and daycare start may be collinear
stsplit bfsplit, after(time=ageatbfstop) at(1)
gen newbf=bfed
  replace newbf=abs(bfsplit-1) if untilwhenc==1
*list patient time event bfsplit bfed untilwhenc newbf
stcox bfed
stcox newbf newdaycare

stepwise, pr(.20) pe(0.10) forward: stcox fallbaby anysibs male bfed newdaycare newbf pg race4
stepwise, pr(.20) pe(0.10) forward: stcox fallbaby numberoffamily male bfed newdaycare newbf pg race4


/*Code below looks at playgroup and playgroup plus daycare.  Neither were significant
use "F:\SCH-Zerr Saliva\Source Data\CoxAnalysisUnsplit120810.dta", clear

gen anystart=min(ageatpgstart, ageatdcstart)
gen anystop=max(ageatpgend, ageatdcend)
gen gpstartsplit=1 if (time-anystart)>0 & (time-anystart)!=.
tab gpstartsplit
gen enddiff2=time-anystop
tab enddiff2
*no group activity stop before acquisition.  No need to do this split

stsplit gpsplit, after(time=anystart) at(1)
gen newgp=anygroup
  replace newgp=gpsplit if gpstartsplit==1
  replace newgp=0 if anystart>time
list patient time event gpsplit anystart ageatpgstart ageatdcstart anygroup newgp
stcox anygroup
stcox anygroup fallbaby
*Association with "any group" i.e. daycare or playgroup is confounded by fallbaby

use "F:\SCH-Zerr Saliva\Source Data\CoxAnalysisUnsplit120810.dta", clear
gen pgdiff=time-ageatpgstart
gen pgstartsplit=1 if pgdiff!=. & pgdiff>0

stsplit pgsplit, after(time=ageatpgstart) at(1)
gen newpg=pg
  replace newpg=pgsplit if pgstartsplit==1
  replace newpg=0 if ageatpgstart>time
list patient time event pgsplit ageatpgstart pg newpg
stcox newpg
stcox newpg fallbaby
*/
