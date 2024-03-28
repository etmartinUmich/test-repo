*********************************************************************************************
*This Coding File Analysis the Taqman Allele tests. Used in Natural History Manuscript
*********************************************************************************************

insheet using "U:/Secure/MartinEpi/Analysis/SCH-Zerr Saliva/Source Data/AlleleFile_nodups2.txt", clear
drop if position==.
sort sample
drop note*

list position sample if firstallele!=secondallele

*If no heterogenous alleles, set all alleles to first

gen allele=first
drop first second

drop if allele==""

list sample position allele if sample==19109|sample==21270
*4742 21270 had discrepant results - set to missing, needs to be reviewed later
*Also 4181 19109
drop if sample==21270 & position==4742 |sample==19109 & position==4181

*No hets listed here, will have to look into that later


reshape wide allele, i(sample) j(position)
drop if sample==.
sort sample



save "Z:/Analysis/SCH-Zerr Saliva/Source Data/Allele_11JUL2013.dta", replace

insheet using "Z:/Analysis/SCH-Zerr Saliva/Source Data/SampleList_formatted20APR.txt", clear
drop if sample==.
sort sample

merge 1:1 sample using "Z:/Analysis/SCH-Zerr Saliva/Source Data/Allele.dta"

gen ndate=date(date, "DMY", 2020)
drop date
gen date=ndate
drop ndate

tab sample if _merge==2
*We have 4 mystery samples that are probably typos:
*179
*6391
*9892
*21295

keep if _merge==3

gen complete=0
  replace complete=1 if allele2842!="" & allele3923!="" & allele4181!="" & allele4229!="" & allele4742!="" & allele4799!="" & allele4862!="" & allele4874!=""
*Abstract analysis

*3482 and 4502 are bad assays. 3482 has a few hits, but no variation found in individual patients.  Drop.
drop allele3482
drop allele4502

egen rownonmiss=rownonmiss(allele*), strok
egen rowmiss=rowmiss(allele*)
gen pcomprow=rownonmiss/(rowmiss+rownonmiss)
keep if pcomprow>0.5

sort patient sample

sort patient
egen tag=tag(patient)
bysort patient: gen total=_N
sort patient date
egen seq=seq(), from(1) by(patient)

*This code will find immediate mismatches but needs work to find mismatches after a period of missing values
*
*by patient: gen mismatch`x'=1 if seq!=1 & allele`x'!=allele`x'[_n-1] & allele`x'!="" & allele`x'[_n-1]!=""
*}

foreach x in 2842 3923 4181 4229 4742 4799 4862 4874 {
egen minmode`x'=mode(allele`x'), minmode by(patient)
gen onemiss`x'=1 if allele`x'!=minmode`x' & allele`x'!=""
by patient: egen mismatch`x'=max(onemiss`x')
}



sort patient
egen ncomp=total(complete), by(patient)

tab1 allele* //to get total number of specimens tests (add all the totals)

*It's reasonable to require at least 4 samples for the IDSA abstract.  Not many mismatch kids are lost (2 kids)

count if total>3 & tag
egen anymismatch=rowmax(mismatch*)
log using "Z:/Analysis/SCH-Zerr Saliva/Source Data/allelesbypt_11JUL.log", text replace
sort patient smo
by patient: list sample smo allele2842 allele3923 allele4181 allele4229 allele4742 allele4799 allele4862 allele4874 complete ncomp total mismatch4229 anymismatch
log close

log using "Z:/Analysis/SCH-Zerr Saliva/Source Data/allelesbypt_completeonly_11JUL.log", text replace
sort patient smo
by patient: list sample smo allele2842 allele3923 allele4181 allele4229 allele4742 allele4799 allele4862 allele4874 ncomp if pcomprow==1 & anymismatch==1
log close

count if tag
tab total if tag
count if tag

gen nmiss=total-ncomp
gen pcomp=1-(nmiss/total)

tab1 allele*, miss
tab1 allele*

tab ncomp if tag


egen ptmismatch=max(anymismatch), by(patient)
tab ptmismatch if tag, miss
tab ptmismatch if ncomp>3 & tag, miss



*Make sure to eliminate samples with no data at all
*How about setting a threshold for the abstract - 50% of sites reporting?

************************************



tab pcomprow, miss
*224 samples
 
count if tag
*71 patients had SNP data at at least 50% completeness

tab total if tag
*Range 1 to 10 samples
summ total if tag, detail
summ total if tag & total>1, detail
*48 had at least 2 samples


egen dtmax=max(date), by(patient)
egen dtmin=min(date), by(patient)
gen dtrange=dtmax-dtmin

summ dtrange if tag & dtrange>0, detail

tab ptmismatch if tag
*12
tab ptmismatch total if tag


*For Jan and Danielle - make a table of Completeness thresholds and numbers then completeness thresholds and SNPS coverage

*For cherry picking - First: Prioritize important sites
*get minor allele frequencies:
tab1 allele*



*For abstract - 

*No. pts. No samples per patient
*No sites, minor major allele % , missing data %
*No pts, with 4 or more samples per pt
*No pts with heterogeneity found based on more than one SNP pattern.
*No types per patient

*TO DO:
*Clean up missing data
*Get interim samples for heterogeneity kids
*Figure out if we have mixed types
*sequence confirmation?


*Get date of allele change to do symptoms analysis.

sort patient date
by patient: gen order = _n
gen sorttotal= _N
gen switch=0

foreach x in 2842 3923 4181 4229 4742 4799 4862 4874 {
sort patient date
by patient: gen lagallele`x'=allele`x'[_n-1] if order!=1
gen switchdt`x'=date if allele`x'!=lagallele`x' & allele`x'!="" & lagallele`x'!=""
gen prevend`x'=date[_n-1] if allele`x'!=lagallele`x' & allele`x'!="" & lagallele`x'!=""
gen datesince`x'=prevend`x'-switchdt`x' 
replace switch=1 if switchdt`x'!=.
list patient sample date order allele`x' lagallele`x' switchdt`x' datesince`x' if mismatch`x'==1
}
*This gives median time between previous resolution and next genotype
list patient sample datesince*
egen datesince=rowmax(datesince2842 datesince3923 datesince4181 datesince4229 datesince4742 datesince4799 datesince4862 datesince4874)
list patient datesince if datesince!=.
summ datesince, detail

sort patient switch date
by patient switch: gen switch_n=_n if switch==1
gen alleleorder=0
replace alleleorder=switch_n if alleleorder==0 & switch_n!=.  

sort patient date   
by patient: replace alleleorder=alleleorder[_n-1] if order!=1 & alleleorder==0 | alleleorder==.

*confirm that the patient with two switches has two separate alleles
sort patient date
list patient date alleleorder allele* if patient==221

sort patient date
gen switchstart=switch

sort patient date
by patient: egen maxswitch=max(alleleorder)
tab maxswitch if tag

sort patient date
by patient: list date alleleorder allele* if maxswitch==1 | maxswitch==2

keep sample patient tag date alleleorder switchstart allele* lagallele* switchdt*
save "Z:/Analysis/SCH-Zerr Saliva/Source Data/alleleswitches.dta", replace
list sample patient date allele2842 switchdt2842
*For the sample with two switches, only the earliest switch was used if you do neworder==1.
describe
