* NOTE: The original analysis included a small number of participant-level
* corrections taken from field notes (individual GDS item values, missing-reason
* codes, and participation-level reclassifications). These are keyed on study IDs
* and are omitted from this public version to protect participant privacy.
* The analysis dataset itself is not public, so this script documents the
* analysis steps but cannot be run without access to the SWISS100 data.

*==============================================================================*
* SWISS100 - GDS item non-response & missing data                     V7
* Requires: dtable (Stata 18+), estout (ssc install estout)
*==============================================================================*

version 19
clear all
set more off

* ----- EDIT THESE TWO PATHS ONCE --------------------------------------------*
cd "PATH/TO/PROJECT/FOLDER"
global data "vulnerability_data_V5_small.dta"
global out  "PATH/TO/PROJECT/FOLDER/output"

use "$data", clear

*##############################################################################*
* PART A - MANUAL RECODES + DERIVATIONS (verbatim from V1)
*##############################################################################*


* ---- drop pre-existing derived score vars (ROBUST: capture) ----------------*
* V5 already contains some of these; capture drop ignores absent ones.
foreach v in wb_gds_nonmiss wb_gds_score wb_gds_nonmiss2 wb_gds_score2 ///
    wb_gds_8i wb_gds_3i wb_gds_rowmean wb_gds_score_meansub wb_gds_rowmean_d ///
    wb_gds_score_meansub2 wb_gds_score_meansubR wb_gds_score_X15 wb_gds_score_X15R ///
    wb_gdsSS_rowmean wb_gdsSS_rowmean_d wb_gdsSS_score_meansub wb_gdsSS_score_meansub2 ///
    wb_gdsSS_score_meansubR wb_gds_score_X4 wb_gds_score_X4R {
    capture drop `v'
}
capture drop *_copy
capture drop *_meansub
capture drop *_meansub2
capture drop *_meansubSS
capture drop *_meansub3

label define id_particip_ 1 "L0 - Proxy replacement" 2 "L1 - Low" 3 "L2 - Intermediate" 4 "L3 - Full", modify

* ---- recompute scores from raw items ---------------------------------------*
* complete-case GDS-15
egen wb_gds_rowmiss=rownonmiss(wb_gds1r wb_gds2 wb_gds3 wb_gds4 wb_gds5r wb_gds6 wb_gds7r wb_gds8 wb_gds9 wb_gds10 wb_gds11r wb_gds12 wb_gds13r wb_gds14 wb_gds15)
egen wb_gds_score=rowtotal(wb_gds1r wb_gds2 wb_gds3 wb_gds4 wb_gds5r wb_gds6 wb_gds7r wb_gds8 wb_gds9 wb_gds10 wb_gds11r wb_gds12 wb_gds13r wb_gds14 wb_gds15) if wb_gds_rowmiss>14
rename wb_gds_rowmiss wb_gds_nonmiss
label var wb_gds_nonmiss "Number of GDS items non missing"

egen wb_gds_nonmiss2=rownonmiss(wb_gds3 wb_gds4 wb_gds8 wb_gds14)
egen wb_gds_score2=rowtotal(wb_gds3 wb_gds4 wb_gds8 wb_gds14) if wb_gds_nonmiss2>3

* 8-item rule (GDS-15) and 3-item rule (GDS-4)
egen wb_gds_8i=rowtotal(wb_gds1r wb_gds2 wb_gds3 wb_gds4 wb_gds5r wb_gds6 wb_gds7r wb_gds8 wb_gds9 wb_gds10 wb_gds11r wb_gds12 wb_gds13r wb_gds14 wb_gds15) if wb_gds_nonmiss>7
label var wb_gds_8i "Short GDS score - 8 items rule"
egen wb_gds_3i=rowtotal(wb_gds3 wb_gds4 wb_gds8 wb_gds14) if wb_gds_nonmiss2>2
label var wb_gds_3i "Short-short GDS score - 3 items rule"

* mean-substitution (truncated) GDS-15
foreach var of varlist wb_gds1r wb_gds2 wb_gds3 wb_gds4 wb_gds5r wb_gds6 wb_gds7r wb_gds8 wb_gds9 wb_gds10 wb_gds11r wb_gds12 wb_gds13r wb_gds14 wb_gds15 {
    gen `var'_copy=`var'
    replace `var'_copy=. if missing(`var')
}
egen wb_gds_rowmean=rowmean(wb_gds1r_copy-wb_gds15_copy)
gen wb_gds_rowmean_d=0 if wb_gds_rowmean<0.5 & wb_gds_rowmean!=.
replace wb_gds_rowmean_d=1 if wb_gds_rowmean>=0.5 & wb_gds_rowmean!=.
foreach var of varlist wb_gds1r wb_gds2 wb_gds3 wb_gds4 wb_gds5r wb_gds6 wb_gds7r wb_gds8 wb_gds9 wb_gds10 wb_gds11r wb_gds12 wb_gds13r wb_gds14 wb_gds15 {
    gen `var'_meansub=`var'
    replace `var'_meansub=wb_gds_rowmean_d if missing(`var'_copy)
}
egen wb_gds_score_meansub=rowtotal(wb_gds1r_meansub-wb_gds15_meansub)

* mean-substitution NO truncation + rounded
foreach var of varlist wb_gds1r wb_gds2 wb_gds3 wb_gds4 wb_gds5r wb_gds6 wb_gds7r wb_gds8 wb_gds9 wb_gds10 wb_gds11r wb_gds12 wb_gds13r wb_gds14 wb_gds15 {
    gen `var'_meansub2=`var'
    replace `var'_meansub2=wb_gds_rowmean if missing(`var'_copy)
}
egen wb_gds_score_meansub2=rowtotal(wb_gds1r_meansub2-wb_gds15_meansub2)
gen wb_gds_score_meansubR=round(wb_gds_score_meansub2,1)
replace wb_gds_score_meansub=.  if wb_gds_nonmiss==0
replace wb_gds_score_meansub2=. if wb_gds_nonmiss==0
replace wb_gds_score_meansubR=. if wb_gds_nonmiss==0

* ipsative row-mean x 15 (Koyama method)
gen wb_gds_score_X15=wb_gds_rowmean*15
gen wb_gds_score_X15R=round(wb_gds_score_X15,1)

* GDS-4 mean-sub (truncated) + no-trunc + ipsative x4
egen wb_gdsSS_rowmean=rowmean(wb_gds3_copy wb_gds4_copy wb_gds8_copy wb_gds14_copy)
gen wb_gdsSS_rowmean_d=0 if wb_gdsSS_rowmean<0.5 & wb_gdsSS_rowmean!=.
replace wb_gdsSS_rowmean_d=1 if wb_gdsSS_rowmean>=0.5 & wb_gdsSS_rowmean!=.
foreach var of varlist wb_gds3 wb_gds4 wb_gds8 wb_gds14 {
    gen `var'_meansubSS=`var'
    replace `var'_meansubSS=wb_gdsSS_rowmean_d if missing(`var'_copy)
}
egen wb_gdsSS_score_meansub=rowtotal(wb_gds3_meansubSS wb_gds4_meansubSS wb_gds8_meansubSS wb_gds14_meansubSS)
foreach var of varlist wb_gds3 wb_gds4 wb_gds8 wb_gds14 {
    gen `var'_meansub3=`var'
    replace `var'_meansub3=wb_gdsSS_rowmean if missing(`var'_copy)
}
egen wb_gdsSS_score_meansub2=rowtotal(wb_gds3_meansub3 wb_gds4_meansub3 wb_gds8_meansub3 wb_gds14_meansub3)
gen wb_gdsSS_score_meansubR=round(wb_gdsSS_score_meansub2,1)
gen wb_gds_score_X4=wb_gdsSS_rowmean*4
gen wb_gds_score_X4R=round(wb_gds_score_X4,1)
replace wb_gdsSS_score_meansub=.  if wb_gds_nonmiss==0
replace wb_gdsSS_score_meansub2=. if wb_gds_nonmiss==0
replace wb_gdsSS_score_meansubR=. if wb_gds_nonmiss==0

label var wb_gds_score          "GDS score Short Form (GDS-15) - full case"
label var wb_gdsSS_score_meansub "GDS-4 Short-Short - Mean-substituted (truncated)"
label var wb_gds_score_meansubR  "GDS-15 Short - Mean-substituted (rounded)"
label var wb_gds_score_X15R      "GDS-15 Short - Row mean X 15 items"
label var wb_gds_score2          "GDS Short-Short Form (GDS-4) - full case"

* depression indicators
gen GDS15_depr=1 if wb_gds_score_X15R>4 & wb_gds_score_X15R~=.
replace GDS15_depr=0 if wb_gds_score_X15R<5 & wb_gds_score_X15R~=.
gen GDS4_depr=1 if wb_gds_score_X4R>1 & wb_gds_score_X4R~=.
replace GDS4_depr=0 if wb_gds_score_X4R<2 & wb_gds_score_X4R~=.

* fill reversed items where missing (for raw-item INR/MbD logic)
replace wb_gds1r=wb_gds1 if missing(wb_gds1r)
replace wb_gds5r=wb_gds5 if missing(wb_gds5r)
replace wb_gds7r=wb_gds7 if missing(wb_gds7r)
replace wb_gds11r=wb_gds11 if missing(wb_gds11r)
replace wb_gds13r=wb_gds13 if missing(wb_gds13r)

* ---- id_particip_new: fill 83 missing as L0, then manual reclassifications --*
gen id_particip_new=id_particip
label values id_particip_new id_particip_
replace id_particip_new=1 if missing(id_particip)

* MbD flag (person-level GDS-15): structurally missing item 1 among L0/L1
gen GDS_mbd=0
replace GDS_mbd=1 if wb_gds1r==. & id_particip_new==1
replace GDS_mbd=1 if wb_gds1r==. & id_particip_new==2


* ---- FIELD TEAM: fill missing id_team from canton --------------------------*
* Each canton belongs to exactly one field team; 83 participants have id_team
* missing, but their canton is known.
capture drop team_fill
bysort id_canton: egen team_fill = mode(id_team), minmode
count if !missing(id_team) & id_team != team_fill
assert r(N)==0                                  // no canton spans two teams
count if missing(id_team)
di as txt "id_team missing before fill: " r(N) "   (expected 83)"
replace id_team = team_fill if missing(id_team)
assert !missing(id_team)
tab id_team, missing                            // expected 100 / 97 / 80, total 277

* per-item MbD flags
foreach var of varlist wb_gds3 wb_gds4 wb_gds8 wb_gds14 wb_gds1r wb_gds2 wb_gds5r wb_gds6 wb_gds7r wb_gds9 wb_gds10 wb_gds11r wb_gds12 wb_gds13r wb_gds15 {
    gen `var'_mbd=0
    replace `var'_mbd=1 if `var'==.
}
gen GDS_mbd_v2=0
foreach var of varlist wb_gds3 wb_gds4 wb_gds8 wb_gds14 wb_gds1r wb_gds2 wb_gds5r wb_gds6 wb_gds7r wb_gds9 wb_gds10 wb_gds11r wb_gds12 wb_gds13r wb_gds15 {
    replace GDS_mbd_v2=1 if `var'_mbd==1
}
gen GDS15_mbd=GDS_mbd
drop GDS_mbd_v2 GDS_mbd
label var GDS15_mbd "GDS-15 missing-by-design (excludes L0,L1)"

gen GDS4_mbd=0
replace GDS4_mbd=1 if id_particip_new==1
label var GDS4_mbd "GDS-4 missing-by-design (excludes L0)"
* ---- INR (spontaneous) person + item level --------------------------------*
* Explicit item lists (no varlist ranges: ranges depend on variable order).
local items15 wb_gds1r wb_gds2 wb_gds3 wb_gds4 wb_gds5r wb_gds6 wb_gds7r ///
              wb_gds8 wb_gds9 wb_gds10 wb_gds11r wb_gds12 wb_gds13r wb_gds14 wb_gds15
local items4  wb_gds3 wb_gds4 wb_gds8 wb_gds14

gen GDS15_inr_any=0
foreach var of local items15 {
    replace GDS15_inr_any=1 if `var'==.a | `var'==.b | `var'==.c | `var'==.d
}
replace GDS15_inr_any=. if GDS15_mbd==1

local inr15list
foreach var of local items15 {
    gen `var'_inr=0
    replace `var'_inr=1 if `var'==.a | `var'==.b | `var'==.c | `var'==.d
    replace `var'_inr=. if `var'_mbd==1
    local inr15list `inr15list' `var'_inr
}
egen GDS15_inr_cnt=rowtotal(`inr15list')
replace GDS15_inr_cnt=. if GDS15_mbd==1

gen GDS4_inr_any=0
foreach var of local items4 {
    replace GDS4_inr_any=1 if `var'==.a | `var'==.b | `var'==.c | `var'==.d
}
replace GDS4_inr_any=. if GDS4_mbd==1
egen GDS4_inr_cnt=rowtotal(wb_gds3_inr wb_gds4_inr wb_gds8_inr wb_gds14_inr)
replace GDS4_inr_cnt=. if GDS4_mbd==1

* ---- eligibility flags and labels -------------------------------------------*
* id_particip_new is coded 1=L0, 2=L1, 3=L2, 4=L3
gen byte elig15 = inlist(id_particip_new,3,4)      // L2-L3
gen byte elig4  = inlist(id_particip_new,2,3,4)    // L1-L3

label var sd_age_OFS            "Age, years"
label var sd_gender_OFS         "Sex"
label var sd_living_private_xcp "Living arrangement"
label var Frailty_cat           "Frailty"
label var cg_mmse_totscore      "MMSE-21"
label var id_team               "Field team (linguistic region)"
label var GDS15_inr_any         "Any INR in GDS-15 (eligible)"
label var GDS4_inr_any          "Any INR in GDS-4 (eligible)"


*##############################################################################*
* SECTION 0 - FREEZE THE FLOW-CHART COUNTS  (Figure 1)
* Each assert stops the file if a count differs from the manuscript.
*##############################################################################*

di as txt "===== FLOW-CHART COUNTS (V7) ====="
count
assert r(N)==277
tab id_particip_new, missing                      // 107 / 38 / 29 / 103
count if id_particip_new==1
assert r(N)==107
count if elig4
assert r(N)==170                                  // GDS-4 eligible (L1-L3)
count if elig15
assert r(N)==132                                  // GDS-15 eligible (L2-L3)

/*Participation level for depression assessment reflects each centenarian's
capacity to self-report at the time of the GDS, and may differ slightly from
participation classifications used elsewhere in the SWISS100 cohort*/

tab GDS15_inr_any, missing                        // 73 no INR / 59 any INR
count if GDS15_inr_any==0
assert r(N)==73
count if GDS15_inr_any==1
assert r(N)==59
tab GDS4_inr_any, missing                         // 142 no INR / 28 any INR
count if GDS4_inr_any==0
assert r(N)==142
count if GDS4_inr_any==1
assert r(N)==28

di as txt "--- Score availability (analytic samples) ---"
count if !missing(wb_gds_score_X15R)              // 165 (L1-L3, L1 extrapolated)
assert r(N)==165
count if !missing(wb_gds_score_X15R) & elig15     // 130 = GDS-15 analytic sample
assert r(N)==130
count if !missing(wb_gdsSS_score_meansub)         // 165 = GDS-4 analytic sample
assert r(N)==165
count if !missing(wb_gds_score)                   // 73 complete cases
assert r(N)==73


*##############################################################################*
* TABLE 1 - Baseline characteristics, overall and by participation level
*##############################################################################*
dtable, ///
    by(id_particip_new, tests total) ///
    continuous(sd_age_OFS cg_mmse_totscore, statistics(mean sd)) ///
    factor(sd_gender_OFS sd_living_private_xcp Frailty_cat id_team, ///
           statistics(fvfrequency fvpercent)) ///
    sformat("(%s)" sd) ///
    title("Table 1. Baseline characteristics by participation level") ///
    export("$out/Table1_baseline_byLevel.docx", replace)
tab id_team id_particip_new, col                  // VD 29/12/14/45, ZH 36/12/7/42, TI 42/14/8/16


*##############################################################################*
* TABLE 2 - INR per GDS-15 item, by participation level + Total
* L1 (id_particip_new==2) was administered only items 3, 4, 8, 14.
* Item order = administration order 1..15.
*##############################################################################*
local order wb_gds1r wb_gds2 wb_gds3 wb_gds4 wb_gds5r wb_gds6 wb_gds7r ///
            wb_gds8 wb_gds9 wb_gds10 wb_gds11r wb_gds12 wb_gds13r wb_gds14 wb_gds15

capture frame drop t2
frame create t2 str12 item int pos int inrL1 int nL1 int inrL2 int nL2 int inrL3 int nL3
local pos 0
foreach v of local order {
    local ++pos
    local row ("`v'") (`pos')
    foreach lev in 2 3 4 {                        // 2=L1, 3=L2, 4=L3
        qui count if `v'_inr==1 & id_particip_new==`lev'
        local a = r(N)
        qui count if `v'_inr==0 & id_particip_new==`lev'
        local b = r(N)
        local row `row' (`a') (`=`a'+`b'')
    }
    frame post t2 `row'
}

frame t2 {
    gen inrLT = inrL1 + inrL2 + inrL3
    gen nLT   = nL1 + nL2 + nL3
    gen inrL23 = inrL2 + inrL3
    gen nL23   = nL2 + nL3
    gen pct_L1 = 100*inrL1/nL1  if nL1>0
    gen pct_L2 = 100*inrL2/nL2
    gen pct_L3 = 100*inrL3/nL3
    gen pct_L23 = 100*inrL23/nL23
    gen pct_LT = 100*inrLT/nLT
    format pct_* %4.1f

    di as txt _n "=== Table 2: INR per item (n INR / n administered, %) ==="
    list item pos inrL1 nL1 pct_L1 inrL2 pct_L2 inrL3 pct_L3 inrLT nLT pct_LT, sep(0) noobs

    * ---- overall INR rates (fixes the 7.2% of V6) ----
    foreach s in inrL1 nL1 inrL2 nL2 inrL3 nL3 {
        qui summ `s'
        scalar S_`s' = r(sum)
    }
    di as txt _n "L1      : " S_inrL1 " of " S_nL1 " = " %4.1f 100*S_inrL1/S_nL1 "%   (expected 20/152 = 13.2%)"
    di as txt    "L2      : " S_inrL2 " of " S_nL2 " = " %4.1f 100*S_inrL2/S_nL2 "%   (expected 67/435 = 15.4%)"
    di as txt    "L3      : " S_inrL3 " of " S_nL3 " = " %4.1f 100*S_inrL3/S_nL3 "%   (expected 103/1545 = 6.7%)"
    di as txt    "L2-L3   : " (S_inrL2+S_inrL3) " of " (S_nL2+S_nL3) " = " %4.1f 100*(S_inrL2+S_inrL3)/(S_nL2+S_nL3) "%   (expected 170/1980 = 8.6%)"
    di as txt    "L1-L3   : " (S_inrL1+S_inrL2+S_inrL3) " of " (S_nL1+S_nL2+S_nL3) " = " %4.1f 100*(S_inrL1+S_inrL2+S_inrL3)/(S_nL1+S_nL2+S_nL3) "%   (expected 190/2132 = 8.9%, Table 2 Total)"

    * ---- Pearson r: item position vs item INR rate ----
    di as txt _n "=== Pearson r, item position vs INR rate ==="
    qui pwcorr pos pct_L1
    di as txt "L1 (4 items)     : r = " %5.3f r(rho) "   (expected 0.818)"
    qui pwcorr pos pct_L2
    di as txt "L2               : r = " %5.3f r(rho) "   (expected 0.629)"
    qui pwcorr pos pct_L3
    di as txt "L3               : r = " %5.3f r(rho) "   (expected 0.688)"
    qui pwcorr pos pct_L23
    di as txt "L2-L3 combined   : r = " %5.3f r(rho) "   (expected 0.762)"
    qui pwcorr pos pct_LT
    di as txt "Total (L1-L3)    : r = " %5.3f r(rho) "   (expected 0.818)"
    spearman pos pct_L23, stats(rho p)

    export excel using "$out/Table2_INR_by_item.xlsx", firstrow(variables) replace
}
frame drop t2

* person-level INR by level (rows "Any INR in GDS-15")
tab GDS15_inr_any id_particip_new if elig15, col     // L2 16 (55.2%), L3 43 (41.7%), Total 59 (44.7%)


*##############################################################################*
* TABLE 3 - Determinants of GDS-15 INR (eligible L2-L3)
* Binary depression (ipsative GDS-15 >=5) in ALL four models.
*##############################################################################*
eststo clear
eststo m1: logit GDS15_inr_any i.GDS15_depr c.cg_mmse_totscore if elig15, or
eststo m2: logit GDS15_inr_any i.GDS15_depr c.cg_mmse_totscore c.sd_age_OFS ///
    i.sd_gender_OFS i.id_team i.sd_living_private_xcp i.Frailty_cat if elig15, or
eststo m3: nbreg GDS15_inr_cnt i.GDS15_depr c.cg_mmse_totscore if elig15, irr
eststo m4: nbreg GDS15_inr_cnt i.GDS15_depr c.cg_mmse_totscore c.sd_age_OFS ///
    i.sd_gender_OFS i.id_team i.sd_living_private_xcp i.Frailty_cat if elig15, irr

esttab m1 m2 m3 m4 using "$out/Table3_determinants.rtf", replace ///
    eform b(2) ci(2) star(* 0.05 ** 0.01) nonumbers ///
    stats(N, labels("Number of observations")) ///
    mtitles("Logit OR (unadj)" "Logit OR (adj)" "NegBin IRR (unadj)" "NegBin IRR (adj)") ///
    title("Table 3. Determinants of GDS-15 item non-response (L2-L3)") ///
    addnotes("OR/IRR >1 = higher likelihood or count of non-response. Depression = ipsative GDS-15 >=5. Ref: no depression; male; VD; institution; robust.")
esttab m1 m2 m3 m4, eform b(2) ci(2) star(* 0.05 ** 0.01) stats(N) ///
    mtitles("OR unadj" "OR adj" "IRR unadj" "IRR adj")
* Expected (manuscript): depression OR 1.39 / 2.28 / IRR 1.35 / 2.09*; MMSE 0.88* / 0.86* / 0.87** / 0.89*; N 130 / 129 / 130 / 129


*##############################################################################*
* TABLE 4 (descriptive) - GDS scores by scoring method (mean, SD, n)
* Final variables used in the manuscript text.
*##############################################################################*
di as txt "=== GDS-15 scores by method ==="
tabstat wb_gds_score wb_gds_8i wb_gds_score_meansubR wb_gds_score_X15R, ///
    stat(mean sd n) col(stat) format(%6.2f)
* complete case 3.77 (2.54) n=73 | Rule-8 3.76 (2.36) n=127 | ipsative 3.75 (2.94) n=165
* sample-mean (rounded, meansubR): value to be inserted in the manuscript

di as txt "=== GDS-4 scores by method ==="
tabstat wb_gds_score2 wb_gds_3i wb_gdsSS_score_meansub wb_gds_score_X4R, ///
    stat(mean sd n) col(stat) format(%6.2f)
* complete case 0.68 (0.96) n=142 | Rule-3 0.65 (0.93) n=158 | sample-mean 0.68 (0.98) n=165
* ipsative: value to be inserted in the manuscript (manuscript currently 0.67 (0.95))

* Not used in the manuscript, for comparison only
tabstat wb_gds_score_meansub wb_gds_score_meansub2 wb_gdsSS_score_meansub2 wb_gdsSS_score_meansubR, ///
    stat(mean sd n) col(stat) format(%6.2f)

* Complete cases: raw sum vs rounded ipsative on the same 73 people
tabstat wb_gds_score wb_gds_score_X15 wb_gds_score_X15R if !missing(wb_gds_score), ///
    stat(mean sd n) col(stat) format(%6.3f)

di as txt "=== Depression prevalence (ipsative) ==="
tab GDS15_depr, missing                            // expected 61 of 165 = 37.0%
tab GDS15_depr if elig15, missing                  // GDS-15 analytic sample (n = 130)
tab GDS4_depr, missing                             // expected about 18.2%
tab GDS15_depr id_particip_new, col

*==============================================================================*
* END OF V7
*==============================================================================*
