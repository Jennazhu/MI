
/*
cd ../..
*/

clear all
global atr "1-data/atrophy"




*** import  data
	use "$atr/dat_wide.dta", clear
	codebook id
	codebook id if !missing(atrophy)
  
  *rename time-invariant (baseline) vars to make cross-temporal MI spec easier to read 
  rename (age  htn  smoke  drink  dm  apoe  educ)  ///
         (age1 htn1 smoke1 drink1 dm1 apoe1 educ1)
  rename (atrophy  male  siterace  wmh)  ///
         (atrophy1 male1 siterace1 wmh1)
  
  
********************************************************************************
*Step 1:  Declare multiple-imputation data
  mi set wide

  
  *examine missingness patterns
  mi misstable summarize globz1 globz2 globz3 globz4 globz5 // summerize missingness of each variable
  mi misstable patterns  globz1 globz2 globz3 globz4 globz5, asis freq // summerize missingness of each variable
  

********************************************************************************
*Step 2: "Register" the missing variables


** Step 2.1: register imputed variables
** All baseline variables should be here (everyone is alive at baseline)
  global impvars ///
          globz1 globz2 globz3 globz4 globz5  ///
	      age1 htn1 smoke1 drink1 dm1 apoe1 educ1 // v2globz
	 
	 * summarize impvars
     summ idnum $impvars
	 * registeer
	  mi register imputed  $impvars

	   
** Step 2.2: register regular variables
	* Variables that do NOT contain missing data should be included after the 
	* "mi register regular" statement 
  global regvars ///
	     atrophy1 male1 siterace_JB1 siterace_FB1 siterace_FW1  wmh1  ///
	     alive1 alive2 alive3 alive4 alive5 ///  
	     dem1 dem2 dem3 dem4 dem5  
	   
	summ idnum $regvars
	mi register regular $regvars
	   
	   

********************************************************************************
*Step 3 : imputation

    * xtset for mi data
	capture mi xtset, clear
    
	*create inclusion variable statements for each visit (here 1 to 5)
  global incl0 // blank, nothing earlier to include than v1
	forvalues i = 1(1)5 {
    local j = `i'-1 // used in the global incl statements
    global incl`i'  (atrophy1*time`i') (atrophy1*male1) (atrophy1*htn1) (atrophy1*educ1) ///
                    (atrophy1*apoe1) (atrophy1*siterace_JB1)  (atrophy1*siterace_FB1)  (atrophy1*siterace_FW1) ///
                    (atrophy1*c.age1) (atrophy1*drink1) ///
                    (dem`i'*time`i') (alive`i'*time`i') ///
                    ${incl`j'}
    di
    di "Global macro var incl" `i' ": " "${incl`i'}"
	}

	*create omitted variable statements for each visit (here 1 to 4)
  *omit anything from a future time
  global omit5 // blank, nothing future to omit from last visit but needed for loop
  forvalues i=5(-1)1 {
    local j = `i'-1 // used in the global omit name
    global omit`j' alive`i' time`i' dem`i'  globz`i' ${omit`i'}
    di
    di "Global macro var omit" `i' ": " "${omit`i'}"
  }
  

  

  * global seed for imputation
  global mseed = 042020 // set MI random seed (change this every once in awhile)
  
  *choose which types of death information to include:
  *global keepdeathcat = 2 // don't use the visit if they died during
  *global keepdeathcat = 1 // don't use the visit after death if they died in-between visits
  global keepdeathcat = 0 // keep an observation in the dataset if they died in-between visits
	* start imputation
	mi impute chained ///
	   (regress, omit( $omit1) ) age1   /// baseline normally distb covs /// 
	   (logit,   omit( $omit1) ) htn1  dm1  /// baseline binary covs
	   (ologit,  omit( $omit1) ) apoe1 educ1 drink1 smoke1  /// baseline ordinal covs
	   (regress if alivecat1 > $keepdeathcat , incl( $incl1 ) omit( $omit1) ) globz1 ///
	   (regress if alivecat2 > $keepdeathcat , incl( $incl2 ) omit( $omit2) ) globz2 ///
	   (regress if alivecat3 > $keepdeathcat , incl( $incl3 ) omit( $omit3) ) globz3 ///
	   (regress if alivecat4 > $keepdeathcat , incl( $incl4 ) omit( $omit4) ) globz4 ///
	   (regress if alivecat5 > $keepdeathcat , incl( $incl5 ) omit( $omit5) ) globz5 ///
		  = atrophy1 male1 i.siterace1 wmh1 ///
        time1 time2 time3 time4 time5  ///
			  alive1 alive2 alive3 alive4 alive5 ///
  			dem1 dem2 dem3 dem4 dem5 /// 
		   , add(30) burnin(50) rseed($mseed) report force dots augment 

    *check:
	forvalues i = 1(1)5 {
    tab alivecat`i' if alivecat`i' > 0
	}
       
	  save "$atr/impute/midat", replace


	*****************************
	** change it into long format
	******************************

	use "$atr/impute/midat", clear

	
	* reshape the data
	mi reshape long vdate time alive alivecat globz , i(idnum) j(visit) 

	* xtset the data
	mi xtset idnum visit
    
	
	save "$atr/impute/midat_long", replace
  
  
  