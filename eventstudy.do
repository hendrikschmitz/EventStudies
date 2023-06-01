
cap program drop eventstudy
program define eventstudy, eclass
version 15.0

syntax [varlist(fv)] [if] [in] , eventvariable(varlist) reference(numlist)  [Absorb(varlist fv) by(varlist) weighted(varlist) cluster(varlist)  keinbinning_pre keinbinning_post  /// 
	window_lower(real -100) window_upper(real -100) window_lower_graph(real -100) window_upper_graph(real -100) iv_endog(varlist) iv_instrument(varlist) true(varlist) /// 
	gr_opts(str) gr_title(str) scproperties(str) rcapproperties(str)  level(real 0.05) supp_reg_output  graphresultssave(str) sunabraham cohort(varlist) control_cohort(varlist) ///
	interact(varlist) interact_missing event_dummies_manuell(str) supp_figure]


	
quie {
	
	numlist "`reference'"
	local nn "`r(numlist)'"

	********************************************************
	* Vorschritt: zunächst die tatsächlichen Werte speichern
	********************************************************

	if "`true'" != "" {
		preserve
		collapse `true', by(`eventvariable')
		tempfile true_dat
		save `true_dat'
		restore
	}

 

	************************************************************
	* Eigentlicher Start: Variablen einlesen und generieren usw.
	************************************************************
	
	preserve
	
	if "`if'" != "" keep `if'

	marksample touse
	gettoken depvar varlist : varlist

	loc eventvariable: list varlist | eventvariable
	loc eventvariable: list eventvariable - varlist
	global xvars `varlist' 

	if "`cluster'" != "" {
		loc cluster: list varlist | cluster
		loc cluster: list cluster - varlist
		loc cluster_reg "cluster(`cluster')"
		loc cluster_sa "vce(cluster `cluster')"
	}


	loc ols reg
	loc abs
	if "`absorb'" != "" {
		loc abs "absorb(`absorb')"
		loc ols reghdfe
	}
	
		
	if "`sunabraham'" == "sunabraham" & ("`abs'" == "" | "`control_cohort'" == "" | "`cohort'" == "" | "$xvars" == "") {
	no di in red "Sun and Abraham requies controls, cohort, control cohort and absorb!"
	exit
	}

	
	***************************
	* Schleife je Sumbsample
	***************************

	if "`by'" != ""  {
		tempvar sample
		egen `sample' = group(`by')
		sum `sample'
		local subsamples = r(max)
	}
	else local subsamples = 1

	tempfile alldataeventstudy
	save `alldataeventstudy'

	*************************************************************************************************************************************************
	forvalues subsample = 1/`subsamples' {
		use `alldataeventstudy', clear
		if `subsamples' != 1	local if_subsample if `sample' == `subsample'

		* Default Effect-Window: Minimum und Maximum der Eventvariable
		quie sum `eventvariable'
		local max = r(max)
		local min = r(min)
		if `window_lower' == -100 local window_lower `min'
		if `window_upper' == -100 local window_upper `max'

		* Event-variablen im Effect-Fenster generieren
		forvalues j = `window_lower'/-1 {
			local k = -1*`j'
			gen e_min`k' = `eventvariable' ==`j'
		}
		forvalues j = 0/`window_upper' {
			gen e`j' = `eventvariable' == `j'
		}	
	

	
		local i = `window_upper'	
		local h = `window_lower'
		*  Für das Binning eine weitere Kategorie erstellen
		if `window_upper' < `max' &  "`keinbinning_post'" == ""  {
			local i = `window_upper' +1
			gen e`i' = `eventvariable' > `window_upper'
			replace e`i' = 0 if `eventvariable' == .
		}

		if `window_lower' > `min' & "`keinbinning_pre'" == ""  {	
			local h = `window_lower' -1
			local hh = -1*`h'
			gen e_min`hh' = `eventvariable' < `window_lower' 
		}

	
		**************************
		* Global mit Eventdummies
		**************************	
		
		global event_vars
		forvalues j = `h'/`i' {
			if `: list j in nn'==0 {    /*`reference2' == -100 als default, wenn zweiter Wert nicht angegeben) */
				local k = -1*`j'
				if `j' < 0 & `j' >= `min'  global event_vars $event_vars e_min`k'  /* Hier auch sicherstellen, dass Event-dummies außerhalb von min und max der Eventvariablen nicht mit aufgenommen werden. */
				if `j' >= 0 & `j' <= `max'  global event_vars $event_vars e`j'
			}
		}
		
	
		if "`interact'" != "" {
			tempvar intvar
			egen `intvar' = group(`interact')
			sum `intvar'
			local n_intvar = r(max)
			tab `intvar', gen(intvar)
			
			
			global event_vars_int
			forvalues cat = 1/`n_intvar'	{
				if "`interact_missing'" == "interact_missing" replace intvar`cat' = 0 if intvar`cat' == .    /* Missings auf 0 setzen. */
				foreach var of varlist $event_vars {
					gen `var'_`cat' = `var' * intvar`cat'
					global event_vars_int $event_vars_int `var'_`cat' 
					}
				}
			global event_vars $event_vars_int
		}	
	
		if "`event_dummies_manuell'" != "" global event_vars `event_dummies_manuell'	
		
		
		*************************
		* Regression
		*************************
		
		if "`supp_reg_output'"=="" local no no

		local weight_reg
		if "`weighted'" != "" local weight_reg [aw = `weighted']

	
		if "`sunabraham'" != "sunabraham" {
			if "`iv_endog'" == "" {
				noi di "`no' `ols'	`depvar' $event_vars $xvars `weight_reg' `if_subsample', `cluster_reg' `abs'"
						`no' `ols'	`depvar' $event_vars $xvars `weight_reg' `if_subsample', `cluster_reg' `abs'			
			}

			if "`iv_endog'" != "" {
				noi di "`depvar' (`iv_endog' = `iv_instrument') $event_vars $xvars `weight_reg' `if_subsample', `cluster_reg'"
				`no' ivregress 2sls 	`depvar' (`iv_endog' = `iv_instrument') $event_vars $xvars `weight_reg' `if_subsample', `cluster_reg'
			}
		}
		
		if "`sunabraham'" == "sunabraham" {
			eventstudyinteract `depvar' $event_vars, cohort(`cohort') control_cohort(`control_cohort') covariates($xvars) `abs' `cluster_sa'
		}

		*est store results
		*`no'  esttab results, se(3) label star(* 0.1 ** 0.05 *** 0.01) b(3)  lines wide replace nogaps title(`depvar') keep($event_vars)


		
		* Ein paar Infos damit man auch ohne Regressionsoutput versteht, was gemacht wurde:
		count if e(sample)
		local obs = r(N)
		count if e(sample) & `eventvariable' == . 
		local obs_control = r(N)

		noi di "Information on the event study regressions:"
		noi di "Reference category: `nn'"
		no display "Number of observations in the regression: " `obs'
		
		if `obs_control' != 0 {
			no display `obs_control' " observations did not receive a treatment and are a pure control group."
			no display "They enter the regressions with all event-time indicators set to 0."
		}


		*************************
		* Grafik
		*************************

		
		mat b = e(b)
		if "`interact'" != "" {
			clear
			set obs 1000
			mat b2 = e(b)'
			svmat2 b2, name(bvector) r(Bnames)
			save "`graphresultssave'", replace
			restore
			no di ""
			no di "No graph because of interactions."
			quie exit
		}
		
		mat V = e(V)
		if "`sunabraham'" == "sunabraham" {
			mat b = e(b_iw)
			mat V = e(V_iw)	
		}

		local rows = (-1*`h') + (`i') + 1  /* alle negativen und positiven + die 0  */
		mat def results = J(`rows',3,.)

		local k = 0  /* Zähler für jede Zeile der zu erstellenden Ergebnismatrix (inklusive Referenzgruppen) */
		local kk = 0  /* Zähler für jede Zeile des genutzten Regressionsoutputs (ohne Referenzgruppen) */ 
		if "`iv_endog'" != "" local kk = 1 /* Endogene Variable bei IV wird oben angezeigt, verschiebt dann die Eventdummies */
		forvalues j = `h'/`i' {
			local k = `k' + 1
			mat results[`k',1] = `j'

			if `: list j in nn'==0 {
				local kk = `kk' + 1
				mat results[`k',2] = b[1,`kk']
				mat results[`k',3] = sqrt(V[`kk',`kk'])
				if "`sunabraham'" == "sunabraham" mat results[`k',3] = sqrt(V[1,`kk'])
			}

			if `: list j in nn'==1 {
				mat results[`k',2] = 0
				mat results[`k',3] = 0
			}
		}

		clear
		set obs 1
	
		svmat results
		rename results1 `eventvariable'
		rename results2 `depvar'
		gen ciu = `depvar' - `=invnormal(1-`level'/2)'*results3
		gen cio = `depvar' + `=invnormal(1-`level'/2)'*results3
		drop results3


		if "`true'" != "" {
			sort `eventvariable'
			merge `eventvariable' using `true_dat' 
			drop _merge
			local grafiktrue (line `true' `eventvariable', lwidth(thick) lpattern(dash) lcolor(red))  
		}


		* Default ist Länge Effect window.
		if `window_lower_graph' == -100  local window_lower_graph `window_lower'
		if `window_upper_graph' == -100  local window_upper_graph `window_upper'

		drop if `eventvariable' < `window_lower_graph' |  `eventvariable' > `window_upper_graph'


		if `subsample' > 1 {	
			local gap = (`subsample'-1)/10*2
			replace `eventvariable'=`eventvariable'+`gap'

			rename `eventvariable' `eventvariable'`subsample'
			rename `depvar' `depvar'`subsample'
			rename ciu ciu`subsample'
			rename cio cio`subsample'
			gen sample = `subsample'

		}
		tempfile graphdata`subsample'
		save `graphdata`subsample''


	* by()
	}
	
	*************************************************************

	if `subsamples' > 1 {
		use `graphdata1'
		forvalues subsample = 2/`subsamples' {
			append using `graphdata`subsample''
		}	
	}

	*if "`scproperties'"=="" local scproperties "mcolor(blue) msize(medium) msymbol(circle)"	
	
	/*
	if	"`scproperties'"=="" local scproperties "mcolor(black*0.6) msize(large) msymbol(triangle)"	
	local scproperties2 "mcolor(black*0.8) msize(large) msymbol(square)"	
	local scproperties3 "mcolor(black)     msize(large) msymbol(diamond)"	
	*/
	if	"`scproperties'"=="" local scproperties "col(navy) msymb(Oh)"	
	local scproperties2 "color(black)  msymbol(Sh)"	
	local scproperties3 "mcolor(black*0.8)   msymbol(diamond)"
	

	if	"`rcapproperties'"=="" local rcapproperties "col(navy)"	
	local rcapproperties2 "lcolor(black)"	
	local rcapproperties3 "lcolor(black*0.8)"	

	if `subsamples' == 1 local scatter (scatter `depvar' `eventvariable', sort `scproperties') (rcap ciu  cio `eventvariable', `rcapproperties') 
	
	if `subsamples' == 2 local scatter (scatter `depvar' `eventvariable', sort `scproperties') (rcap ciu  cio `eventvariable', `rcapproperties')  ///  ///
			|| (scatter `depvar'2 `eventvariable'2, sort `scproperties2')   (rcap ciu2  cio2 `eventvariable'2, `rcapproperties2')

	if `subsamples' == 3 local scatter (scatter `depvar' `eventvariable', sort `scproperties') (rcap ciu  cio `eventvariable', `rcapproperties')  ///
		|| (scatter `depvar'2 `eventvariable'2, sort `scproperties2')  (rcap ciu2  cio2 `eventvariable'2, `rcapproperties2') ///
		|| (scatter `depvar'3 `eventvariable'3, sort `scproperties3')  (rcap ciu3  cio3 `eventvariable'3, `rcapproperties3') 	
	
	la var `eventvariable' "t"
	cap la var `eventvariable'2 "`by'"
	
	* Default xline
	if strmatch("`gr_opts'", "*xline*") != 1 local xline xline(-0.5, lwidth(medium) lpattern(dash) lcolor(black))
	if strmatch("`gr_opts'", "*xtitle*") != 1 local xtitle xtitle(Event time , size(medium))
	if strmatch("`gr_opts'", "*xlabel*") != 1 local xlabel xlabel(`window_lower_graph'(1)`window_upper_graph', labsize(medium)) xscale(range(`window_lower_graph' (1) `window_upper_graph')) 
	
	if "`supp_figure'" == "" twoway `scatter' `grafiktrue', scheme(s1mono) `gr_opts' ytitle(Effect) yline(0, lwidth(medium) lpattern(dash) lcolor(black))  `xline' `xlabel' `xtitle' `gr_title'
	
	
	
	* Grafikergebnisse speichern
	if "`graphresultssave'" != "" save "`graphresultssave'", replace		
		
	restore
}

end


