*! psacalc


program psacalc2, rclass

	version 11.2
	syntax anything [, mcontrol(string) rmax(real 1.0) delta(real 1.0) beta(real 0) model(string)]

	* Get treatment var and type of calculation

	tokenize `anything'
	local type "`1'"
	local treatment "`2'"

	if "`3'"!="" {
		di as err _n "Too many variables"
		exit 103
	}

	tempvar error_hat treat_res touse

	tempname regmodel beta_tilde r_tilde beta_o r_o sigma_yy sigma_xx t_x scbeta scdelta scrmax

	* If model is specified, estimate
	if "`model'"!="" {
		qui `model'
	}

	* If model is not specified, uses previous model
	* Error message if not found
	else {
		 if "`e(cmd)'" == "" {
			di as err _n "Previous estimation results not found"
			exit 301
		}
	}

	* Get command
	loc command = e(cmd)

	* Mark sample from if in in previous model
	gen `touse'=e(sample)
	* Get weights if previous regression had them
	if "`e(wtype)'"!="" {
		loc weight "`=e(wtype)' `=e(wexp)'"
	}

	* If xtreg, get type
	if "`command'"=="xtreg" loc fetype=e(model)

	* If areg, get absorbed
	if "`command'"=="areg" {
		local absorb=e(absvar)
		local absorb = "absorb(`absorb')"
	}

		* If reghdfe, get absorbed
	if "`command'"=="reghdfe" {
		local absorb=e(absvars)
		local absorb = "absorb(`absorb')"
    local residuals= "residuals" // reghdfe requires this option to post est 'predict, d'
	}


	est store `regmodel'

	/*=========================================================================
						1: Capture errors
	===========================================================================*/

	local names: colnames e(b)
	loc words: word count `names'

	// Delta ignored if requesting delta and entered delta!=1
	if "`type'"=="delta" & `delta'!=1 di as error "User entered delta ignored"

	// Beta ignored if requesting beta and entered beta!=1
	if "`type'"=="beta" & `beta'!=0 di as error "User entered beta ignored"

	// Does not have a constant
	if strpos("`names'","_cons")==0  {
		di as error "Model does not have a constant"
		exit 5001
	}

	// No controls TODO add check for FE
	if `words'==2 & "`absorb'"=="" /// fail if no controls or absorbed effects
		& strpos("`names'","_cons")!=0 & strpos("`names'","`treatment'")!=0  {
		di as error "Model has neither controls nor absorbed effects treated as controls"
		exit 5001
	}


	// No previous regression, and failed to specify linear model
	if !inlist(e(cmd),"regress","xtreg","areg","reghdfe") {
		di as error _n "This command only works after linear regression"
		di as error _n "Estimate a linear regression before running this command"
		di as error _n "Or specify a full linear regression model in options"
		exit 5001
	}

	// Do not allow xtreg, re , xtreg, pa or xtreg, mle
	if e(cmd)=="xtreg" & !inlist("`fetype'","fe","be") {
		di as error _n "Only fe or be allowed for xtreg"
		exit 5001
	}

	// rmax out of bounds
	if `rmax'>1 {
		di as err _n "The maximum possible R-squared is 1"
		exit 5001
	}

	scalar `beta_tilde'=_b[`treatment']
	* rsq from xtreg, fe is also e(r2), within
	* rsq from xtreg, be is also e(r2), between
	scalar `r_tilde'=e(r2)

	// rmax less than controlled R-squared
	if `rmax'<e(r2) {
		di as err _n "Maximum R-squared provided is less than controlled R-squared"
		exit 5001
	}

	// bad type
	if "`type'"~="beta" & "`type'"~="delta" {
		di as err _n "Invalid type: must be beta or delta"
		exit 5001
	}

	// specified treatment not in the model

	if strpos("`names'","`treatment'")==0 {

		// treatment is not in the model at all
		if strpos(e(cmdline),"`treatment'")==0 {
			di as err _n "Specified treatment variable missing from model"
			exit 5001
		}

		// treatment is not an independent variable in the model
		else {
			di as err _n "Specified treatment variable is not an independent variable in the model"
			exit 5001
		}
	}

	/// treatment in mcontrols
	if "`mcontrol'"!="" {
		loc mclist="`mcontrol'"
		loc tlist ="`treatment'"
		if "`mcontrol'"!="" {
			loc mcmt: list mclist - tlist
			if "`mcmt'"!="`mcontrol'" {
				di as err "Treatment can not be an unrelated control"
				exit 5001
			}
		}

		/// mcontrols not in initial controls
		loc clist="`names'"
		loc cmmc: list clist - mclist
		if "`cmmc'"=="`clist'" {
			di as err "Unrelated control not in regression"
			exit 5001
		}
	}

	/*=========================================================================
						2: Main program arc
	===========================================================================*/

	*------------------------ 2.1: Get remaining inputs from regressions ----------------------------------



	loc depvar=e(depvar)
	loc indepvars: colnames e(b)
	loc cons "_cons"
	loc indepvars: list indepvars - cons

	if "`command'"=="areg" | "`command'"=="reghdfe" loc command2 = "reg" // treat "areg" and "reghdfe" the same
	else loc command2="`command'"

	quietly `command2' `depvar' `treatment' `mcontrol' [`weight'] if `touse', `fetype'

	scalar `beta_o'=_b[`treatment']
	scalar `r_o'=e(r2)

	// Define variance Terms

	quietly sum `e(depvar)' if `touse'
	scalar `sigma_yy'=r(Var)
	* Residualize X on m before taking this variance

	if "`command'"=="xtreg" & ("`fetype'"=="fe") loc res="e"
	else if "`command'"=="xtreg" & ("`fetype'"=="be") loc res="ue"
	else loc res "resid"

	* Check if treatment is factor variable, build if necessary


	cap reg `treatment'
	if _rc==198	{
		tempvar treat2
		gen `treat2'=`treatment'
	}
	else {
		loc treat2="`treatment'"
	}

	* Weigths go here because they affect the calculation of the mean
	* Fixed effects go here because we are using the variance of the treatment in the differenced regression

	quietly `command2' `treat2' `mcontrol' [`weight'] if `touse', `fetype'
	quietly predict `treat_res' if `touse', `res'

	quietly sum `treat_res'  if `touse'
	scalar `sigma_xx'=r(Var)

	* Rest includes mcontrol
	* For areg, this regression must absorb fe
  * For reghdfe,
	local rest: list indepvars - treatment
	quietly `command' `treat2' `rest' [`weight'] if `touse', `fetype' `absorb' `residuals'
	quietly predict `error_hat' if `touse', `res'
	quietly sum `error_hat' if `touse'
	quietly {
		scalar `t_x'=r(Var)
	}

	scalar `scbeta'=`beta'
	scalar `scdelta'=`delta'
	scalar `scrmax'=`rmax'

	qui est restore `regmodel'
*------------------------ 2.2: Solve for outputs in mata ----------------------------------

	mata: beta_o=main("`beta_o'","`beta_tilde'","`r_o'","`r_tilde'","`sigma_yy'","`sigma_xx'","`t_x'","`scdelta'","`scbeta'","`scrmax'")
	loc boundx=boundx

	if `delta'==1 {
		loc j=2
		foreach x in betax altroot1 distx dist1 {
			loc `x'=`x'
		}
	}

	else if `delta'~=1 {
		foreach x in betax altroot1 altroot2 distx dist1 dist2 {
			loc `x'=`x'
		}
		if `altroot1'==. & `altroot2'==. loc j=1
		else if (`altroot1'==. & `altroot2'!=. | `altroot2'!=. & `altroot2'==.) loc j=2
		else loc j=3
	}

	/*=========================================================================
					3: Output
	===========================================================================*/

		// Case 2: beta
	if "`type'"=="beta" {

		di _n as txt ///
		_col(18) "{hline 4} Treatment Effect Estimate {hline 4}" _n ///
		_col(14) "{c |}" _col(20) "Estimate" _col(39) "Sq. difference" _col(59) "Bias changes" _n ///
		_col(14) "{c |}" _col(37) "from controlled beta" _col(61) "direction" _n ///
		_col(1) "{hline 13}{c +}{hline 64}" _n ///
		_col(1) "Beta" _col(14) "{c |}" _col(18) as result %11.5f `betax' _col(39) %11.3g `distx' _col(63) "`markx'"  _n ///
		_col(1) as txt "Alt. sol. 1" _col(14) "{c |}" _col(18) as result %11.5f `altroot1' _col(39) %11.3g `dist1' _col(63) "`mark1'"  _n ///
		_col(1) as txt "Alt. sol. 2" _col(14) "{c |}" _col(18) as result %11.5f `altroot2' _col(39) %11.3g `dist2' _col(63) "`mark2'"   _n ///
		as txt _col(1) "{hline 13}{c +}{hline 64}"

		di _n as txt ///
		_col(18) "{hline 4} Inputs from Regressions {hline 4}" _n ///
		_col(14) "{c |}" _col(21) "Coeff." _col(39) "R-Squared" _n ///
		_col(1) "{hline 13}{c +}{hline 64}" _n ///
		_col(1) "Uncontrolled" _col(14) "{c |}" _col(18) as res %12.5f `beta_o' _col(39) %5.3f `r_o' _n ///
		_col(1) as txt "Controlled" _col(14) "{c |}" _col(18) as res %12.5f `beta_tilde' _col(39) %5.3f `r_tilde' _n ///
		_col(1) as txt "{hline 13}{c +}{hline 64}"

		di _n as txt ///
		_col(18) "{hline 4} Other Inputs {hline 4}" _n ///
		_col(1) "{hline 13}{c +}{hline 64}" _n ///
		_col(1) "R_max" _col(14) "{c |}" _col(18) %5.3f `rmax' _n ///
		_col(1) "Delta" _col(14) "{c |}" _col(18) %5.3f `delta' _n ///
		_col(1) "Unr. Controls" _col(14) "{c |}" _col(18) "`mcontrol'"  _n ///
		_col(1) "{hline 13}{c +}{hline 64}"

		return scalar beta=betax
		return scalar dist=distx
		if `j'==2 | `j'==3 {
			return scalar altsol1=altroot1
			return scalar altdist1=dist1
		}
		if `j'==3 {
			return scalar altsol2=altroot2
			return scalar altdist2=dist2
		}
		return scalar root_count=`j'
		return scalar delta=`delta'

	}

		// Case 3: delta
	if "`type'"=="delta" {

		di _n as txt ///
		_col(18) "{hline 4} Bound Estimate {hline 4}" _n ///
		_col(1) "{hline 13}{c +}{hline 64}" _n ///
		_col(1) "delta" _col(14) "{c |}" _col(18) as result %11.5f `boundx' _n ///
		as txt _col(1) "{hline 13}{c +}{hline 64}"

		di _n as txt ///
		_col(18) "{hline 4} Inputs from Regressions {hline 4}" _n ///
		_col(14) "{c |}" _col(21) "Coeff." _col(49) "R-Squared" _n ///
		_col(1) "{hline 13}{c +}{hline 64}" _n ///
		_col(1) "Uncontrolled" _col(14) "{c |}" _col(18) as res %12.5f `beta_o' _col(49) %5.3f `r_o' _n ///
		_col(1) as txt "Controlled" _col(14) "{c |}" _col(18) as res %12.5f `beta_tilde' _col(49) %5.3f `r_tilde' _n ///
		_col(1) as txt "{hline 13}{c +}{hline 64}"

		di _n as txt ///
		_col(18) "{hline 4} Other Inputs {hline 4}" _n ///
		_col(1) "{hline 13}{c +}{hline 64}" _n ///
		_col(1) "R_max" _col(14) "{c |}" _col(18) %5.3f `rmax' _n ///
		_col(1) "Beta" _col(14) "{c |}" _col(18) %9.6f `beta' _n ///
		_col(1) "Unr. Controls" _col(14) "{c |}" _col(18) "`mcontrol'"  _n ///
		_col(1) "{hline 13}{c +}{hline 64}"

		return scalar delta=boundx
		return scalar beta=`beta'
	}

	return scalar rmax=`rmax'
	return local cmd=e(cmd)
	return local depvar="`depvar'"
	return local treatment="`treatment'"
	return local indepvars="`indepvars'"
	return local mcontrol="`mcontrol'"
	return local type="`type'"




end

mata

	// Solve for Delta for Beta = beta-hat

	function bound( real scalar ds,
					real scalar bo_m_bt,
					real scalar rt_m_ro_t_syy,
					real scalar rm_m_rt_t_syy,
					real scalar beta_tilde,
					real scalar beta,
					real scalar t_x,
					real scalar sigma_xx)
	{
		real scalar bt_m_b, num1, num2, num3, num4, den1, den2, den3, den4, num, den

		bt_m_b = beta_tilde - beta


		num1 = (bt_m_b)*rt_m_ro_t_syy*t_x
		num2 = (bt_m_b)*sigma_xx*t_x*(bo_m_bt)^2
		num3 = 2*(bt_m_b)^2*(t_x*bo_m_bt*sigma_xx)
		num4 = ((bt_m_b)^3)*((t_x*sigma_xx-t_x^2))

		den1 = rm_m_rt_t_syy*bo_m_bt*sigma_xx
		den2 = bt_m_b*rm_m_rt_t_syy*(sigma_xx-t_x)
		den3 = (bt_m_b^2)*(t_x*bo_m_bt*sigma_xx)
		den4 = (bt_m_b^3)*(t_x*sigma_xx-t_x^2)

		num=num1+num2+num3+num4
		den=den1+den2+den3+den4

		ds=num/den
	}

	/// Solutions for beta if delta==1
	/// Quadratic function
	/// v = (- cap_theta +/- sqrt ( cap_theta + d1_1))/(d1_2)

	function d1quadsol( real scalar rm_m_rt_t_syy,
						real scalar rt_m_ro_t_syy,
						real scalar bo_m_bt,
						real scalar sigma_xx,
						real scalar t_x,
						real scalar beta_o,
						real scalar betax,
						real scalar altsol1,
						real scalar beta_tilde,
						real scalar distx,
						real scalar dist1,
						real scalar markx,
						real scalar mark1
	)
	{
		real scalar cap_theta, d1_1, d1_2, sol1, sol2, beta1, beta2, solc

		cap_theta = rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*(bo_m_bt^2)
		d1_1 = 4*rm_m_rt_t_syy*(bo_m_bt^2)*(sigma_xx^2)*t_x
		d1_2 = -2*t_x*bo_m_bt*sigma_xx

		sol1 = (-1*cap_theta-sqrt((cap_theta)^2+d1_1))/(d1_2)
		sol2 = (-1*cap_theta+sqrt((cap_theta)^2+d1_1))/(d1_2)

		beta1 = beta_tilde - sol1
		beta2 = beta_tilde - sol2

		if ( (beta1-beta_tilde)^2 < (beta2-beta_tilde)^2) {
			betax = beta1
			altsol1 = beta2
		}
		else {
			betax = beta2
			altsol1 = beta1
		}
		/// Change to altsol if bias of beta_o vs beta_tilde is of different sign that bias of beta_tilde - beta_x
		if ( sign(betax-beta_tilde)!=sign(beta_tilde-beta_o) ) {
			solc = betax
			betax = altsol1
			altsol1 = solc
		}

		/// Mark solutions that violate assumption 3

		if ( sign(betax-beta_tilde)!=sign(beta_tilde-beta_o) ) markx=1
		if ( sign(altsol1-beta_tilde)!=sign(beta_tilde-beta_o) ) mark1=1


		distx = (betax - beta_tilde)^2
		dist1 = (altsol1 - beta_tilde)^2

	}

	/// Solve for cubic roots if delta is not 1

	function dnot1cubsol(	real scalar bo_m_bt,
						real scalar sigma_xx,
						real scalar delta,
						real scalar t_x,
						real scalar rm_m_rt_t_syy,
						real scalar rt_m_ro_t_syy,
						real scalar betax,
						real scalar altsol1,
						real scalar altsol2,
						real scalar distx,
						real scalar dist1,
						real scalar dist2,
						real scalar beta_tilde,
						real scalar beta_o,
						real scalar markx,
						real scalar mark1,
						real scalar mark2

	)
	{
		real scalar A, B, C, Q, R, D, discrim, theta, sol1, sol2, sol3, w, t1, t2, crt1, crt2
		real matrix sols, dists, min
		/// ax3 + bx2 + cx +d = 0
		///Recasted as
		/// x3 + b/a x2 + c/a x + d/a = 0
		/// Call this
		/// x3 + A x2 + B x + C = 0


		A= (t_x*bo_m_bt*sigma_xx*(delta-2))/((delta-1)*(t_x*sigma_xx-t_x^2))
		B= (delta*rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*bo_m_bt'^2)/((delta-1)*(t_x*sigma_xx-t_x^2))
		C = (rm_m_rt_t_syy*delta*bo_m_bt*sigma_xx)/((delta-1)*(t_x*sigma_xx-t_x^2))



		Q = (A^2-3*B)/9
		R = (2*A^3-9*A*B+27*C)/54
		D = R^2-Q^3
		discrim = R^2-Q^3



		if (discrim <0) {

			theta = acos(R/sqrt(Q^3))


			sol1 = -2*sqrt(Q)*cos(theta/3)-(A/3)
			sol2 = -2*sqrt(Q)*cos((theta+2*pi())/3)-(A/3)
			sol3 = -2*sqrt(Q)*cos((theta-2*pi())/3)-(A/3)



			sols= (sol1,sol2,sol3)
			sols =  beta_tilde :-sols

			dists = sols :- beta_tilde
			dists = dists:^2


			/// Change to alternative solutions if first solution violates assumption 3
			for (i=1; i<=3; i++) {
				if ( sign(sols[i]-beta_tilde)!=sign(beta_tilde-beta_o) ) {
					dists[i]=max(dists)+1
				}
			}

			min=J(1,3,.)
			w=J(1,3,.)

			minindex(dists,3,min,w)



			betax = sols[min[1]]
			altsol1 = sols[min[2]]
			altsol2 = sols[min[3]]

			distx= dists[min[1]]
			dist1= dists[min[2]]
			dist2= dists[min[3]]

			/// Mark solutions that violate assumption 3

			if ( sign(betax-beta_tilde)!=sign(beta_tilde-beta_o) ) markx=1
			if ( sign(altsol1-beta_tilde)!=sign(beta_tilde-beta_o) ) mark1=1
			if ( sign(altsol2-beta_tilde)!=sign(beta_tilde-beta_o) ) mark2=1

		}



		else {


			t1=-1*R+sqrt(D)
			t2=-1*R-sqrt(D)

			/// Watch the cuberoots

			crt1=sign(t1) * abs(t1)^(1/3)
			crt2=sign(t2) * abs(t2)^(1/3)

			sol1 = crt1+crt2-(A/3)
			betax=beta_tilde-sol1

		}

	}


	void main(  string scalar Beta_o,
				string scalar Beta_tilde,
				string scalar R_o,
				string scalar R_tilde,
				string scalar Sigma_yy,
				string scalar Sigma_xx,
				string scalar T_x,
				string scalar Delta,
				string scalar Beta,
				string scalar R_max
	)
	{
		real scalar beta_o, beta_tilde, r_o, r_tilde, sigma_yy, sigma_xx, t_x, delta, beta, r_max

	beta_o=st_numscalar(Beta_o)
	beta_tilde=st_numscalar(Beta_tilde)
	r_o=st_numscalar(R_o)
	r_tilde=st_numscalar(R_tilde)
	sigma_yy=st_numscalar(Sigma_yy)
	sigma_xx=st_numscalar(Sigma_xx)
	t_x=st_numscalar(T_x)
	delta=st_numscalar(Delta)
	beta=st_numscalar(Beta)
	r_max=st_numscalar(R_max)

	bo_m_bt = beta_o - beta_tilde
	rt_m_ro_t_syy = (r_tilde - r_o)*sigma_yy
	rm_m_rt_t_syy = (r_max - r_tilde)*sigma_yy


	/// Solve delta bound
	bound(ds=.,bo_m_bt,rt_m_ro_t_syy,rm_m_rt_t_syy,beta_tilde,beta,t_x,sigma_xx)
	st_numscalar("boundx",ds)

	/// Solve quadratic if delta = 1
	if (delta==1)  d1quadsol(rm_m_rt_t_syy,rt_m_ro_t_syy,bo_m_bt,sigma_xx,t_x,beta_o,betax=.,altsol1=.,beta_tilde,distx=., dist1=.,markx=.,mark1=.)
	/// Solve cubic if delta !=1
	else dnot1cubsol(bo_m_bt,sigma_xx,delta,t_x,rm_m_rt_t_syy,rt_m_ro_t_syy,betax=.,altsol1=.,altsol2=.,distx=.,dist1=.,dist2=.,beta_tilde,beta_o,markx=.,mark1=.,mark2=.)

	st_numscalar("betax",betax)
	st_numscalar("altroot1",altsol1)
	st_numscalar("altroot2",altsol2)
	if (markx==1) st_local("markx","Yes")
	if (mark1==1) st_local("mark1","Yes")
	if (mark2==1) st_local("mark2","Yes")

	st_numscalar("distx",distx)
	st_numscalar("dist1",dist1)
	st_numscalar("dist2",dist2)

	}




end
