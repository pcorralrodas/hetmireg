*! version 0.2 March27-2018
*! Added first, and made GLS function faster, and more memory friendly
*! Paul Corral - pcorralrodas@worldbank.org 
*! World Bank Group - Poverty and Equity Global Practice 
*! Global Solutions Group on Welfare Measurement and Statistical Capacity


cap program drop hetmireg
program define hetmireg, eclass byable(recall)
	version 11, missing
	#delimit ;
	syntax varlist(min=2 numeric fv) [if] [in] [aw], 
	SIMs(integer)
	uniqid(varname)
	errdraw(string)
	by(varname)
	[Robust 
	yhat(varlist numeric min=1 fv)
	yhat2(varlist numeric min=1 fv)
	lny 
	seed(numlist)
	SIMName(string)
	het(varlist numeric min=1 fv)
	mlong
	first
	];
#delimit cr	
set more off

	tokenize `varlist'
	
// Local for dependent variable
local depvar `1'

// obtain the independent variables
macro shift 
local indeps `*'
local Ind `*'

//Remove collinear exolanatory vars
_rmcoll `indeps', forcedrop
local indeps  `r(varlist)'

tempname beta sigma sigma2 vcov df var_r alfa a_vcov dfa thetarg theexp a 
tempvar error1 touse touse2 lnyhat usehet bigb sige2 nn

//INDICATE HETEROSKEDASTICITY MODELS
if ("`het'"!="" | "`yhat'"!="" | "`yhat2'"!=""){
	local hetero = 1
	if ("`het'"!="" & "`yhat'"=="" & "`yhat2'"=="") local zcond = "100"
	if ("`het'"!="" & "`yhat'"!="" & "`yhat2'"=="") local zcond = "110"
	if ("`het'"!="" & "`yhat'"!="" & "`yhat2'"!="") local zcond = "111"
	if ("`het'"=="" & "`yhat'"!="" & "`yhat2'"=="") local zcond = "010"
	if ("`het'"=="" & "`yhat'"=="" & "`yhat2'"!="") local zcond = "001"
	if ("`het'"=="" & "`yhat'"!="" & "`yhat2'"!="") local zcond = "011"
	if ("`het'"!="" & "`yhat'"=="" & "`yhat2'"!="") local zcond = "101"
}
else local hetero = 0

if (`hetero'==1){
	qui:gen `usehet'=1
	foreach x of varlist `depvar' `indeps' `het' `yhat' `yhat2'{
		qui:replace `usehet'=0 if missing(`x')
	}
		
}

	//Weights
	local wvar : word 2 of `exp'
	if "`wvar'"=="" {
		tempvar w
		gen `w' = 1
		local wvar `w'
	}

local errdraw=trim(lower("`errdraw'"))

if (("`by'"=="" & "`target'"=="") ) {
display "OK"
	display as error "You must specify a variable to differentiate the sources, or the target data path"
	error 1000000000
	exit
}

if (("`by'"!="" & "`target'"!="")) {
display "OK"
	display as error "You must specify a variable to differentiate the sources, or the target data path, not both"
	error 1000000000
	exit
}

if ("`errdraw'"!="normal" & "`errdraw'"!="empirical"){
	display as error "You must specify either normal, or empirical drawing of residuals"
	error 1000000000
	exit
}

if ("`seed'"=="") set seed `c(rngstate)'
else set seed `seed'

if ("`mlong"==""){
	sort `by' `uniqid'
	
	if ("`by'"!=""){
		qui:tab `by' 
		if r(r)!=2{
		display as error "You must specify two groups, source and target data"
		error 498
		exit
		}
	}
}



if (`hetero'!=1){ 
	reg `depvar' `indeps' [aw=`wvar'], `robust'
	qui: gen `touse'=e(sample)
	qui: predict `error1' if `touse', res
	mata: st_view(dres=.,.,"`error1'", "`touse'")
	mat `beta' = e(b)
	mat `vcov' = e(V)
	mat `sigma'= e(rmse)
	//degrees of freedom
	scalar `df' = e(df_r)
	mata: beta   =st_matrix("`beta'")
	mata: vcov   =st_matrix("`vcov'")
	mata: sigma2 =(st_matrix("`sigma'"))^2
	mata: df     =st_numscalar("`df'")
}
else{
	display in yellow "First Stage of Heteroscedastic Model"
	reg `depvar' `indeps' if `usehet'==1 [aw=`wvar']
	qui: gen `touse'=e(sample)
	qui: predict `lnyhat' if `touse'==1, xb
	qui: predict `error1' if `touse'==1, res
	scalar `df' = e(df_r)
		mata: st_view(lnyhat=.,.,"`lnyhat'", "`touse'")	
		mata: st_view(ech=.,.,"`error1'", "`touse'")
		mata: maxA      = 1.05*max(ech:^2)
		mata: st_local("maxA",strofreal(maxA))
		mata: ech[.,.]    = ln((ech:^2):/(maxA:-(ech:^2)))
		mata: st_view(_x  = .,.,"`indeps'", "`usehet'")		
		mata: st_view(_y  = .,.,"`depvar'", "`usehet'")
	if ("`yhat'"!=""){
		foreach x of local yhat{
		tempname `x'_yhat
		qui: gen ``x'_yhat' = `lnyhat'*`x' if `touse'==1
		qui: replace ``x'_yhat' = `x' if `touse'==0
		local yhat_e `yhat_e' ``x'_yhat'
		}
	}
	if ("`yhat2'"!=""){
		foreach x of local yhat2{
		tempname `x'_yhat2
		qui: gen ``x'_yhat2' = `lnyhat'*`x'*`lnyhat' if `touse'==1
		qui: replace ``x'_yhat2' = `x' if `touse'==0
		local yhat2_e `yhat2_e' ``x'_yhat2'
		}
	}
	display in yellow "Alpha Model for Heteroscedasticity"
	reg `error1' `het' `yhat_e' `yhat2_e' if `usehet'==1 [aw=`wvar']
	scalar `dfa'=e(df_r)
	mata: dfa  =st_numscalar("`dfa'")
	qui: predict double `bigb' if `usehet'==1, xb
	qui: mat `alfa'  =e(b)
	qui: mat `a_vcov'=e(V)
	qui: mata: alfa   =st_matrix("`alfa'")
	qui: mata: a_vcov =st_matrix("`a_vcov'")
	qui: replace `bigb'=exp(`bigb')
	qui: scalar `var_r' = e(rmse)^2
	mata: s2 = st_numscalar("`var_r'")
	qui:gen `sige2'=(`maxA'*(`bigb'/(1+`bigb'))) + (0.5*`var_r')*((`maxA'*(`bigb'*(1-`bigb')))/((1+`bigb')^3))  
	mata: st_view(ech=.,.,"`sige2'","`usehet'")
	mata: st_view(wt=.,.,"`wvar'", "`usehet'")
	mata: sig_n=0
	mata: info=1,rows(ech)
}
if ("`target'"!=""){
	use "`target'", clear
	qui:gen `touse2' = 1
}
qui: levelsof `by' if `depvar'!=., local(`thetarg')
else                gen `touse2' = `touse'==0 & `by'!=``thetarg''

//Simulation on to the target data
//...indicate observations to use

foreach x in `indeps' `het' `yhat' `yhat2'{
	qui:replace `touse2'=0 if `touse2'==1 & `x'==.
}

if ("`simname'"=="") local simname yhat_
if ("`mlong'"==""){
	forval z=1/`sims'{
		qui: gen double `simname'`z'=.
		local thev `thev' `simname'`z'
	}
}
else{	
	local thev `depvar'
}




if (`hetero'==0){
	if ("`first'"==""){
		if "`mlong'"=="mlong"{
			if ("`lny'"=="lny") mata: _1y = vec(exp(_mi_sim("`indeps'","`touse2'",beta,vcov,df,sigma2,dres,`sims')))
			else                mata: _1y = vec((_mi_sim("`indeps'","`touse2'",beta,vcov,df,sigma2,dres,`sims')))		
			//local sims = `sims'+1
			tempfile _nnio
			gen `theexp'=`by'!=``thetarg''
			qui:count 
			local hu0 = r(N)
			preserve
				qui:keep if `theexp'==1
				tempfile _nnio
				qui:save `_nnio'
				qui:count
				local hu=r(N)
			restore
				local `a'=1
				while (``a''<=`sims'){
					qui:append using `_nnio', gen(`nn')
					qui:replace _mi_m=``a'' if `nn'==1
					local `a' = ``a''+1
					drop `nn'
				}
			qui:count
			qui: char _dta[_mi_M] `sims'
			qui: char _dta[_mi_n] `hu'
			qui: char _dta[_mi_N] `hu0'
			replace `touse2'=`touse2'==1 & _mi_m>0
			mata: st_store(.,st_varindex(tokens("`thev'")),"`touse2'",_1y)
		}
		else{
			if ("`lny'"=="lny") mata: st_store(.,st_varindex(tokens("`thev'")),"`touse2'",exp(_mi_sim("`indeps'","`touse2'",beta,vcov,df,sigma2,dres,`sims')))
			else                mata: st_store(.,st_varindex(tokens("`thev'")),"`touse2'",(_mi_sim("`indeps'","`touse2'",beta,vcov,df,sigma2,dres,`sims')))
		}
	}
}
else{
	mata: b_gls=_f_hh_gls4(_y,_x,wt, ech, sig_n, info,1 ,0)
		//HT insert here matrix from outside... use mkmat
	mata: st_matrix("`beta'",*b_gls[1,1])
	mata: st_matrix("`vcov'",(*b_gls[1,2]))
	//Need to do drawings of betas and other parameters for heterosk imputations
	qui mat rownames `beta' = `indeps' _cons
	qui mat `beta' = `beta''
	qui mat rownames `vcov' = `indeps' _cons
	qui mat colnames `vcov' = `indeps' _cons
	ereturn post `beta' `vcov', depname(`depvar') esample(`touse') 
	display in yellow "GLS coefficients"
	ereturn display
	if ("`first'"==""){
		if "`mlong'"=="mlong"{
			if ("`lny'"=="lny") mata: _1y = vec(exp(_mi_simhet("`indeps'","`touse2'", *b_gls[1,1]',*b_gls[1,2],s2,*b_gls[1,3],`sims',"`het' `yhat_e' `yhat2_e'",alfa, a_vcov, maxA, dfa)))
			else                mata: _1y = vec(_mi_simhet("`indeps'","`touse2'", *b_gls[1,1]',*b_gls[1,2],s2,*b_gls[1,3],`sims',"`het' `yhat_e' `yhat2_e'",alfa, a_vcov, maxA, dfa))		
			//local sims = `sims'+1
			tempfile _nnio
			gen `theexp'=`by'!=``thetarg''
			qui:count 
			local hu0 = r(N)
			preserve
				qui:keep if `theexp'==1
				tempfile _nnio
				qui:save `_nnio'
				qui:count
				local hu=r(N)
			restore
				local `a'=1
				while (``a''<=`sims'){
					qui:append using `_nnio', gen(`nn')
					qui:replace _mi_m=``a'' if `nn'==1
					local `a' = ``a''+1
					drop `nn'
				}
			count
			qui: char _dta[_mi_M] `sims'
			qui: char _dta[_mi_n] `hu'
			qui: char _dta[_mi_N] `hu0'
			replace `touse2'=`touse2'==1 & _mi_m>0
			mata: st_store(.,st_varindex(tokens("`thev'")),"`touse2'",_1y)
		}
		else{
			if ("`lny'"=="lny") mata: st_store(.,st_varindex(tokens("`thev'")),"`touse2'",exp(_mi_simhet("`indeps'","`touse2'", *b_gls[1,1]',*b_gls[1,2],s2,*b_gls[1,3],`sims',"`het' `yhat_e' `yhat2_e'",alfa, a_vcov, maxA, dfa)))
			else                mata: st_store(.,st_varindex(tokens("`thev'")),"`touse2'",(_mi_simhet("`indeps'","`touse2'", *b_gls[1,1]',*b_gls[1,2],s2,*b_gls[1,3],`sims',"`het' `yhat_e' `yhat2_e'",alfa, a_vcov, maxA, dfa)))
		}
	}
}
ereturn local _thev  `"`thev'"'
ereturn local df = `df'
if (`hetero'==1){
	ereturn matrix a_b  = `alfa'
	ereturn matrix a_V  =`a_vcov'
	ereturn local maxA = `maxA'
	ereturn local var_r = `var_r'
	ereturn local dfa = `dfa'
}
end

mata
//function to draw multivariate normal distribution
function _f_dnorm(real scalar n, real matrix M, real matrix V) {
	return(M :+ invnormal(uniform(n,cols(V)))*cholesky(V)')
}
function _f_sampleepsi(real scalar n, real scalar dim, real matrix eps) {				  
	sige2 = J(dim,n,0)
	N = rows(eps)
	if (cols(eps)==1) for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),1]
	else              for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),i]
	//for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(rows(eps)*runiform(dim,1)),i]
	return(sige2)	
}

function _mi_sim(string scalar xname, 
				string scalar touse, 
				real matrix beta, 
				real matrix vcov, 
				real scalar df, 
				real scalar sigma2, 
				real matrix epsi,
				real scalar nsim)
{
	st_view(x=.,.,tokens(xname),touse)
	if (st_local("errdraw")=="normal"){
		 sigma2 = (sigma2*df):/rchi2(1,nsim,df)
		 y = rnormal(1,1,quadcross((x, J(rows(x),1,1))',_f_dnorm(nsim,beta,vcov)'),sqrt(sigma2))
	}
	else y = quadcross((x, J(rows(x),1,1))',_f_dnorm(nsim,beta,vcov)') + _f_sampleepsi(nsim, rows(x), epsi)
	return(y)
}

//Function for heteroscedastic draws
function _mi_simhet(string scalar xname, 
					string scalar touse, 
					real matrix beta, 
					real matrix vcov, 
					real scalar sigma2, 
					real matrix epsi,
					real scalar nsim,
					string scalar zname,
					real matrix alfa,
					real matrix a_vcov,
					real scalar maxA,
					real scalar dfa)
{
	st_view(x=.,.,tokens(xname),touse)
	ztype = st_local("zcond")	
	if (st_local("errdraw")=="normal"){
		y = quadcross((x, J(rows(x),1,1))',_f_dnorm(nsim,beta,vcov)')		
		y = y + rnormal(1,1,0,sqrt(_sig2_het(sigma2,maxA,alfa,a_vcov,dfa,nsim,ztype, touse, zname, y)))		
	}
	else{
		y = quadcross((x, J(rows(x),1,1))',_f_dnorm(nsim,beta,vcov)') + _f_sampleepsi(nsim, rows(x), epsi)
	}
	return(y)
}

function _sig2_het(s2, maxA, alfa, a_vcov, dfa, nsim,ztype, touse, zname, y){
	alfa = _f_dnorm(nsim,alfa,a_vcov)'
	//Bring in the Zs
	if (st_local("het")!="") st_view(z=.,.,tokens(st_local("het")),touse)
	if (st_local("yhat_e")!="") st_view(z1=.,.,tokens(st_local("yhat_e")),touse)
	if (st_local("yhat2_e")!="") st_view(z2=.,.,tokens(st_local("yhat2_e")),touse)
	
	if (ztype=="100") bigb = (z,J(rows(z),1,1))*alfa
	if (ztype=="110") bigb = z*alfa[|1,.\cols(z),.|] + z1*alfa[|(cols(z)+1),. \ (cols(z1)+cols(z)),.|]:*y + J(rows(y),1,1)*alfa[rows(alfa),.]
	if (ztype=="111") bigb = z*alfa[|1,.\cols(z),.|] + z1*alfa[|(cols(z)+1),. \ (cols(z1)+cols(z)),.|]:*y + z2*alfa[|(cols(z)+cols(z1)+1),. \ (cols(z)+cols(z1)+cols(z2)),.|]:*(y:^2) + J(rows(y),1,1)*alfa[rows(alfa),.]
	if (ztype=="101") bigb = z*alfa[|1,.\cols(z),.|] + z2*alfa[|(cols(z)+1),. \ (cols(z2)+cols(z)),.|]:*(y:^2) + J(rows(y),1,1)*alfa[rows(alfa),.]
	if (ztype=="010") bigb = z1*alfa[|1,. \ cols(z1),.|]:*y + J(rows(y),1,1)*alfa[rows(alfa),.]
	if (ztype=="011") bigb = z1*alfa[|1,. \ cols(z1),.|]:*y + z2*alfa[|(cols(z1)+1),. \ (cols(z1)+cols(z2)),.|]:*(y:^2) + J(rows(y),1,1)*alfa[rows(alfa),.]
	if (ztype=="001") bigb = z2*alfa[|1,. \ cols(z2),.|]:*(y:^2) + J(rows(y),1,1)*alfa[rows(alfa),.]
	s2 = (s2*dfa):/rchi2(1,nsim,dfa)	
	bigb=exp(bigb)
	bigb=(maxA:*(bigb:/(1:+bigb))):+ (0.5:*s2):*((maxA:*(bigb:*(1:-bigb))):/((1:+bigb):^3))
	return(bigb)
}



function _f_hh_gls4(real matrix y, 
					real matrix x, 
					real matrix wt, 
					real matrix sig_e, 
					real scalar sig_n, 
					real matrix info, 
					real scalar EB, 
					real scalar bs) 
{
	pointer(real matrix) rowvector glsout3
	glsout3 = J(1,3,NULL)
	x=x,J(rows(x),1,1)
	
	//This can still be sped up!!
		cv= wt:/sig_e			
		
	//Capital sigma matrix for GLS	
	xy=quadcross(x,cv,y)
	ixx=invsym(quadcross(x,cv,x))
	b=quadcross(ixx,xy)

		xtwewx = quadcross(x,(cv:*wt),x)
		
	//VCov
	vcov2 = quadcross(quadcross(ixx,xtwewx)',ixx)
	
	glsout3[1,1] = &(b)
	glsout3[1,2] = &(vcov2)
	glsout3[1,3] = &(y -quadcross(x',b))
	return(glsout3)

}

end
