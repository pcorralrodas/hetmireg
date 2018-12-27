{smcl}
{* *! version 1.0.0  10sept2017}{...}
{cmd:help hetmireg}
{hline}

{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:hetmireg} {hline 1}}Multiple imputation with empirical and heteroscedastic errors{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 23 2}
{opt hetmireg} {depvar} [{indepvars}] {ifin} [aw] {cmd:,}
{opt sim:s(integer)}
{opt uniqid(varname)}
{opt errdraw(string)}
{opt by(varname)}
[{opt R:obust}
{opt het({indepvars})}
{opt yhat(string)}
{opt yhat2(string)}
{opt lny}
{opt seed(numlist)}
{opt simn:ame(string)}]

{title:Description}

{pstd}
{cmd:hetmireg} Supports multiple imputation by including empirical drawing of errors, as 
well as heteroscedatic estimation of models.

{title:Options}

{phang}
{opt sims(integer)} Specifies the number of simulations desired

{phang}
{opt uniqid(varname)} Variable used for ensuring replicability of results, it must be a unique id.

{phang}
{opt errdraw(string)} Specifies method for drawing errors in the simulations. Options are: normal, or empirical.

{phang}
{opt by(varname)} Variable used to identify source and target dataset.

{phang}
{opt robust} Specifies that robust standard errors must be obtained in the OLS.

{phang}
{opt het(varlist)} Variables used for heteroscedastic specification.

{phang}
{opt yhat(varlist)} Variables to interact with yhat for heteroscedastic model

{phang}
{opt yhat2(varlist)} Variables to interact with yhat squared for heteroscedastic model

{phang}
{opt lny} Indicates that dependent variable is in log terms

{phang}
{opt seed(numlist)} Seed to ensure replicability, if not specified it uses the current seed's state.

{phang}
{opt simname(string)} Prefix for simulated vectors.

{phang}
{opt first} Specifies that only model be shown, no imputations are done.

{title:Example}

{title:Authors}

{pstd}
Paul Corral{break}
The World Bank - Poverty and Equity Global Practice {break}
Washington, DC{break}
pcorralrodas@worldbank.org{p_end}
