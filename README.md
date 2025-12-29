# hetmireg — Multiple imputation with empirical and heteroscedastic errors

## Overview

`hetmireg` implements a multiple-imputation workflow for welfare (or other continuous outcomes) that supports:

* **Empirical error draws** (non-parametric re-sampling from residuals) as an alternative to normality.
* **Heteroscedastic error modeling** via a user-specified variance (or scale) specification.
* A **source/target** setup (via `by()`) to support cross-sample imputation contexts.

This command is designed for simulation-based imputation where the analyst wants to relax homoscedastic-normal error assumptions and better reflect observed residual behavior.

---

## Installation

* Place the package files (`hetmireg.ado`, help file, etc.) on your Stata ado-path, or install from your chosen repository.

---

## Syntax

```stata
hetmireg depvar [indepvars] [if] [in] [aw] , sims(integer) ///
    uniqid(varname) errdraw(string) by(varname) ///
    [ robust het(varlist) yhat(varlist) yhat2(varlist) lny ///
      seed(numlist) simname(string) first ]
```

---

## Required options

* `sims(#)`
  Number of simulations / imputations to generate.

* `uniqid(varname)`
  Unique identifier used to ensure replicability (must uniquely identify observations).

* `errdraw(string)`
  Error draw method. Supported values:

  * `normal` — draws errors from a normal distribution
  * `empirical` — draws errors empirically (re-sampling residuals)

* `by(varname)`
  Identifies source vs. target dataset (or domains) used in the imputation setup.

---

## Optional options

* `robust`
  Uses robust standard errors in the OLS estimation stage.

* `het(varlist)`
  Variables used for the heteroscedastic specification.

* `yhat(varlist)`
  Variables to interact with the fitted value (`yhat`) in the heteroscedastic model.

* `yhat2(varlist)`
  Variables to interact with the fitted value squared (`yhat^2`) in the heteroscedastic model.

* `lny`
  Indicates the dependent variable is in logs (log-scale model / log-dependent variable handling).

* `seed(numlist)`
  Sets the seed for replicability. If omitted, `hetmireg` uses Stata’s current RNG state.

* `simname(string)`
  Prefix for the simulated vectors created by the command.

* `first`
  Estimation-only mode: displays the fitted model but does not generate imputations.

---

## Output

Depending on options, `hetmireg` will:

* estimate the baseline model (OLS with optional robust SEs),
* optionally estimate heteroscedastic components,
* generate `sims(#)` simulated/imputed vectors using the selected error draw approach,
* name simulated outputs using `simname()` (if provided).

---

## Example

```stata
* Example skeleton (adjust variable names and options to your application)
hetmireg lny_welfare age educ hhsize, sims(100) ///
    uniqid(hhid) errdraw(empirical) by(sample) ///
    het(region urban) yhat(urban) yhat2(urban) ///
    seed(12345) simname(mi_)
```

---

## Notes and good practice

* Ensure `uniqid()` is truly unique and stable across runs; replicability depends on it.
* When using `errdraw(empirical)`, confirm residuals are defined on the appropriate scale (especially when `lny` is used).
* In cross-sample settings, verify `by()` correctly identifies which observations are used for model estimation versus imputation targets.

---
## Methodological background

hetmireg follows the simulation-based imputation logic of Elbers, Lanjouw, and Lanjouw (2002), extending the core approach by allowing (i) empirical residual draws and (ii) a heteroscedastic error specification within the multiple-imputation routine.

Reference

Elbers, C., Lanjouw, J. O., & Lanjouw, P. (2002). Micro-level estimation of welfare. World Bank Policy Research Working Paper.

## Authors

Paul Corral
The World Bank — Poverty and Equity Global Practice
Washington, DC
[pcorralrodas@worldbank.org](mailto:pcorralrodas@worldbank.org)
