# `psacalc` with support for `reghdfe`
A minor revision to `psacalc` to add support for `reghdfe`.

`psacalc` by Emily Oster is a Stata command (available from [SSC](https://ideas.repec.org/s/boc/bocode.html): `ssc install psacalc` which implements the method outlined in Oster (2017) to calculate treatment effects and relative degree of selection under proportional selection of observables and unobservables (see `help psacalc` after installation for details).

The current version of `psacalc` (22 Oct 2020) supports linear models estimated using `regress`, `areg`, and `xtreg`; however, the popular [`reghdfe`](http://scorreia.com/software/reghdfe/) is not supported. [See this Statalist discussion.](https://www.statalist.org/forums/forum/general-stata-discussion/general/1523251-combining-ivreghdfe-and-psacalc)

In this repo I propose a simple edit to the `psacalc.ado` file available from SSC which adds support for `reghdfe`. You can easily see what I have changed from the version on SSC by inspecting the commit history of this repo.

When using `areg` or `reghdfe`, my edits also allow the user to treat _all_ covariates as nuisance parameters thus calculating treatment effects and relative degree of selection based only on the fixed effects. For example:

```Stata
clear all

// make some data
set obs 2000
set seed 101
generate f = rnormal() // the individual fixed effect
generate id = _n
expand 10 // create 10 years of data
bysort id: generate year=_n
generate x1 = rnormal()+2*f
generate x2 = 2*rnormal()+f
generate u = 3*rnormal()+5*f
generate y = x1 + x2 + u

// psacalc2 is allowed after reghdfe
reghdfe y x1 x2, absorb(id)
psacalc2 delta x1

// psacalc2 is also allowed if the only controls are those in `absorb()`
reghdfe y x1, absorb(id)
psacalc2 delta x1

// psacalc2 treating all parameters as nuisance but those in `absorb()`
reghdfe y x1 x2, absorb(id)
psacalc2 delta x1, mcontrol(x2)
```

All my edits are based on the assumption that the <img src="https://render.githubusercontent.com/render/math?math=R^2">s calculated by `areg` and `reghdfe` are comparable researchers should always verify that any estimation software is doing what they think it is before proceeding. Which is to say, any errors in my edits are my own, but you own any errors of your own.

I welcome any comments or improvements.

### Future work:

- Allow `psacalc` to treat the fixed-effects absorbed by `reghdfe` as nuisance parameters. `psacalc` currently allows this for any control, but only for fixed effects estimated with `xtreg`.

### References:
Oster, Emily (2017). "Unobservable Selection and Coefficient Stability: Theory and Validation", _Journal of Business Economics and Statistics_,  [DOI: 10.1080/07350015.2016.1227711]( https://doi.org/10.1080/07350015.2016.1227711)
