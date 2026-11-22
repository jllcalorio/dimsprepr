## Submission

This is a new release.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'John Lennon Calorio <jllcalorio@gmail.com>'
  New submission

## Test environments

* local Windows 11 x64, R 4.5.2: `R CMD check --as-cran --run-donttest`,
  with remote URL checks, with all suggested packages installed and again
  with `_R_CHECK_DEPENDS_ONLY_=true` (no suggested packages available).

## Notes for the reviewer

* All examples run in under 5 seconds, including the `\donttest{}` parts.
* Suggested packages (mgcv, sva, pmp, vsn, mice, imputeLCMD) are only needed
  for optional methods and are always used conditionally, in examples,
  tests, and the vignette.
* The methods are described in the references listed in the Description
  field and in each help page.
