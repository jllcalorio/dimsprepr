## Submission

This is a new release.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'John Lennon Calorio <jllcalorio@gmail.com>'
  New submission

  Possibly misspelled words in DESCRIPTION: Broadhurst, Dieterle, Kirwan,
  Vynck, et al, QComics, Metabolomics, metabolomics, spectrometry

  These are correctly spelled: author surnames of the cited methods, the
  name of the QComics guidelines, and domain terms.

## Test environments

* win-builder: R 4.6.2 (release) and R-devel (2026-11-17 r90926), both
  1 NOTE (above).
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
