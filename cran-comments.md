## Submission

This is a new submission.

sorenR is an experimental, pure-R implementation. It depends only on base and
recommended packages at run time and needs no network access, external
software or compiled code. All examples, tests and the vignette use small
synthetic data and run in a few seconds.

## Test environments

* local: macOS (aarch64), R 4.6.1, `R CMD check --as-cran --no-manual`
  (no LaTeX installation is available locally, so the PDF manual was not
  built there)
* GitHub Actions: macos-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1), with `--as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission
