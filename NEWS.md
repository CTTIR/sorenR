# sorenR 0.1.0

* Prepare the package for CRAN: full-sentence description, `URL`,
  `BugReports` and `Language` fields, a spelling word list and
  `cran-comments.md`.
* Add runnable examples to every exported function.
* Record the installed package version in fitted models instead of a fixed
  string.
* Run `R CMD check` on a CRAN-like GitHub Actions matrix (macOS, Windows,
  and Ubuntu with R devel, release and oldrel-1).

# sorenR 0.0.1

* Show the logo in the new dimension reduction and visualization color
  (deep teal #146554), shared with songR; `tools/branding/render.R` recolors
  the approved hexmakR master and keeps its name outlines.
* Standardize the named v5 package logo across README, documentation pages
  and available application headers, with matching light and dark artwork.
* Initial experimental pure-R prototype graph implementation.
* Frozen radial-basis transform and bounded replay/coverage memory.
* Transactional updates with a checked finite-anchor mean displacement budget.
* Fixed, adaptive, frozen and unrestricted modes. Adaptive scores are uncalibrated.
