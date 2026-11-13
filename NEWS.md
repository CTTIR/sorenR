# sorenR 0.2.0

This version replaces the update mechanism and the layout. Models saved with
0.1.0 cannot be loaded.

* Updates no longer move existing prototypes. New prototypes are appended only
  where observations fall outside conformally calibrated support, and layout
  coordinates stay within a lifetime radius (`stability_radius`) of their
  creation coordinates. `soren_diagnostics()$certificate` reports the bound.
* `soren_transform(reference = "origin")` evaluates the originally fitted map
  from any later model version.
* Support diagnostics now hold a conformal p-value, a novelty flag, the share
  of weight on added prototypes and the displacement from the original map.
* The layout is stress majorization on prototype graph path lengths or
  Euclidean distances, replacing the attraction and repulsion objective.
* The map uses one bandwidth per prototype and optional frozen local linear
  maps (`local_linear`).
* `soren_control()` changed: `stability_budget`, `learning_rate`, `repulsion`,
  `clip` and `growth_threshold` were removed, as was the `"adaptive"` mode;
  `stability_radius`, `novelty_alpha`, `growth_intrusion`, `calibration_size`,
  `local_linear` and `geodesic` were added. Defaults now use 64 prototypes.
* Update records use the statuses `accepted`, `unchanged`, `frozen` and
  `rollback`.
* Distances and mapping are vectorized; the neighborhood audit maps in blocks.
* Validation errors name the part of the model that failed.

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
