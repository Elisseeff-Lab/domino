# Since some of the vignettes either take some time or prompt strange
# issues when built with pkgdown, we can precompile them instead
# (We do have to make sure they get updated though, so here is a script
# to assist with that):

knitr::knit("domino2.Rmd.orig", "domino2.Rmd")
knitr::knit("plotting_vignette.Rmd.orig", "plotting_vignette.Rmd")
knitr::knit("domino_object_vignette.Rmd.orig", "domino_object_vignette.Rmd")
