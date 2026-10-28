# Interactive app: shiny::runApp() from this folder, or deploy with rsconnect::deployApp().
# For deployment, install ceakit from GitHub first:
#   remotes::install_github("fali-ebm/Galleri", subdir = "ceakit")
# so the server can install it too.

# Lets deployment tools detect packages; never executed.
if (FALSE) { library(shiny); library(bslib); library(DT); library(markdown); library(mgcv); library(ceakit) }

ceakit::cea_app(".")   # must be the last line
