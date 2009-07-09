`p.control` <-
function(DL=FALSE,diff.varying=FALSE,ridge.formulas=list(cold="xlag^2",heat="xlag^2")){
#auxiliary function penalty options
              list(DL=DL,diff.varying=diff.varying, ridge.formulas=ridge.formulas)
              }

