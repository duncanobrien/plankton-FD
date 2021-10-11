## Round 2 Function (renamed 'roundTO' to prevent homonyms) ##

roundTO <- function(x,y){
#roundTO rounds number to nearest multiple of arbitrary precision.
#Z <- roundTO(X,Y) rounds X to nearest multiple of Y.
round(x/y)*y
}
