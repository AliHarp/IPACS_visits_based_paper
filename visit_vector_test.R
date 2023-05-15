visit_vector <- function(){
  v <- (seq(4.25,8.25,length.out = 8)) 
  if (is.integer(v) && length(v) == 0L){
    v <- 0}
  return(v)}
visit_vector <- visit_vector()
#remove rounding bias
for (i in 1:length(visit_vector)){
  visit_vector[i] <- round(visit_vector[i]/0.5)*0.5
}
print(visit_vector)