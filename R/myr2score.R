myr2score=function(m,o){
mean_y<-mean(m)
ss_tot <- sum((m - mean_y) ^ 2)
ss_res = sum((m - o) ^ 2)
r2 <- 1 - (ss_res / ss_tot)
}
