library(fitdistrplus)
# set files / data to use
files <- c(
   'https://climexp.knmi.nl/data/xgdcnCA001114745.dat', # 1921 - 1949
   'https://climexp.knmi.nl/data/xgdcnCA001114740.dat', # 1944 - 1969
   'https://climexp.knmi.nl/data/xgdcnCA001114741.dat', # 1970 - 1991
   'https://climexp.knmi.nl/data/xgdcnCA007034395.dat', # 1973 - 1994
   'https://climexp.knmi.nl/data/xgdcnCA001114739.dat', # 1990 - 2012
   'https://climexp.knmi.nl/data/xgdcnCA001114746.dat' # 2006 - present
)
# set up variable to store the data
d <- NULL
# now iterate over the files and read the data
for( idx in 1:length(files) ){
   f <- files[idx]
   tmpin <-read.fwf( url( f ), skip=19, widths=c(5,3,3,10) )
   tmpin$file <- rep( idx, nrow(tmpin))
   if( is.null((d))){
      d <- tmpin
   }else{
      d <- rbind( d, tmpin)
   }
}
# set variable names, get date and day of year
names( d ) <- c('year','month','day','tmax','file')
d$date <- as.POSIXct( paste(d$year,'-',d$month,'-',d$day,sep='') )
d$jday <- format.POSIXct(d$date,'%j')
d$jday <- as.numeric(d$jday)
d$tmax <- d$tmax + 273.15 # change to K to avoid negative numbers for fitting gamma dist.


# now do the calculation
window <- 2 # use 5 day window (date ± 2 days)
ulim_all <- rep(NA, 365) # vector to store Tmax upper limit for all stations
ulim_current <- rep(NA, 365) # vector to store Tmax upper limit for current / last stations

dist <- 'gamma'

# iterate over days getting 99.99% quantile
for( jdayidx in seq(window + 1, 365 - window, 1) ){
   # all data
   ss <- subset(d, abs(jday - jdayidx) <= window & year < 2021)
   fg <- fitdist( ss$tmax, distr=dist, method='mle')
   if( dist == 'norm' ){
      ulim_all[ jdayidx ] <- qnorm( 0.9999, mean = fg$estimate[1], sd = fg$estimate[2] )
   }else if ( dist == 'gamma'){
      ulim_all[ jdayidx ] <- qgamma( 0.9999, shape = fg$estimate[1], rate = fg$estimate[2] )   
   }else{
      ulim_all[ jdayidx ] <- quantile( ss$tmax, 0.9999)
   }
   # current station
   ss <- subset(d, abs(jday - jdayidx) <= window & year < 2021 & file == 6)
   fg <- fitdist( ss$tmax, distr=dist, method='mle')
   if( dist == 'norm' ){
      ulim_current[ jdayidx ] <- qnorm( 0.9999, mean = fg$estimate[1], sd = fg$estimate[2] )
   }else if ( dist == 'gamma'){
      ulim_current[ jdayidx ] <- qgamma( 0.9999, shape = fg$estimate[1], rate = fg$estimate[2] )   
   }else{
      ulim_current[ jdayidx ] <- quantile( ss$tmax, 0.9999)  
   }
}
# now plot
plot( 1:365, ulim_all-273.15 , type='l', xlab = "Day of year", ylab = "Tmax", ylim=c(0,60), main='Tmax extremes (99.99 percentile)\n5 day window (day of year ± 2 days)')
ss <- subset(d, year == 2021)
lines( 1:365, ulim_current-273.15, col='green' )
lines( (tmax-273.15) ~ jday, ss, col='red')
legend( x = 0, y= 60, legend = c("1921 - 2020 (all stations)",'2006 - 2020 (current station)', '2021'), col=c('black','green','red'), lty = c(1,1))