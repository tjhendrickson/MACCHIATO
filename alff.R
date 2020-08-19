library(signal)
n = 500 # Length of time series
TR = 0.8 # Repetition time
Xt = arima.sim(n = n, model = list(ar=.5)) # Simulate an autoregressive time series
# Apply a bandpass filter, 2nd order butterworth, bandpass (0.01, 0.08)
bf = butter(2,
            c(0.01,0.08)/(1/TR/2),
            type="pass")
Xt.filter = signal::filter(bf,Xt)

# Estimate the power spectrum using a smoothed periodogram. Smoothing span is 12% the length of the time series. (This is arbitrary, but works well in practice.)
power = spec.pgram(Xt.filter,
                   detrend=TRUE,
                   spans = floor((.12*n)/2)*rep(1,2),
                   plot=FALSE) 
hertz = (1/TR/2)*power$freq # Convert normalized frequencies to Hertz
# Calculate alff
alff = mean(sqrt(power$spec[which(hertz <= 0.08 & hertz >= 0.01)]))

# If you don't want to smooth across frequencies:
#power = spec.pgram(Xt,
#                   plot=FALSE) 

# fALFF does not work with the filtered time series
power = spec.pgram(Xt,
                   detrend=TRUE,
                   spans = floor((.12*n)/2)*rep(1,2),
                   plot=FALSE) 
falff = sum(sqrt(power$spec[which(hertz <= 0.08 & hertz >= 0.01)]))/sum(sqrt(power$spec))