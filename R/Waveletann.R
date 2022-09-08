#' Wavelet Transform Using Maximal Overlap Discrete Wavelet Transform (MODWT) Algorithm
#'
#' @param ts Univariate time series
#' @param Wvlevels The level of wavelet decomposition
#' @param Filter Wavelet filter
#' @param bndry The boundary condition of wavelet decomposition
#' @param FFlag The FastFlag condition of wavelet decomposition: True or False
#' @import fracdiff forecast stats wavelets Metrics
#' @return
#' \itemize{
#'   \item WaveletSeries - The wavelet transform of the series
#' }
#' @export
#'
#' @examples
#' data<-rnorm(100,mean=100,sd=50)
#' WaveletFitting(ts=data,Wvlevels=3,Filter='haar',bndry='periodic',FFlag=TRUE)
#' @references
#' \itemize{
#'\item Aminghafari, M. and Poggi, J.M. 2007. Forecasting time series using wavelets. Internationa Journal of Wavelets, Multiresolution and Inforamtion Processing, 5, 709 to 724

#' \item Percival D. B. and Walden A. T. 2000. Wavelet Methods for Time-Series Analysis. Cambridge Univ. Press, U.K.

#' \item Paul R. K., Prajneshu and Ghosh H. 2013. Wavelet Frequency Domain Approach for Modelling and Forecasting of Indian Monsoon Rainfall Time-Series Data. Journal of the Indian society of agricultural statistics, 67, 319 to 327.
#' }

WaveletFitting <- function(ts,Wvlevels,Filter='haar',bndry='periodic',FFlag=TRUE)
{
  mraout <- wavelets::modwt(ts, filter=Filter, n.levels=Wvlevels,boundary=bndry, fast=FFlag)
  WaveletSeries <- cbind(do.call(cbind,mraout@W),mraout@V[[Wvlevels]])
  return(list(WaveletSeries=WaveletSeries,WVSeries=mraout))
}

#' Wavelet-ANN Hybrid Model for Forecasting
#'
#' @param ts Univariate time series
#' @param Waveletlevels The level of wavelet decomposition
#' @param Filter Wavelet filter
#' @param boundary The boundary condition of wavelet decomposition
#' @param FastFlag The FastFlag condition of wavelet decomposition: True or False
#' @param nonseaslag Number of non seasonal lag
#' @param seaslag Number of non seasonal lag
#' @param hidden Size of the hidden layer
#' @param NForecast The forecast horizon: A positive integer
#' @import fracdiff forecast stats wavelets Metrics
#'
#' @return
#' \itemize{
#'   \item Finalforecast - Forecasted value
#'   \item FinalPrediction - Predicted value of train data
#'   \item Accuracy - RMSE and MAPE for train data
#' }
#' @export
#'
#' @examples
#' N <- 100
#'PHI <- 0.2
#'THETA <- 0.1
#' SD <- 1
#'M <- 0
#'D <- 0.2
#' Seed <- 123

#'set.seed(Seed)
#'Sim.Series <- fracdiff::fracdiff.sim(n = N,ar=c(PHI),ma=c(THETA),d=D,rand.gen =rnorm,sd=SD,mu=M)
#'simts <- as.ts(Sim.Series$series)
#'WaveletForecast<-WaveletFittingann(ts=simts,Waveletlevels=3,Filter='d4',
#'nonseaslag=5,hidden=3,NForecast=5)
#' @references
#' \itemize{
#'\item Aminghafari, M. and Poggi, J.M. 2012. Nonstationary time series forecasting using wavelets and kernel smoothing. Communications in Statistics-Theory and Methods, 41(3),485-499.
#' \item Paul, R.K. A and Anjoy, P. 2018. Modeling fractionally integrated maximum temperature series in India in presence of structural break. Theory and Applied Climatology 134, 241–249.
#' }
WaveletFittingann<- function(ts,Waveletlevels,Filter='haar',boundary='periodic',FastFlag=TRUE,nonseaslag,seaslag=1,hidden,NForecast)

{
  Actual<-ts
  WS <- WaveletFitting(ts=ts,Wvlevels=Waveletlevels,Filter=Filter,bndry=boundary,FFlag=FastFlag)$WaveletSeries
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  #-----------------------------------------------------------#
  # Fitting of ANN model to the Wavelet Coef                #
  #-----------------------------------------------------------#
  for(WVLevel in 1:ncol(WS))
  {
    ts <- NULL
    ts <- WS[,WVLevel]
    WaveletANNFit <- forecast::nnetar(y=as.ts(ts), p=nonseaslag, P=seaslag, size=hidden)

    WaveletANNPredict <- WaveletANNFit$fitted
    WaveletANNForecast <- forecast::forecast(WaveletANNFit,h=NForecast)
    AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletANNPredict)
    AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletANNForecast$mean))
  }
  Finalforecast <- rowSums(AllWaveletForecast,na.rm = T)
  FinalPrediction <- rowSums(AllWaveletPrediction,na.rm = T)
  RMSE<-Metrics::rmse(Actual,FinalPrediction)
  MAPE<-Metrics::mape(Actual,FinalPrediction)
  Accuracy<-cbind(RMSE,MAPE)
  return(list(Finalforecast=Finalforecast,FinalPrediction=FinalPrediction, Accuracy=Accuracy))
}

