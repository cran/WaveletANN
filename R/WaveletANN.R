#=====================================================================================#
# PURPOSE : Application 0f Wavelet-ANN hybrid model for forecasting time series     #
# AUTHOR  : Ranjit Kumar Paul                                    #
# DATE    : 16 March, 2019                                                          #
# VERSION : Ver 0.1.0                                                                 #
#=====================================================================================#

#---------------------------------------------------------------------------------------#
# Computing Wavelet Coefficients using MODWT algorithm using haar filter                #
#---------------------------------------------------------------------------------------#

WaveletFitting <- function(ts,Wvlevels,bndry,FFlag)
{
  mraout <- wavelets::modwt(ts, filter='haar', n.levels=Wvlevels,boundary=bndry, fast=FFlag)
  WaveletSeries <- cbind(do.call(cbind,mraout@W),mraout@V[[Wvlevels]])
  return(list(WaveletSeries=WaveletSeries,WVSeries=mraout))
}

WaveletFittingann<- function(ts,Waveletlevels,boundary,FastFlag,nonseaslag,seaslag,NForecast)

{
  WS <- WaveletFitting(ts=ts,Wvlevels=Waveletlevels,bndry=boundary,FFlag=FastFlag)$WaveletSeries
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  #-----------------------------------------------------------#
  # Fitting of ANN model to the Wavelet Coef                #
  #-----------------------------------------------------------#
  for(WVLevel in 1:ncol(WS))
  {
    ts <- NULL
    ts <- WS[,WVLevel]
    WaveletANNFit <- forecast::nnetar(y=as.ts(ts))

    WaveletANNPredict <- WaveletANNFit$fitted
    WaveletANNForecast <- forecast::forecast(WaveletANNFit,h=NForecast)
    AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletANNPredict)
    AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletANNForecast$mean))
  }
  Finalforecast <- rowSums(AllWaveletForecast,na.rm = T)
  FinalPrediction <- rowSums(AllWaveletPrediction,na.rm = T)
  return(list(Finalforecast=Finalforecast,FinalPrediction=FinalPrediction))
}
