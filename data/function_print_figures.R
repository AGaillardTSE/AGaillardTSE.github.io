 
install.packages("RCurl",lib="/media/user", repos='http://cran.us.r-project.org') 
library("RCurl", lib.loc="/media/user") 

ftpUpload('data_LWS_adjusted.csv', 
          "ftp://univnantes:VW29DTUNBf4DLdx@ftpperso.free.fr/data_LWS/data_wealth_concentration_2.csv",port=21,verbose=TRUE) 
 
ftpUpload('data_LWS_adjusted_JPN.csv', 
          "ftp://univnantes:VW29DTUNBf4DLdx@ftpperso.free.fr/data_LWS/data_wealth_concentration_2_JPN.csv",port=21,verbose=TRUE) 
 
