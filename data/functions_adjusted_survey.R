compute_share <- function(qw,asset,surweight){
  asset_Q = wtd.quantile(asset, probs=qw, na.rm = TRUE, weights=surweight)
  return(sum(surweight*asset*(asset>=asset_Q),na.rm=T)/sum(surweight*asset,na.rm=T))
}


compute_total <- function(qw,asset,surweight){
  asset_Q = wtd.quantile(asset, probs=qw, na.rm = TRUE, weights=surweight)
  return(sum(surweight*asset*(asset>=asset_Q),na.rm=T))
}


compute_DATA <- function(asset,surweight){
  return(sum(surweight*asset,na.rm=T))
}


# ### function to get the quantile
# pareto_shape = 1.55
# qw=0.99
get_quantile_wealth <- function(sum_top,networth,surweight,surweight_tail,min_threshold,pareto_shape,qw){
  
  yy = (qw - (1-sum_top))/sum_top
  
  if(yy <= 0){
    wealth_q = wtd.quantile(networth, probs=qw, na.rm = TRUE, weights=surweight)
  }else{
    wealth_q = qpareto(yy,min_threshold,pareto_shape) 
  }
  
  
  ## total survey wealth
  tot_surv_w = sum(networth*surweight, na.rm = TRUE)
  ## total tail wealth above 99%
  tot_surv_w_q = sum(networth[networth >= wealth_q]*surweight[networth >= wealth_q], na.rm = TRUE)
  ## share of wealth from survey
  frac_q_surv = tot_surv_w_q / tot_surv_w * 100
  
  ## using the pareto tail
  tot_surv_w_non_pareto_part = sum(networth[networth < min_threshold]*surweight[networth < min_threshold], na.rm = TRUE)
  tot_pareto_part_w = min_threshold*pareto_shape*sum(surweight_tail)/(pareto_shape-1)
  tot_pareto_w = tot_surv_w_non_pareto_part+tot_pareto_part_w
  
  if(wealth_q > min_threshold){
    CDF = 1 - (min_threshold/wealth_q)^(pareto_shape)
    frac_w_q_est = ((tot_pareto_part_w*(1-CDF)^(1-(1/pareto_shape)))/tot_pareto_w)*100
  }else{
    frac_w_q_est = (1-(sum(networth[networth < wealth_q]*surweight[networth < wealth_q], na.rm = TRUE))/tot_pareto_w)*100
  }
  
  return(c(frac_w_q_est,frac_q_surv))
}






#### Pareto discretize function
pareto_discretize_fun <- function(pareto_shape,wmin){
  
  ## step 1
  ## construct vector of wealth equally spaced on the CDF of the pareto
  # F(x) = z = 1-(wmin/w)^a
  z = seq(0,0.9999,by=0.0001)
  wealth = wmin/(1-z)^(1/pareto_shape)
  
  ## step 2
  # construct (w1+w2)/2, (w2+w3)/2,...(w9,999+w10,000)/2, w10,000
  #wealth_vector = (wealth[1:9999] + wealth[2:10000])/2
  wealth_vector = wealth
  
  
  ## step 3 construct probability
  probability_vector = rep(1,10000)*0.0001
  
  return(as.data.frame(cbind(wealth_vector, probability_vector)))
}






### function to get the adjusted pareto tail 
# nb_rich_list     = nrow(bill_data_current_y)
# exchange_rate_US = OECD_data$ex_rate_US[OECD_data$year == list_year[p] & OECD_data$isocode == list_countries[p]]
# networth = networth_UA
# min_threshold = wmin
# country_name = "FRA"
# year_selected = 2014
pareto_estimate <- function(country_name,min_threshold,year_selected,networth,surweight,exchange_rate_US,nb_rich_list){
  
  HFCS_all = as.data.frame(cbind(surweight,networth))
  HFCS_temp = HFCS_all[which(HFCS_all$networth >= min_threshold),]
  
  # table(HFCS_temp$sa0100 )
  # table(HFCS_temp$ra0010 )
  # table(HFCS_temp$sa0010, HFCS_temp$ra0300)
  # HFCS_temp = HFCS_temp[,c("surweight","networth")]
  # HFCS_temp_networth = aggregate(networth ~ sa0010,data=HFCS_temp,FUN=mean)
  # HFCS_temp_surweight = aggregate(surweight ~ sa0010,data=HFCS_temp,FUN=sum)
  # HFCS_temp = merge(HFCS_temp_networth,HFCS_temp_surweight,by="sa0010")
  # HFCS_temp = HFCS_temp[,-1]
  HFCS_temp = HFCS_temp[order(HFCS_temp$networth),]
  rlist = bill_data[which(bill_data$isocode == country_name & bill_data$year == as.numeric(year_selected)),]
  
  ## computing empirical ccdf
  HFCS_temp$cumul_weight = (HFCS_temp$surweight)[length(HFCS_temp$surweight)]
  for(i in (nrow(HFCS_temp)-1):1){
    HFCS_temp$cumul_weight[i] = HFCS_temp$surweight[i] + HFCS_temp$cumul_weight[i+1]
  }
  total_weight = (HFCS_temp$cumul_weight)[1] + nrow(rlist)
  (HFCS_temp$surweight)[length(HFCS_temp$cumul_weight)]
  HFCS_temp$cumul_weight = HFCS_temp$cumul_weight/total_weight
  range(HFCS_temp$cumul_weight)
  # plot(HFCS_temp$networth,HFCS_temp$cumul_weight,log="yx")
  
  
  
  ## Forbes observations
  if(nb_rich_list > 1){
  if(country_name == "FR"){rlist[rlist == "7,668.00"] <- 7.668}
  
    rich = as.numeric(as.character(rlist$networth))*1000000000*exchange_rate_US
    rich = sort(rich)
    rich = as.data.frame(rich)
    rich$weight = 1
    rich$cumul_weight = 1
    for(i in (nrow(rich)-1):1){
      rich$cumul_weight[i] = rich$weight[i+1] + rich$cumul_weight[i+1] 
    }
    rich$cumul_weight = rich$cumul_weight/total_weight
    colnames(rich) <- c("networth","surweight","cumul_weight")
  
    tot_list = rbind(HFCS_temp,rich)
    tot_list$cumul_weight = tot_list$cumul_weight + 0.00000001 
    range(tot_list$networth)
    range(tot_list$cumul_weight)
  }
  
  
  # par(mar=c(4,4,2,2))
  # plot(HFCS_temp$networth/wmin,HFCS_temp$cumul_weight,log="yx", col="darkblue", frame.plot=FALSE, ylim=c(0.0000000001,200),xlim=c(1,10000), xlab="wealth (in million euros)", ylab="Empirical CCDF")
  # par(new=TRUE)
  # plot(rich$networth/wmin,rich$cumul_weight,log="yx", col="darkred", frame.plot=FALSE, pch="+", xlab="", ylab="", axes=FALSE, ylim=c(0.0000000001,200),xlim=c(1,10000))
  
  
  ###### ESTIMATION OF THE POWER LAW using the data.
  ## literature: Gabaix (2009), Clauset et al. (2009)
  ## we estimate :: ln(frequency) = -alpha ln(wi/wmin)
  
  ## regression for France:
  if(nb_rich_list > 1){model_for = (lm(log(cumul_weight) ~ log(networth/wmin) - 1, data= tot_list))}
  model_nofor = (lm(log(cumul_weight) ~ log(networth/wmin) - 1, data= HFCS_temp))
  
  # ## mle regressio
  # model_mle = dpareto.ll(tot_list$networth*tot_list$cumul_weight,theta=0.5)
  # model_mle = epareto(tot_list$networth*tot_list$cumul_weight)
  # model_mle = pareto2_estimate_mle(tot_list$networth*tot_list$cumul_weight,s = min_threshold, smin=0.0001,smax=3)
  # model_mle = fgpd(tot_list$networth*tot_list$cumul_weight,pvector=c(sigmau= 1000000,1.5))
  
  
  # reg_wealth = seq(wmin,10000000000, by=wmin)
  # reg_freq = (reg_wealth/wmin)^(model_nofor$coefficients)
  # reg_freq_for = (reg_wealth/wmin)^(model_for$coefficients)
  
  # par(mar=c(4,4,2,2))
  # plot(HFCS_temp$networth/wmin,HFCS_temp$cumul_weight,log="yx", col="darkblue", frame.plot=FALSE, ylim=c(0.0000001,20),xlim=c(1,10000), xlab="wealth (in million euros)", ylab="Empirical CCDF")
  # par(new=TRUE)
  # plot(rich$networth/wmin,rich$cumul_weight,log="yx", col="darkred", frame.plot=FALSE, pch="+", xlab="", ylab="", axes=FALSE, ylim=c(0.0000001,20),xlim=c(1,10000))
  # par(new=TRUE)
  # plot(reg_wealth/wmin,reg_freq,log="yx", col="black", frame.plot=FALSE, type="l", xlab="", ylab="", axes=FALSE, ylim=c(0.0000001,20),xlim=c(1,10000))
  # par(new=TRUE)
  # plot(reg_wealth/wmin,reg_freq_for,log="yx", col="black", frame.plot=FALSE, type="l", lty=5,  xlab="", ylab="", axes=FALSE, ylim=c(0.0000001,20),xlim=c(1,10000))
  
  if(nb_rich_list > 1){estimated_pareto = c(summary(model_nofor)$coef[1],summary(model_for)$coef[1],summary(model_nofor)$coef[2],summary(model_for)$coef[2])
  }else{estimated_pareto = c(summary(model_nofor)$coef[1],summary(model_nofor)$coef[1],summary(model_nofor)$coef[2],summary(model_nofor)$coef[2])}
  
  
  ## get the total weight in the pareto tails:
  sum_top = sum(HFCS_temp$surweight)/sum(HFCS_all$surweight)
  surweight_tail = (HFCS_temp$surweight)

  if(nb_rich_list > 1){
    frac_99 = get_quantile_wealth(sum_top,networth,surweight,surweight_tail,min_threshold,-summary(model_for)$coef[1],0.99)
    frac_95 = get_quantile_wealth(sum_top,networth,surweight,surweight_tail,min_threshold,-summary(model_for)$coef[1],0.95)
  }else{
    frac_99 = get_quantile_wealth(sum_top,networth,surweight,surweight_tail,min_threshold,-summary(model_nofor)$coef[1],0.99)
    frac_95 = get_quantile_wealth(sum_top,networth,surweight,surweight_tail,min_threshold,-summary(model_nofor)$coef[1],0.95)    
  }
  
  return(c(estimated_pareto,frac_99,frac_95))
}
