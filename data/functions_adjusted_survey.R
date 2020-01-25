wtd.mean <- function(x, weights=NULL, normwt='ignored', na.rm=TRUE)
{
  if(! length(weights)) return(mean(x, na.rm=na.rm))
  if(na.rm) {
    s <- ! is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }

  sum(weights * x) / sum(weights)
}



wtd.var <- function(x, weights=NULL, normwt=FALSE, na.rm=TRUE,
                    method = c('unbiased', 'ML'))
  ## By Benjamin Tyner <btyner@gmail.com> 2017-0-12
{
  method <- match.arg(method)
  if(! length(weights)) {
    if(na.rm) x <- x[!is.na(x)]
    return(var(x))
  }
  
  if(na.rm) {
    s       <- !is.na(x + weights)
    x       <- x[s]
    weights <- weights[s]
  }

  if(normwt)
    weights <- weights * length(x) / sum(weights)

  if(normwt || method == 'ML')
    return(as.numeric(stats::cov.wt(cbind(x), weights, method = method)$cov))

  # the remainder is for the special case of unbiased frequency weights
  sw  <- sum(weights)
  if(sw <= 1)
      warning("only one effective observation; variance estimate undefined")
  
  xbar <- sum(weights * x) / sw
  sum(weights*((x - xbar)^2)) / (sw - 1)
}

wtd.quantile <- function(x, weights=NULL, probs=c(0, .25, .5, .75, 1), 
                         type=c('quantile','(i-1)/(n-1)','i/(n+1)','i/n'), 
                         normwt=FALSE, na.rm=TRUE)
{
  if(! length(weights))
    return(quantile(x, probs=probs, na.rm=na.rm))

  type <- match.arg(type)
  if(any(probs < 0 | probs > 1))
    stop("Probabilities must be between 0 and 1 inclusive")

  nams <- paste(format(round(probs * 100, if(length(probs) > 1) 
                             2 - log10(diff(range(probs))) else 2)), 
                "%", sep = "")

  i <- is.na(weights) | weights == 0
  if(any(i)) {
    x <- x[! i]
    weights <- weights[! i]
    }
  if(type == 'quantile') {
    w <- wtd.table(x, weights, na.rm=na.rm, normwt=normwt, type='list')
    x     <- w$x
    wts   <- w$sum.of.weights
    n     <- sum(wts)
    order <- 1 + (n - 1) * probs
    low   <- pmax(floor(order), 1)
    high  <- pmin(low + 1, n)
    order <- order %% 1
    ## Find low and high order statistics
    ## These are minimum values of x such that the cum. freqs >= c(low,high)
    allq <- approx(cumsum(wts), x, xout=c(low,high), 
                   method='constant', f=1, rule=2)$y
    k <- length(probs)
    quantiles <- (1 - order)*allq[1:k] + order*allq[-(1:k)]
    names(quantiles) <- nams
    return(quantiles)
  } 
  w <- wtd.Ecdf(x, weights, na.rm=na.rm, type=type, normwt=normwt)
  structure(approx(w$ecdf, w$x, xout=probs, rule=2)$y, 
            names=nams)
}


wtd.Ecdf <- function(x, weights=NULL, 
                     type=c('i/n','(i-1)/(n-1)','i/(n+1)'), 
                     normwt=FALSE, na.rm=TRUE)
{
  type <- match.arg(type)
  switch(type,
         '(i-1)/(n-1)'={a <- b <- -1},
         'i/(n+1)'    ={a <- 0; b <- 1},
         'i/n'        ={a <- b <- 0})

  if(! length(weights)) {
    ##.Options$digits <- 7  ## to get good resolution for names(table(x))
    oldopt <- options('digits')
    options(digits=7)
    on.exit(options(oldopt))
    cumu <- table(x)    ## R does not give names for cumsum
    isdate <- testDateTime(x)  ## 31aug02
    ax <- attributes(x)
    ax$names <- NULL
    x <- as.numeric(names(cumu))
    if(isdate) attributes(x) <- c(attributes(x),ax)
    cumu <- cumsum(cumu)
    cdf <- (cumu + a)/(cumu[length(cumu)] + b)
    if(cdf[1]>0) {
      x <- c(x[1], x);
      cdf <- c(0,cdf)
    }

    return(list(x = x, ecdf=cdf))
  }

  w <- wtd.table(x, weights, normwt=normwt, na.rm=na.rm)
  cumu <- cumsum(w$sum.of.weights)
  cdf <- (cumu + a)/(cumu[length(cumu)] + b)
  list(x = c(if(cdf[1]>0) w$x[1], w$x), ecdf=c(if(cdf[1]>0)0, cdf))
}


wtd.table <- function(x, weights=NULL, type=c('list','table'), 
                      normwt=FALSE, na.rm=TRUE)
{
  type <- match.arg(type)
  if(! length(weights))
    weights <- rep(1, length(x))

  isdate <- testDateTime(x)  ## 31aug02 + next 2
  ax <- attributes(x)
  ax$names <- NULL
  
  if(is.character(x)) x <- as.factor(x)
  lev <- levels(x)
  x <- unclass(x)
  
  if(na.rm) {
    s <- ! is.na(x + weights)
    x <- x[s, drop=FALSE]    ## drop is for factor class
    weights <- weights[s]
  }

  n <- length(x)
  if(normwt)
    weights <- weights * length(x) / sum(weights)

  i <- order(x)  # R does not preserve levels here
  x <- x[i]; weights <- weights[i]

  if(anyDuplicated(x)) {  ## diff(x) == 0 faster but doesn't handle Inf
    weights <- tapply(weights, x, sum)
    if(length(lev)) {
      levused <- lev[sort(unique(x))]
      if((length(weights) > length(levused)) &&
         any(is.na(weights)))
        weights <- weights[! is.na(weights)]

      if(length(weights) != length(levused))
        stop('program logic error')

      names(weights) <- levused
    }

    if(! length(names(weights)))
      stop('program logic error')

    if(type=='table')
      return(weights)

    x <- all.is.numeric(names(weights), 'vector')
    if(isdate)
      attributes(x) <- c(attributes(x),ax)

    names(weights) <- NULL
    return(list(x=x, sum.of.weights=weights))
  }

  xx <- x
  if(isdate)
    attributes(xx) <- c(attributes(xx),ax)

  if(type=='list')
    list(x=if(length(lev))lev[x]
           else xx, 
         sum.of.weights=weights)
  else {
    names(weights) <- if(length(lev)) lev[x]
                      else xx
    weights
  }
}


wtd.rank <- function(x, weights=NULL, normwt=FALSE, na.rm=TRUE)
{
  if(! length(weights))
    return(rank(x, na.last=if(na.rm) NA else TRUE))

  tab <- wtd.table(x, weights, normwt=normwt, na.rm=na.rm)
  
  freqs <- tab$sum.of.weights
  ## rank of x = # <= x - .5 (# = x, minus 1)
  r <- cumsum(freqs) - .5*(freqs-1)
  ## Now r gives ranks for all unique x values.  Do table look-up
  ## to spread these ranks around for all x values.  r is in order of x
  approx(tab$x, r, xout=x)$y
}


wtd.loess.noiter <- function(x, y, weights=rep(1,n),
                             span=2/3, degree=1, cell=.13333, 
                             type=c('all','ordered all','evaluate'), 
                             evaluation=100, na.rm=TRUE) {
  type <- match.arg(type)
  n <- length(y)
  if(na.rm) {
    s <- ! is.na(x + y + weights)
    x <- x[s]; y <- y[s]; weights <- weights[s]; n <- length(y)
  }
  
  max.kd <- max(200, n)
  # y <- stats:::simpleLoess(y, x, weights=weights, span=span,
  #                          degree=degree, cell=cell)$fitted
  y <- fitted(loess(y ~ x, weights=weights, span=span, degree=degree,
		control=loess.control(cell=cell, iterations=1)))

  switch(type,
         all=list(x=x, y=y),
         'ordered all'={
           i <- order(x);
           list(x=x[i],y=y[i])
         },
         evaluate={
           r <- range(x, na.rm=na.rm)
           approx(x, y, xout=seq(r[1], r[2], length=evaluation))
         })
}

num.denom.setup <- function(num, denom)
{
  n <- length(num)
  if(length(denom) != n)
    stop('lengths of num and denom must match')
  
  s <- (1:n)[! is.na(num + denom) & denom != 0]
  num <- num[s];
  denom <- denom[s]
  
  subs <- s[num > 0]
  y <- rep(1, length(subs))
  wt <- num[num > 0]
  other <- denom - num
  subs <- c(subs, s[other > 0])
  wt <- c(wt, other[other > 0])
  y <- c(y, rep(0, sum(other>0)))
  list(subs=subs, weights=wt, y=y)
}
                
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
