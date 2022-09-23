.logistic_growth_stat <- function(u, ta = NULL, generation_time_scale = Tg)
{
  if ( is.null(ta))
    ta = .get_comparator_sample(u) 
  if ( is.null( ta )){ #if still null, cant get good comparison 
    return( 
      list( lgr = NA, lgrp = .5, gam_r = NA
            , AIC = NA
            , AIC_gam = NA
            , growthrates = setNames( c(NA, NA), c('Logistic', 'GAM') )
            , relative_model_support = NA
            , plot = NULL 
      )
    )
  }
  tu = descendantSids[[u]] 
  sta = sts[ ta ]
  stu = sts[ tu ]
  X = data.frame( time = c( sta, stu ), type = c( rep('control', length(ta)), rep('clade',length(tu)) ) )
  X = na.omit( X ) 
  m = glm( type=='clade' ~ time, data = X, family = binomial( link = 'logit' ))
  s = summary( m ) 
  rv = unname( coef( m )[2] * generation_time_scale ) 
  p = NA 
  if ( is.na( rv )){
    message( 'NA growth stat, node: ', u  )	
  } else{ 
    p = s$coefficients[2, 4 ]
  }
  ## time dep growth ; needs a larger sample size 
  X$estimated = predict( m )
  if ( compute_gam & (length( tu ) > 50) ){
    m1 = mgcv::gam( type=='clade' ~ s(time, bs = "bs", k = 4, m=1) , family = binomial(link='logit') , data = X)
    X$estimated = predict( m1 )
    
    tout = seq( min(X$time) , max(X$time), length=5)
    tout1 = tout[4] + diff(tout)[1]/2 
    dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
    r = dlo*Tg / ((max(tout) - tout1))  
    aic = c(AIC(m), AIC(m1))
  } else{
    r = NA 
    aic = c(AIC(m), Inf)
  }
  
  list( lgr = rv, lgrp = p, gam_r =r 
        , AIC = aic[1]
        , AIC_gam = aic[2]
        , growthrates = setNames( c(rv, r), c('Logistic', 'GAM') )
        , relative_model_support = setNames( exp(  min(aic) - aic  )/2 ,  c('Logistic', 'GAM') )
        , plot = .freq_figure( X ) 
  )
}




m1 = mgcv::gam( type=='clade' ~ s(time, bs = "bs", k = 4, m=1) , family = binomial(link='logit') , data = X)
X$estimated = predict( m1 )

tout = seq( min(X$time) , max(X$time), length=5)
tout1 = tout[4] + diff(tout)[1]/2 
dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
r = dlo*Tg / ((max(tout) - tout1))  
aic = c(AIC(m), AIC(m1))

k1 <- list( lgr = rv, lgrp = p, 
         gam_r =r 
      , AIC = aic[1]
      , AIC_gam = aic[2]
      , growthrates = setNames( c(rv, r), c('Logistic', 'GAM') )
      , relative_model_support = setNames( exp(  min(aic) - aic  )/2 ,  c('Logistic', 'GAM') )
      , plot = .freq_figure( X ) 
      )

k2 <- list( lgr = rv, lgrp = p, 
            gam_r =r 
            , AIC = aic[1]
            , AIC_gam = aic[2]
            , growthrates = setNames( c(rv, r), c('Logistic', 'GAM') )
            , relative_model_support = setNames( exp(  min(aic) - aic  )/2 ,  c('Logistic', 'GAM') )
            , plot = .freq_figure( X ) 
)




# fread vs benchmark testing ####

benchmark(
  "bs" = {
    
  m1 = mgcv::gam( type=='clade' ~ s(time, bs = "bs", k = 4, m=1) , family = binomial(link='logit') , data = X)
  X$estimated = predict( m1 )
  tout = seq( min(X$time) , max(X$time), length=5)
  tout1 = tout[4] + diff(tout)[1]/2 
  dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
  r = dlo*Tg / ((max(tout) - tout1))  
  aic = c(AIC(m), AIC(m1))
  
  },
"cr" = {
  m1 = mgcv::gam( type=='clade' ~ s(time, bs = "cr", k = 4, m=1) , family = binomial(link='logit') , data = X)
  X$estimated = predict( m1 )
  tout = seq( min(X$time) , max(X$time), length=5)
  tout1 = tout[4] + diff(tout)[1]/2 
  dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
  r = dlo*Tg / ((max(tout) - tout1))  
  aic = c(AIC(m), AIC(m1))
  
},
"cs" = {
  m1 = mgcv::gam( type=='clade' ~ s(time, bs = "cs", k = 4, m=1) , family = binomial(link='logit') , data = X)
  X$estimated = predict( m1 )
  tout = seq( min(X$time) , max(X$time), length=5)
  tout1 = tout[4] + diff(tout)[1]/2 
  dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
  r = dlo*Tg / ((max(tout) - tout1))  
  aic = c(AIC(m), AIC(m1))
  
},


"bs_bam" = {
  
  m1 = mgcv::bam( type=='clade' ~ s(time, bs = "bs", k = 4, m=1) , family = binomial(link='logit') , data = X)
  X$estimated = predict( m1 )
  tout = seq( min(X$time) , max(X$time), length=5)
  tout1 = tout[4] + diff(tout)[1]/2 
  dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
  r = dlo*Tg / ((max(tout) - tout1))  
  aic = c(AIC(m), AIC(m1))
  
},


"cr_bam" = {
  
  m1 = mgcv::bam( type=='clade' ~ s(time, bs = "cr", k = 4, m=1) , family = binomial(link='logit') , data = X)
  X$estimated = predict( m1 )
  tout = seq( min(X$time) , max(X$time), length=5)
  tout1 = tout[4] + diff(tout)[1]/2 
  dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
  r = dlo*Tg / ((max(tout) - tout1))  
  aic = c(AIC(m), AIC(m1))
  
},


"cs_bam" = {
  
  m1 = mgcv::bam( type=='clade' ~ s(time, bs = "cs", k = 4, m=1) , family = binomial(link='logit') , data = X)
  X$estimated = predict( m1 )
  tout = seq( min(X$time) , max(X$time), length=5)
  tout1 = tout[4] + diff(tout)[1]/2 
  dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
  r = dlo*Tg / ((max(tout) - tout1))  
  aic = c(AIC(m), AIC(m1))
  
},

replications = 5,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self")
)


#     test replications elapsed relative user.self sys.self
#2   fread            3   27.26    1.000     36.31     4.06
#1 readcsv            3  680.91   24.978    669.21     3.87


"bs_bam" = {
  
  m1 = mgcv::bam( type=='clade' ~ s(time, bs = "bs", k = 4, m=1) , family = binomial(link='logit') , data = X)
  X$estimated = predict( m1 )
  tout = seq( min(X$time) , max(X$time), length=5)
  tout1 = tout[4] + diff(tout)[1]/2 
  dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
  r = dlo*Tg / ((max(tout) - tout1))  
  aic = c(AIC(m), AIC(m1))
  
},


"bs_bam" = {
  
  m1 = mgcv::bam( type=='clade' ~ s(time, bs = "cr", k = 4, m=1) , family = binomial(link='logit') , data = X)
  X$estimated = predict( m1 )
  tout = seq( min(X$time) , max(X$time), length=5)
  tout1 = tout[4] + diff(tout)[1]/2 
  dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
  r = dlo*Tg / ((max(tout) - tout1))  
  aic = c(AIC(m), AIC(m1))
  
},


"bs_bam" = {
  
  m1 = mgcv::bam( type=='clade' ~ s(time, bs = "cs", k = 4, m=1) , family = binomial(link='logit') , data = X)
  X$estimated = predict( m1 )
  tout = seq( min(X$time) , max(X$time), length=5)
  tout1 = tout[4] + diff(tout)[1]/2 
  dlo = diff( predict( m1, newdata = data.frame( type=NA, time = c( tout1 , max(tout))) )  )
  r = dlo*Tg / ((max(tout) - tout1))  
  aic = c(AIC(m), AIC(m1))
  
},

