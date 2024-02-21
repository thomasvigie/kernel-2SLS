# Hypothesis tests & power -----------

t_test <- function(betahat, se, para, level = 0.05)
{
  if(level <0 || level >1)
    error("the level should be between 0 and 1 !!")
  tvalue <- (betahat - para)/se 
  decision <- abs(tvalue) > qnorm(1 - level/2)
  decision <- ifelse(decision == TRUE, "reject", "not reject") 
  return(decision)
}
# Example
# t_test(betahat = est, se = se, para = 1, level = 0.05)


H0_test <- function(parameter, estimate, std_error, level = 0.05, df)
{
  name <- paste("H0(", parameter, ")", sep = "")
  df <- df %>% mutate(name = t_test(estimate, se = std_error, para = parameter, level = level))
  df[, ncol(df)] <- ifelse(df[, ncol(df)] == "reject", 1, 0) 
  names(df)[ncol(df)] <- name
  return(df[names(df)==name])
}

# listu <- grid %>% map(H0_test, estimate = listu$TV, std_error = listu$TV_se, level = 0.05, df = listu)

test_to_dat <- function(parameter, estimate, std_error, level = 0.05, df, nsim)
{
  dat <- parameter %>% map(H0_test, estimate = estimate, std_error = std_error, level = level, df = df)
  dat <- purrr::reduce(dat, cbind)
  
  dat <- cbind(df [, which(grepl("sample|varu|varv|covuv|fs_fun|bandwidth", names(df)) == TRUE)], dat)
  tests <- dat [, which(grepl("sample|varu|varv|covuv|fs_fun|bandwidth|H0", names(dat)) == TRUE)] %>%
    group_by( sample, varu, varv, covuv, fs_fun, bandwidth ) %>% 
    summarise_all( freq, S = nsim )  
  # nsim <- nrow(df) / nrow(tests)
  #tests[, which(grepl("H0", names(tests)) == TRUE)] <- tests[, which(grepl("H0", names(tests)) == TRUE)] / nsim
  return(tests)
}

# Power graphs (general estimators) ------------------------
# This function computes the rejection frequencies for each design, each bandwidth value and some Hypothesis H0
rej <- function(data, H0, est, se, level = 0.05, name = NULL)
{
  nsim <- max(data$simu)
  # Use left_join or right_join when adding the new power columns
  tests <- data %>% mutate( power = ifelse(t_test(betahat = est, se = se, para = H0, level = level) == "reject", 1, 0)) 
  
  rejections <- tests %>%dplyr::select(sample, varu, varv, covuv, fs_fun, bandwidth, ends_with("power") )%>%
    group_by(sample, varu, varv, covuv, fs_fun, bandwidth)%>%
    summarize(power = sum(power)/nsim)
  
  if(is.null(name) == FALSE)
  {
    names(rejections)[names(rejections) == 'power'] <- name
  }
  
  rejections <- rejections %>% mutate(H0 = H0)
  return(rejections)
}

# Example
# powa_ols <- rej(data = listu, H0 = 1, est = listu$ols, se = listu$ols_se, level = 0.05)
# powa_tsls <- rej(data = listu, H0 = 1, est = listu$tsls, se = listu$tsls_se, level = 0.05)


power_graph <- function(j, data, H0_range, level = 0.05, by = c("sample", "bandwidth"))
{
  
  # b <- length(bandrange) 
  # statistics <- c("Bias", "Variance", "MSE")
  #  da <- list()
  
  # for (i in 1:b){
  # for ( i in 3){      # when making the zooms of the MSE ones
  #   data_to_graph <- MC[[i]] [ , -(1:5)]
  #   data_to_graph <- dplyr::filter(MC[[i]], sample == unique_designs[j,]$sample, 
  #                                  varu == unique_designs[j ,]$varu, 
  #                                  fs_fun == unique_designs[j ,]$fs_fun,
  #                                  varv == unique_designs[j,]$varv,
  #                                  covuv == unique_designs[j,]$covuv)
  #   data_to_graph <- data_to_graph %>% mutate(stat = rep(statistics[i], b)) %>% mutate(x = c(1:b))
  #   
  #   da[[i]] <- data_to_graph
  # }
  # da <- reduce(da, rbind)
  # da$stat <- factor(da$stat, levels = c("Bias", "Variance", "MSE"))
  
  # That selects a design and a bandwidth
  unique_designs <- unique(data[, c(2:6)])    # As the same designs show up multiple times due to the different values of the bandwidth, I report the unique designs
  dat <- list()
  if( by == "sample")
  {
    # Power curves for each sample, one bandwidth
    for (i in 1:length(N))
    {
      datum <- dplyr::filter( data, sample == N[i],
                              varu == unique_designs[j ,]$varu,
                              fs_fun == unique_designs[j ,]$fs_fun,
                              varv == unique_designs[j,]$varv,
                              covuv == unique_designs[j,]$covuv,
                              bandwidth == bandrange[1] )  # We loop over sample sizes, not bandwidth in this chunk of code
      
      # dat_ols <- H0_range %>% map_dfr(rej, data = datum, est = datum$ols, se = datum$ols_se, level = level, name = "ols")
      # dat_tsls <- H0_range %>% map_dfr(rej, data = datum, est = datum$tsls, se = datum$tsls_se, level = level, name = "tsls")
      # dat_np <- H0_range %>% map_dfr(rej, data = datum, est = datum$np_endo, se = sqrt(datum$np_endo_var), level = level, name = "np")
      # dat_ll <- H0_range %>% map_dfr(rej, data = datum, est = datum$ll_endo, se = sqrt(datum$ll_endo_var), level = level, name = "ll")
      # dat_lc <- H0_range %>% map_dfr(rej, data = datum, est = datum$lc_endo, se = sqrt(datum$lc_endo_var), level = level, name = "lc")
      
      # If we look at the exogenous variables simulation
      dat_ols <- H0_range %>% map_dfr(rej, data = datum, est = datum$ols_x, se = sqrt(datum$ols_x_var), level = level, name = "ols")
      dat_tsls <- H0_range %>% map_dfr(rej, data = datum, est = datum$tsls_x, se = sqrt(datum$tsls_x_var), level = level, name = "tsls")
      dat_np <- H0_range %>% map_dfr(rej, data = datum, est = datum$np_x, se = sqrt(datum$np_x_var), level = level, name = "np")
      dat_ll <- H0_range %>% map_dfr(rej, data = datum, est = datum$ll_x, se = sqrt(datum$ll_x_var), level = level, name = "ll")
      dat_lc <- H0_range %>% map_dfr(rej, data = datum, est = datum$lc_x, se = sqrt(datum$lc_x_var), level = level, name = "lc")
      
      
      dat[[i]] <- dat_ols %>% left_join(dat_tsls, by = NULL) %>%
        left_join(dat_np, by = NULL) %>%
        left_join(dat_ll, by = NULL) %>%
        left_join(dat_lc, by = NULL) 
    }
    name_graph <- paste("power graph, ", unique_designs[j, ]$fs_fun,
                        # "n=", unique_designs[j, ]$sample,
                        ", varu=", unique_designs[j, ]$varu,
                        ", varv=", unique_designs[j, ]$varv,
                        ", covuv=", unique_designs[j, ]$covuv, 
                        ", bandwidth=", bandrange[1],
                        ", samples",
                        ".jpg", sep = "")
  }
  else if (by == "bandwidth")
  {
    # Power curves for each bandwidth, one sample size
    for (i in 1:length(bandrange))
    {
      datum <- dplyr::filter( data, sample == unique_designs[j,]$sample,
                              varu == unique_designs[j ,]$varu,
                              fs_fun == unique_designs[j ,]$fs_fun,
                              varv == unique_designs[j,]$varv,
                              covuv == unique_designs[j,]$covuv,
                              bandwidth == bandrange[i] )  # The bandwidth is fixed a well as the power graphs considered don't depend on the bandwidth
      
      # dat_ols <- H0_range %>% map_dfr(rej, data = datum, est = datum$ols, se = datum$ols_se, level = level, name = "ols")
      # dat_tsls <- H0_range %>% map_dfr(rej, data = datum, est = datum$tsls, se = datum$tsls_se, level = level, name = "tsls")
      # dat_np <- H0_range %>% map_dfr(rej, data = datum, est = datum$np_endo, se = sqrt(datum$np_endo_var), level = level, name = "np")
      # dat_ll <- H0_range %>% map_dfr(rej, data = datum, est = datum$ll_endo, se = sqrt(datum$ll_endo_var), level = level, name = "ll")
      # dat_lc <- H0_range %>% map_dfr(rej, data = datum, est = datum$lc_endo, se = sqrt(datum$lc_endo_var), level = level, name = "lc")
      # If we look at the exogenous variables simulation
      dat_ols <- H0_range %>% map_dfr(rej, data = datum, est = datum$ols_x, se = sqrt(datum$ols_x_var), level = level, name = "ols")
      dat_tsls <- H0_range %>% map_dfr(rej, data = datum, est = datum$tsls_x, se = sqrt(datum$tsls_x_var), level = level, name = "tsls")
      dat_np <- H0_range %>% map_dfr(rej, data = datum, est = datum$np_x, se = sqrt(datum$np_x_var), level = level, name = "np")
      dat_ll <- H0_range %>% map_dfr(rej, data = datum, est = datum$ll_x, se = sqrt(datum$ll_x_var), level = level, name = "ll")
      dat_lc <- H0_range %>% map_dfr(rej, data = datum, est = datum$lc_x, se = sqrt(datum$lc_x_var), level = level, name = "lc")
      
      
      
      dat[[i]] <- dat_ols %>% left_join(dat_tsls, by = NULL) %>%
        left_join(dat_np, by = NULL) %>%
        left_join(dat_ll, by = NULL) %>%
        left_join(dat_lc, by = NULL)
    }
    name_graph <- paste("power graph, ", unique_designs[j, ]$fs_fun,
                        ", n=", unique_designs[j, ]$sample,
                        ", varu=", unique_designs[j, ]$varu,
                        ", varv=", unique_designs[j, ]$varv,
                        ", covuv=", unique_designs[j, ]$covuv,
                        ", bandwidths",
                        ".jpg", sep = "")
  }
  dat <- reduce(dat, rbind)
  
  #unique_designs <- unique(dat[, 1:5])    # As the same designs show up multiple times due to the different values of the bandwidth, I report the unique designs
  
  #  data_to_graph <- dplyr::filter(dat, sample == unique_designs[j,]$sample, varu == unique_designs[j ,]$varu, fs_fun == unique_designs[j ,]$fs_fun,
  # varv == unique_designs[j,]$varv,
  # covuv == unique_designs[j,]$covuv)
  
  
  p <-  ggplot(data = dat, aes(x = H0)) +
    geom_line (aes( y = ols, color = "ols"), size = 1 ) +
    geom_line (aes( y = tsls, color = "tsls"), size = 1 ) +
    geom_line (aes( y = np, color = "np_power"), size = 1 ) +
    # geom_line (aes( y = ll, color = "ll"), size = 1 ) +
    # geom_line (aes( y = lc, color = "lc"), size = 1 ) +
    
    
    # xlab("H0") +
    xlab(TeX('$H_0$')) +
    ylab("Probability of rejection") +
    theme( axis.line = element_line(colour = "darkblue", 
                                    size = 1, linetype = "solid")) +
    # scale_x_discrete(breaks = as.character(1:b),
    #                  labels = as.character(bandrange) ) +
    # scale_x_discrete(limits = as.character(1:b), breaks = as.character(1:b) ,
    #                  labels =  c(".001", ".005", ".01", ".05", ".1",
    #                              ".5", "1", 
    #                              #"2", 
    #                              #"5",
    #                              "10",
    #                              #"20",
    #                              "100")
    # ) +
  scale_color_manual( breaks = c(  #"lc",
    # "ll",
    "np_power",
    "ols",
    "tsls" ),       # Put them in alphabetical order,it is easier for colours.
    labels = list( #TeX('$\\hat{\\beta}_{lc}$') ,
      #TeX('$\\hat{\\beta}_{ll}$') ,
      TeX('$\\hat{\\beta}_{NP}$'),
      TeX('$\\hat{\\beta}_{OLS}$') ,
      TeX('$\\hat{\\beta}_{2SLS}$')
    ),
    
    values = c(
      # "#663300", 
      #"#3300CC", 
      # "yellow",
      # "#00BFC4",   # Some kind of green-blue     # for beta_lc
      # "blue",          # for beta_lc
      # "orangered3",       # for beta_ll         
      # "orange",       # For for betahat_NP      
      "yellow2",       # For for betahat_NP      
      
      # "purple", 
      # "#00BFC4",   
      # "cyan",       # for OLS
      # "yellow2",    # 2SLS
      "olivedrab3",   # OLS
      "#E69F00"    # 2SLS
      # "#00BA38", # "lawngreen",  
      # # "#00BFC4", 
      # "#00BFC4" # "#3300CC" 
      # "#C77CFF"    # for WMDF
    )) +      
    labs(color = "Estimators") +
    scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1)) +
    theme(axis.text.x = element_text(angle = 45))
  
  if(by == "sample")  
    
  {    p <- p + facet_wrap(~ sample, scales = "free", labeller = as_labeller(c(`100` = "n = 100", `1000` = "n = 1000") ) ) } 
  else if (by == "bandwidth")    
  {    
    p <- p + facet_wrap(~ bandwidth, scales = "free", labeller = as_labeller(c(`0.001` = "h = 0.001", `0.005` = "h = 0.005", `0.01` = "h = 0.01", `0.05` = "h = 0.05", `0.1` = "h = 0.1", `0.5` = "h = 0.5", `1` = "h = 1", `10` = "h = 10", `100` = "h = 100") ) )    
  }     
  p <- p + theme(axis.title.x = element_text(family = "serif"),   # Changes fonts into times new roman
                 axis.title.y = element_text(family = "serif"),
                 legend.text = element_text(family = "serif"),
                 legend.title = element_text(family = "serif"),
                 strip.text.x = element_text(family = "serif")) +    
    theme(aspect.ratio=1)     # for the non-zoomed ones
  if(by == "sample")  
    
  {    ggsave(filename = name_graph, device = "jpg", height = 6, units = "cm" ) } 
  else if (by == "bandwidth")    
  {    
    ggsave(filename = name_graph, device = "jpg")
  }      
  # ggsave(filename = name_graph, device = "jpg", height = 6, units = "cm" )
  # ggsave(filename = name_graph, device = "jpg")
  
  # on.exit(dev.off())
  graphics.off()
  return("done !")
  
}