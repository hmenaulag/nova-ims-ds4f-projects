#-------------------------------------------------------------------------------------------------
# NOVA IMS Executive Education - DATA SCIENCE FOR FINANCE
# Asset Pricing & Portfolio Management
# Group Project - Elements:
#                 1 - Amina Baklizi – 20230515@novaims.unl.pt
#                 2 - Malik Harrak – 20231140@novaims.unl.pt
#                 3 - Hugo Laginha – 20231130@novaims.unl.pt
#                 4 - Saad Islam – 20230513@novaims.unl.pt
# Track no. 2 - Backtesting Risk-Based Portfolios
#-------------------------------------------------------------------------------------------------

rm(list=ls(all.names = TRUE))
graphics.off()
close.screen(all.screens = TRUE)
erase.screen()
#windows.options(record=TRUE)
options(scipen=999)

#WINDOWS Users - UPDATE FOR PC
# ipath <- 'C:/Users/Jorge Bravo/Desktop/Teaching/Asset Pricing Portfolio Theory'
#UBUNTU user
ipath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(ipath)

library(pacman)

p_load(xts)                   # to manipulate time series of stock data
p_load(quantmod)              # to download stock data
p_load(PerformanceAnalytics)  # to compute performance measures
p_load(portfolioBacktest)     # to perform backtest
p_load(CVXR)                  # Convex Optimization in R
p_load(riskParityPortfolio)   # RPP
p_load(HierPortfolios)        # HRP
p_load(pracma)
p_load(DT)
p_load(ggplot2)               # plot
p_load(dplyr)
p_load(fs)

#------------------------------------------------------------------------------------------------
# Overview: 
# The project consists of empirically testing the performance of selected portfolio risk based
# investment strategies.
#------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------
# Task 1: 
# Download monthly market data for S&P500 listed stocks from 01-01-2010 to 31-12-2022.
#------------------------------------------------------------------------------------------------

# Download of Yahoo Finance tickers
data(SP500_symbols)

# Definition of time frame
start_date <- "2010-01-01"
end_date <- "2022-12-31"

# Download historical data from Yahoo Finance
SP500 <- stockDataDownload(
  stock_symbols = SP500_symbols,
  index_symbol = '^GSPC',
  from = start_date,
  to = end_date,
  periodicity = "monthly",  # Corrected argument name
  rm_stocks_with_na = TRUE
)

save(SP500, file = "stockdata_from_2010-01-01_to_2022-12-31_monthly.Rdata")
load("stockdata_from_2010-01-01_to_2022-12-31_monthly.Rdata")
#------------------------------------------------------------------------------------------------
# Task 2: 
# Generate 100 random resamples of 20 S&P500 listed stocks and 3 consecutive years.
#------------------------------------------------------------------------------------------------
set.seed(11111)
samples <- financialDataResample(SP500,
                                 N_sample = 20,
                                 T_sample = 12 * 3,
                                 num_datasets = 100,
                                 rm_stocks_with_na = TRUE
                                 )
View(samples)

#------------------------------------------------------------------------------------------------
# Task 3: 
# Empirically investigate the performance of the following traditional and risk-based portfolios.
# &
# Task 4: 
# Critically discuss the results considering alternative risk-adjusted performance metrics (e.g.,
# returns, volatility, Sharpe ratio, Sterling ratio, drawdown).
#------------------------------------------------------------------------------------------------

# A. Global Minimum-Variance Portfolio with no short-selling
GMVP_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(X)
  SigInv <- pinv(Sigma)
  Ones <- rep(1, nrow(Sigma))
  # design GMVP
  w <- SigInv%*%Ones /as.numeric(t(Ones)%*%SigInv%*% Ones)
  w <- abs(w) / sum(abs(w))
  names(w) <- colnames(Sigma)
  return(w)
}

# B. Inverse Volatility Portfolio
IVP_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(X)  # compute SCM
  # design IVP
  w <- riskParityPortfolio(Sigma, formulation='diag')$w
  return(w)
}
# C. Risk parity portfolio
RPP_fun <- function(dataset, ...){
  X <- diff(log(dataset$adjusted))[-1]
  Sigma <- cov(X)
  rpp <- riskParityPortfolio(Sigma)
  return(rpp$w)
}
# D. Most diversified portfolio
MDP_fun <- function(dataset,...) {
  X <- diff(log(dataset$adjusted))[-1]
  Sigma <- cov(X)
  mu = sqrt(diag(Sigma))
  w_ <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(quad_form(w_, Sigma)),
                  constraints = list(w_ >= 0, t(mu) %*% w_ == 1))
  result <- CVXR::solve(prob)
  w <- as.vector(abs(result$getValue(w_))/sum(abs(result$getValue(w_))))
  names(w) <- colnames(Sigma)
  return(w)
}
# E. Maximum decorrelation portfolio
MDCP_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]
  Sigma <- cov(X)
  C <- diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma)))
  colnames(C) <- colnames(Sigma)
  w <- Variable(nrow(C))
  prob <- Problem(Minimize(quad_form(w, C)),
                  constraints = list(w >= 0, sum(w) == 1))
  result <- CVXR::solve(prob)
  w_ <- w_ <- as.vector(abs(result$getValue(w))/sum(abs(result$getValue(w))))
  names(w_) <- colnames(C)
  return(w_)
}

# F. Hierarchical Risk Parity Portfolio
HRPP_ward <- function(dataset, ...){
  X <- diff(log(dataset$adjusted))[-1]
  Sigma <- cov(X)
  hrppward = HRP_Portfolio(as.matrix(Sigma), linkage = "ward", graph = T)
  w <- hrppward$weights
  names(w) <- colnames(Sigma)
  return(w)
}

# G. Markowitz’s mean-variance portfolio (MVP) with no short-selling
MMVP_fun <- function(dataset, lambda=0.5, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  mu    <- colMeans(X)  # compute mean vector
  Sigma <- cov(X)       # compute the SCM
  # design mean-variance portfolio
  w <- Variable(nrow(Sigma))
  prob <- Problem(Maximize(t(mu) %*% w - lambda*quad_form(w, Sigma)),
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  w_ <- as.vector(abs(result$getValue(w))/sum(abs(result$getValue(w))))
  names(w_) <- colnames(Sigma)
  return(w_)
}

# portfolio list
portfolios <- list("GMVP"   = GMVP_fun,
                   "IVP"    = IVP_fun,
                   "RPP"    = RPP_fun,
                   "MDP"    = MDP_fun,
                   "MDCP"   = MDCP_fun,
                   "HRPP"   = HRPP_ward,
                   "MMVP"   = MMVP_fun
)

# backtest (optimize every quarter, rebalance every month)
bt <- portfolioBacktest(portfolios, samples, 
                        benchmarks = c("1/N", "index"),   # benchmark portfolios
                        rebalance_every = 1,
                        optimize_every = 3,
                        lookback = 12,         # Length of the lookback rolling window in periods
                        shortselling = FALSE,
                        cost = list(buy = 0e-4, sell = 0e-4, short = 0e-4, long_leverage = 0e-4),
                        show_progress_bar = TRUE)

names(bt)

#Portfolios Weights (for all dates)
bt$GMVP$`dataset 1`$w_bop
bt$IVP$`dataset 1`$w_bop
bt$RPP$`dataset 1`$w_bop
bt$MDP$`dataset 1`$w_bop
bt$MDCP$`dataset 1`$w_bop
bt$HRPP$`dataset 1`$w_bop
bt$MMVP$`dataset 1`$w_bop
bt$`1/N`$`dataset 1`$w_bop
bt$index$`dataset 1`$w_bop


# The sum of the weights should be equal to 1 (invest all of the money)
sum(bt$GMVP$`dataset 1`$w_bop[1,])
sum(bt$IVP$`dataset 1`$w_bop[1,])
sum(bt$RPP$`dataset 1`$w_bop[1,])
sum(bt$MDP$`dataset 1`$w_bop[1,])
sum(bt$MDCP$`dataset 1`$w_bop[1,])
sum(bt$HRPP$`dataset 1`$w_bop[1,])
sum(bt$MMVP$`dataset 1`$w_bop[1,])
sum(bt$`1/N`$`dataset 1`$w_bop[1,])
sum(bt$index$`dataset 1`$w_bop[1,])

# Select several performance measures for each portfolio (dataset 1):
# Performance Criteria
bt$GMVP$`dataset 1`$performance
bt$IVP$`dataset 1`$performance
bt$RPP$`dataset 1`$performance
bt$MDP$`dataset 1`$performance
bt$MDCP$`dataset 1`$performance
bt$HRPP$`dataset 1`$performance
bt$MMVP$`dataset 1`$performance
bt$`1/N`$`dataset 1`$performance
bt$index$`dataset 1`$performance

# Performance Metrics per portfolio
backtestSelector(bt, portfolio_name = "GMVP", measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))
backtestSelector(bt, portfolio_name = "IVP", measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))
backtestSelector(bt, portfolio_name = "RPP", measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))
backtestSelector(bt, portfolio_name = "MDP", measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))
backtestSelector(bt, portfolio_name = "MDCP", measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))
backtestSelector(bt, portfolio_name = "HRPP", measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))
backtestSelector(bt, portfolio_name = "MMVP", measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))
backtestSelector(bt, portfolio_name = "1/N", measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))
backtestSelector(bt, portfolio_name = "index", measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))

# Summary
res_sum <- backtestSummary(bt)
names(res_sum)
res_sum$performance_summary
summaryTable(res_sum, type = "DT", order_col = "Sharpe ratio", order_dir = "desc")
# Graphic representation of Summary
summaryBarPlot(res_sum, measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))

# Box Plots
backtestBoxPlot(bt, measure = "Sharpe ratio")
backtestBoxPlot(bt, measure = "Sterling ratio")
backtestBoxPlot(bt, measure = "max drawdown")
backtestBoxPlot(bt, measure = "annual return")
backtestBoxPlot(bt, measure = "annual volatility")

# Cumulative Returns, Drawdown, Sharpe Ratio & Weight Allocation graphics:
chartCumRet_list <- list()
chartDrawdown_list <- list()
chartSharpe_list <- list()

for (i in 1:length(samples)) {
  chartCumRet_list[[i]] <- backtestChartCumReturn(bt, portfolios = names(bt), dataset_num = i)
  chartDrawdown_list[[i]] <- backtestChartDrawdown(bt, portfolios = names(bt), dataset_num = i)
  chartSharpe_list[[i]] <- backtestChartSharpeRatio(bt, portfolios = names(bt), dataset_num = i, lookback = 12, by = 1, gap = 1, bars_per_year = 12)
}
chartCumRet_list[[1]]
chartDrawdown_list[[1]]
chartSharpe_list [[1]]

backtestChartStackedBar(bt, "GMVP", dataset_num = 1, type = "ggplot2", legend = TRUE)
backtestChartStackedBar(bt, "IVP", dataset_num = 1, type = "ggplot2", legend = TRUE)
backtestChartStackedBar(bt, "RPP", dataset_num = 1, type = "ggplot2", legend = TRUE)
backtestChartStackedBar(bt, "MDP", dataset_num = 1, type = "ggplot2", legend = TRUE)
backtestChartStackedBar(bt, "MDCP", dataset_num = 1, type = "ggplot2", legend = TRUE)
backtestChartStackedBar(bt, "HRPP", dataset_num = 1, type = "ggplot2", legend = TRUE)
backtestChartStackedBar(bt, "MMVP", dataset_num = 1, type = "ggplot2", legend = TRUE)
backtestChartStackedBar(bt, "1/N", dataset_num = 1, type = "ggplot2", legend = TRUE)


# Relative Risk Contribution and weight allocation
barplotPortfolioRisk(rep(0.05,20), cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  ggtitle('Naive Portfolio - 1/N') +
  scale_y_continuous(labels = scales::percent)
barplotPortfolioRisk(GMVP_fun(samples$`dataset 1`), cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  ggtitle('Global Minimum Variance Portfolio (GMVP)') +
  scale_y_continuous(labels = scales::percent)
barplotPortfolioRisk(IVP_fun(samples$`dataset 1`), cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  ggtitle('Inverse volatility portfolio (IVP)') +
  scale_y_continuous(labels = scales::percent)
barplotPortfolioRisk(RPP_fun(samples$`dataset 1`), cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  ggtitle('Risk Parity portfolio (RPP)') +
  scale_y_continuous(labels = scales::percent)
barplotPortfolioRisk(MDP_fun(samples$`dataset 1`), cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  ggtitle('Most Diversified Portfolio (MDP)') +
  scale_y_continuous(labels = scales::percent)
barplotPortfolioRisk(MDCP_fun(samples$`dataset 1`), cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  ggtitle('Maximum Decorrelation Portfolio (MDCP)') +
  scale_y_continuous(labels = scales::percent)
barplotPortfolioRisk(HRPP_ward(samples$`dataset 1`), cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  ggtitle('Hierarchical Risk Parity Portfolio (HRPP ward)') +
  scale_y_continuous(labels = scales::percent)
barplotPortfolioRisk(MMVP_fun(samples$`dataset 1`), cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  ggtitle('Markowitz Mean Variance Portfolio (MMVP)') +
  scale_y_continuous(labels = scales::percent)

# Relative Risk Contribution and weight allocation (grouped by stocks)
w_all <- cbind("GMVP" = GMVP_fun(samples$`dataset 1`),
               "IVP"  = IVP_fun(samples$`dataset 1`),
               "RPP"  = RPP_fun(samples$`dataset 1`),
               "MDP" = MDP_fun(samples$`dataset 1`),
               "MDCP" = MDCP_fun(samples$`dataset 1`),
               "HRPP ward" = HRPP_ward(samples$`dataset 1`),
               "MMVP" = MMVP_fun(samples$`dataset 1`))

barplotPortfolioRisk(w_all, cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  scale_y_continuous(labels = scales::percent)
#-----------------------------------------------------------------------------------------------------
# Task 5. 
# Use the same data set to investigate the performance of portfolio management strategies using 
# alternative risk measures, e.g., downside risk, value-at-risk (VaR), Conditional VaR (CVaR) or 
# expected shortfall (ES).
#----------------------------------------------------------------------------------------------------
MVaR95P_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  mu    <- colMeans(X)  
  Sigma <- cov(X)
  n <- ncol(Sigma)
  epsilon <- 1e-6  # Small positive constant
  Sigma_reg <- Sigma + epsilon * diag(n)
  u = chol(Sigma_reg)
  
  w <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(1.645 * norm(u %*% w, "2")),
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  w_ <- as.vector(abs(result$getValue(w))/sum(abs(result$getValue(w))))
  return(w_)
}

portfolio_MVaR <- list("MVaR95" = MVaR95P_fun)
bt_MVAR <- portfolioBacktest(portfolio_MVaR, samples, 
                        rebalance_every = 1,
                        optimize_every = 3,
                        lookback = 12,         # Length of the lookback rolling window in periods
                        shortselling = FALSE,
                        cost = list(buy = 0e-4, sell = 0e-4, short = 0e-4, long_leverage = 0e-4),
                        show_progress_bar = TRUE)
bt <- append(bt, bt_MVAR)
names(bt)

#Portfolios Weights (for all dates)
bt$MVaR95$`dataset 1`$w_bop

# The sum of the weights should be equal to 1 (invest all of the money)
sum(bt$MVaR95$`dataset 1`$w_bop[1,])

# Select several performance measures for each portfolio (dataset 1):
# Performance Criteria
bt$MVaR95$`dataset 1`$performance

# Performance Metrics per portfolio
backtestSelector(bt, portfolio_name = "MVaR95", measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))

# Summary
res_sum <- backtestSummary(bt)
names(res_sum)
res_sum$performance_summary
summaryTable(res_sum, type = "DT", order_col = "Sharpe ratio", order_dir = "desc")

# Graphic representation of summary
summaryBarPlot(res_sum, measures = c("Sharpe ratio", "Sterling ratio", "max drawdown", "annual return", "annual volatility"))

# Box Plots
backtestBoxPlot(bt, measure = "Sharpe ratio")
backtestBoxPlot(bt, measure = "Sterling ratio")
backtestBoxPlot(bt, measure = "max drawdown")
backtestBoxPlot(bt, measure = "annual return")
backtestBoxPlot(bt, measure = "annual volatility")


# Cumulative Returns, Drawdown, Sharpe Ratio & Weight Allocation graphics:

for (i in 1:length(samples)) {
  chartCumRet_list[[i]] <- backtestChartCumReturn(bt, portfolios = names(bt), dataset_num = i)
  chartDrawdown_list[[i]] <- backtestChartDrawdown(bt, portfolios = names(bt), dataset_num = i)
  chartSharpe_list[[i]] <- backtestChartSharpeRatio(bt, portfolios = names(bt), dataset_num = i, lookback = 12, by = 1, gap = 1, bars_per_year = 12)
}
chartCumRet_list[[1]]
chartDrawdown_list[[1]]
chartSharpe_list [[1]]

backtestChartStackedBar(bt, "MVaR95", dataset_num = 1, type = "simple", legend = TRUE)


# Relative Risk Contribution and weight allocation
barplotPortfolioRisk(MVaR95P_fun(samples$`dataset 1`), cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  ggtitle('Risk Decomposition using Var (MVaR95P)') +
  scale_y_continuous(labels = scales::percent)

# Relative Risk Contribution and weight allocation (grouped by stocks)
w_all <- cbind("GMVP" = GMVP_fun(samples$`dataset 1`),
               "IVP"  = IVP_fun(samples$`dataset 1`),
               "RPP"  = RPP_fun(samples$`dataset 1`),
               "MDP" = MDP_fun(samples$`dataset 1`),
               "MDCP" = MDCP_fun(samples$`dataset 1`),
               "HRPP ward" = HRPP_ward(samples$`dataset 1`),
               "MMVP" = MMVP_fun(samples$`dataset 1`),
               "MVaR95" = MVaR95P_fun(samples$`dataset 1`))

barplotPortfolioRisk(w_all, cov(diff(log(samples$`dataset 1`$adjusted))[-1])) +
  scale_y_continuous(labels = scales::percent)
#---------------------------------------------------------------------------------
# EXTRAS
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Transaction Costs
#---------------------------------------------------------------------------------
install.packages("ggfortify")
library(ggfortify)
# backtest with costs of 60 bps
bt_tc <- portfolioBacktest(GMVP_fun, samples,
                           rebalance_every = 1,
                           optimize_every = 3,
                           lookback = 12,         # Length of the lookback rolling window in periods
                           shortselling = FALSE,
                           cost = list(buy = 60e-4, sell = 60e-4),
                           show_progress_bar = TRUE)

# plot wealth time series
wealth <- cbind(bt$GMVP$`dataset 1`$wealth, bt_tc$fun1$`dataset 1`$wealth)
colnames(wealth) <- c("without transaction costs", "with transaction costs")

autoplot(wealth, facets = FALSE, main = "Wealth") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.8, 0.2)) +
  scale_color_manual(values = c("red", "black"))