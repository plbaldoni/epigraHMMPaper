# Table: simulation results

library(data.table)
library(ggplot2)
library(magrittr)
library(plyr)
library(ggpubr)
library(RColorBrewer)
library(kableExtra)

### Simulation scnearios
# Marker:   H3K27me3    H3K36me3
# Groups:   2, 3, 4
# Replicates: 1, 2, 4
# Nwindow: 1e5, 5e5, 1e6
# Nsim: 100

scenario = expand.grid(
  Marker = c('H3K27me3', 'H3K36me3'),
  Groups = c(2, 3, 4),
  Replicates = c(1, 2, 4),
  Nwindow = c(1e5, 5e5, 1e6),
  Nsim = 100,
  Pct = round(100 * c(0.7, 0.8, 0.9, 1))
)
scenario$Label = with(scenario,
                      paste(
                        Marker,
                        paste0(Groups, 'G'),
                        paste0(Replicates, 'R'),
                        paste0(Nwindow, 'W'),
                        paste0(Pct, 'SNR'),
                        sep = "_"
                      ))
scenario <-
  scenario[order(scenario$Groups,
                 scenario$Replicates,
                 scenario$Nwindow,
                 -scenario$Pct), ]
# print(scenario)

### True parameter functions

genpar = function(Marker, Groups, pct) {
  ## Initial probability (always begin with background) and Transition Matrix
  pi = c(1, 0)
  ## Number of states
  N.states = 2 ^ Groups
  if (Marker == 'H3K36me3') {
    ## See 01_H3K36me3_Parameters.R
    ## Transition Probability Matrix
    d1 = 1 - 1 / (131.92967 + 1)
    d2 = 1 - 1 / (43.75663 + 1)
    gamma = matrix(0, N.states, N.states)
    diag(gamma) = c(d1, rep(d2, N.states - 1))
    gamma[1, 2:N.states] = (1 - d1) / (N.states - 1)
    gamma[gamma == 0] = (1 - d2) / (N.states - 1)
    
    ## ChIP parameters
    par.background = c(0.7560612, 0, 3.230474)
    par.enrichment = c(2.4169519, 0, 2.891908)
    
    ## Offsets
    #I am using the same offset for replicates from the same group
    offsets = NULL
    offsets[paste0('Group', 1:Groups)] = list(NULL)
    offsets[paste0('Group', 1:Groups)] = rep(0, Groups)
  }
  if (Marker == 'H3K27me3') {
    ## See 02_H3K27me3_Parameters.R
    ## Transition Probability Matrix
    d1 = 1 - 1 / (143.57818 + 1)
    d2 = 1 - 1 / (49.47422 + 1)
    gamma = matrix(0, N.states, N.states)
    diag(gamma) = c(d1, rep(d2, N.states - 1))
    gamma[1, 2:N.states] = (1 - d1) / (N.states - 1)
    gamma[gamma == 0] = (1 - d2) / (N.states - 1)
    
    ## ChIP parameters
    par.background = c(1.11594274, 0, 3.599259)
    par.enrichment = c(2.2810865, 0, 4.076325)
    
    ## Offsets
    #I am using the same offset for replicates from the same group
    offsets = NULL
    offsets[paste0('Group', 1:Groups)] = list(NULL)
    offsets[paste0('Group', 1:Groups)] = rep(0, Groups)
  }
  # Adjusting for SNR pct
  snr = exp(par.enrichment[1]) / exp(par.background[1])
  var.enrichment = exp(par.enrichment[1]) * (1 + exp(par.enrichment[1]) /
                                               par.enrichment[3])
  meanvar.enrichment = exp(par.enrichment[1]) / (var.enrichment)
  par.enrichment.new = c(log(pct * snr * exp(par.background[1])),
                         0,
                         (pct * snr * exp(par.background[1])) * meanvar.enrichment / (1 - meanvar.enrichment))
  par.enrichment = par.enrichment.new
  
  return(
    list(
      'pi' = pi,
      'gamma' = gamma,
      'par.background' = par.background,
      'par.enrichment' = par.enrichment,
      'offsets' = offsets
    )
  )
}

### Loading Predicted Pattern

pct = c(70, 80, 90, 100)

par.lst = list()
for (i in pct) {
  outRData <-
    paste0('../../Sim/ReadCounts/Metrics_',
           i,
           'SNR/',
           'Metrics_',
           i,
           'SNR_par.RData')
  if (file.exists(outRData)) {
    load(outRData)
    par.lst[[paste0(i)]] = par.list
  }
}
par.lst <- unlist(par.lst, recursive = F)
names(par.lst) <- gsub('.*\\.', '', names(par.lst))

par <- lapply(
  par.lst,
  FUN = function(x) {
    ### Making each sublist a data.table
    newx = lapply(
      x,
      FUN = function(y) {
        return(as.data.table(y, keep.rownames = "Parameter"))
      }
    )
    
    newx = rbindlist(newx, idcol = 'Simulation')
    setnames(newx, c('Simulation', 'Parameter', 'Estimate'))
    return(newx)
  }
)
par = rbindlist(par, idcol = 'Label')

label <- as.data.table(do.call(rbind, strsplit(par$Label, "_")))
setnames(label, c('Marker', 'Groups', 'Replicates', 'Nwindow', 'Pct'))
label[, Groups := as.numeric(gsub("G", "", Groups))]
label[, Replicates := as.numeric(gsub("R", "", Replicates))]
label[, Nwindow := as.numeric(gsub("W", "", Nwindow))]
label[, Pct := as.numeric(gsub("S.*", "", Pct))]

par <- data.table(par[, Label := NULL], label)

### Bringing the true parameters

par$True <- apply(
  par,
  1,
  FUN = function(x) {
    truepar = genpar(Marker = x["Marker"],
                     Groups = as.numeric(x['Groups']),
                     pct = as.numeric(x['Pct']) / 100)
    
    if (x['Parameter'] == "HMM.Mean.Int") {
      return((truepar$par.background[1]))
    }
    if (x['Parameter'] == "HMM.Mean.Diff") {
      return((truepar$par.enrichment[1] - truepar$par.background[1]))
    }
    if (x['Parameter'] == "HMM.Disp.Int") {
      return(log(truepar$par.background[3]))
    }
    if (x['Parameter'] == "HMM.Disp.Diff") {
      return(log(truepar$par.enrichment[3] / truepar$par.background[3]))
    }
  }
)

### Tabulating the results

sumpar <- par[, RBias := (Estimate - True) / True][, .(
  Mean = mean(RBias),
  LB = quantile(RBias, probs = 0.025),
  UB = quantile(RBias, probs = 0.975)
), by = c('Marker',
          'Groups',
          'Replicates',
          'Nwindow',
          'Pct',
          'Parameter',
          'True')]
sumpar$Parameter %<>% as.factor() %<>%  factor(levels = c(
  'HMM.Mean.Int',
  'HMM.Mean.Diff',
  'HMM.Disp.Int',
  'HMM.Disp.Diff'
))

sumpar.subset <-
  sumpar[Marker == 'H3K27me3' &
           Nwindow == 1e6 &
           Pct == 100, -c('Marker', 'Nwindow', "Pct")]
sumpar.subset <-
  data.table::dcast(sumpar.subset,
                    Groups + Parameter + True ~ Replicates,
                    value.var = c("Mean", "LB", "UB"))

sumpar.subset[, CI_1 := paste0('(', sprintf("%.3f", round(LB_1, 3)), ', ', sprintf("%.3f", round(UB_1, 3)), ')')]
sumpar.subset[, CI_2 := paste0('(', sprintf("%.3f", round(LB_2, 3)), ', ', sprintf("%.3f", round(UB_2, 3)), ')')]
sumpar.subset[, CI_4 := paste0('(', sprintf("%.3f", round(LB_4, 3)), ', ', sprintf("%.3f", round(UB_4, 3)), ')')]
sumpar.subset <-
  sumpar.subset[, -c('LB_1', 'LB_2', 'LB_4', 'UB_1', 'UB_2', 'UB_4')]

setcolorder(sumpar.subset,
            c(
              'Groups',
              'Parameter',
              'True',
              paste0(c('Mean_', 'CI_'), 1),
              paste0(c('Mean_', 'CI_'), 2),
              paste0(c('Mean_', 'CI_'), 4)
            ))

sumpar.subset$Parameter %<>% mapvalues(
  from = c(
    'HMM.Mean.Int',
    'HMM.Mean.Diff',
    'HMM.Disp.Int',
    'HMM.Disp.Diff'
  ),
  to = c(
    '$\\beta_{1}$',
    '$\\beta_{3}$',
    '$\\lambda_{1}$',
    '$\\lambda_{3}$'
  )
)

sumpar.subset$Groups %<>% mapvalues(from = 2:4, to = c('Two', 'Three', 'Four'))

### Tabulating

sink("Table_Simulation_ReadCounts.tex")
sumpar.subset %>%
  dplyr::mutate_all(linebreak) %>%
  kable(
    format = 'latex',
    digits = 3,
    col.names = c(
      'Conditions',
      'Parameter',
      'True Value',
      'Relative Bias',
      '($P_{2.5}$, $P_{97.5}$)',
      'Relative Bias',
      '($P_{2.5}$, $P_{97.5}$)',
      'Relative Bias',
      '($P_{2.5}$, $P_{97.5}$)'
    ),
    align = c('l', 'l', 'r', 'r', 'r', 'r', 'r', 'r', 'r'),
    booktabs = T,
    escape = F,
    caption = 'Read count simulation. True values and average relative bias of parameter estimates (and $2.5^{th}$, $97.5^{th}$ percentiles) across a hundred simulated data sets are shown for H3K27me3 data with $10^6$ genomic windows and ENCODE-estimated SNR.\\label{Table:table_simulation}'
  ) %>%
  add_header_above(c(
    ' ' = 3,
    'One Replicate' = 2,
    'Two Replicates' = 2,
    'Four Replicates' = 2
  )) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  collapse_rows(columns = 1, latex_hline = 'major')
sink()