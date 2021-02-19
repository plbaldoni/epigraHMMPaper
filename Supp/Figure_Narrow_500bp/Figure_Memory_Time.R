# Libraries

library(data.table)
library(Polychrome)
library(magrittr)
library(ggplot2)
library(microbenchmark)

# Setup

colors = c(kelly.colors(22)[c('lightblue',
                              'red',
                              'yellow',
                              'purplishpink',
                              'yellowgreen',
                              'buff',
                              'orange')])
methods <-
  c(
    "epigraHMM",
    "ChIPComp +\nMACS2",
    "csaw",
    "DiffBind +\nMACS2",
    "diffReps",
    'RSEG',
    'THOR'
  )
names(colors) <- methods

################################################################################
### Memory
################################################################################

files.csv <-
  list.files('../../Public',
             'slurm*.*.csv$',
             recursive = TRUE,
             full.names = TRUE)

# Filtering

files.csv <- files.csv[grep('Encode_twocells|MACS2', files.csv)]

# Looping through files

mem.dt <- rbindlist(lapply(strsplit(files.csv, '/'), function(x) {
  # Does not run for MACS2
  if (length(grep('MACS2', x)) == 0) {
    subx <- strsplit(x[7], '_|bp')[[1]]
    
    csvx <- fread(paste(x, collapse = '/'))[2, ]$MaxRSS
    
    # Do not proceed if job is not completed
    if (!is.na(csvx)) {
      # If ChIPComp or DiffBind, need to check MACS2
      if (x[4] %in% c('ChIPComp', 'DiffBind')) {
        csvmacs2 <-
          rbindlist(lapply(files.csv[intersect(grep(paste0('MACS2/', x[5]), files.csv), grep(paste0('Helas3|Hepg2'), files.csv))], function(y) {
            fread(y)[JobName == 'batch', ]
          }))$MaxRSS
        
        mem <-
          max(as.numeric(strsplit(csvx, 'K')[[1]]), max(as.numeric(unlist(
            strsplit(as.character(csvmacs2), 'K')
          ))))
      } else{
        mem <- as.numeric(strsplit(csvx, 'K')[[1]])
      }
      
      return(data.table(
        Method = x[4],
        Mark = x[5],
        Data = x[6],
        BP = subx[5],
        Memory = mem * 1e-6
      )) # From Kb to Gb
    }
  }
}))

################################################################################
### Time
################################################################################

files.rda <-
  list.files('../../Public',
             '*Time*.*.RData',
             recursive = TRUE,
             full.names = TRUE)

files.rda <- files.rda[grep('Encode_twocells|MACS2', files.rda)]

# Looping

time.dt <- rbindlist(lapply(strsplit(files.rda, '/'), function(x) {
  # Does not run for MACS2
  if (length(grep('MACS2', x)) == 0) {
    subx <- strsplit(tail(x, 1), '_|bp')[[1]]
    
    load(paste(x, collapse = '/'))
    time.method <-
      summary(cptime[[1]], 's')$mean / 60 # Time in minutes
    rm(cptime)
    
    # If ChIPComp or DiffBind, need to check MACS2
    if (x[4] %in% c('ChIPComp', 'DiffBind')) {
      time.macs <- c(0, 0)
      files.macs <-
        files.rda[intersect(grep(paste0('MACS2/', x[5]), files.rda), grep(paste0('Helas3|Hepg2'), files.rda))]
      for (i in 1:2) {
        load(files.macs[i])
        time.macs[i] <-
          summary(cptime[[1]], 's')$mean / 60 # Time in minutes
        rm(cptime)
      }
      
      time.total <- time.method + sum(time.macs)
    } else{
      time.total <- time.method
    }
    
    return(data.table(
      Method = x[4],
      Mark = x[5],
      Data = x[6],
      BP = subx[6],
      Time = time.total
    ))
  }
}))

################################################################################
### Plotting
################################################################################

dt <-
  merge(mem.dt,
        time.dt,
        by = c('Method', 'Mark', 'Data', 'BP'),
        all.x = TRUE)

# Organizing factors

dt$BP %<>% factor(levels = c(250, 500, 750, 1000)) %<>% plyr::mapvalues(from = c(250, 500, 750, 1000), to = paste0(c(250, 500, 750, 1000), ' bp'))

dt$Method %<>% plyr::mapvalues(
  from = c('ChIPComp', 'DiffBind'),
  to = c('ChIPComp +\nMACS2', 'DiffBind +\nMACS2')
)

# Plot

dat_text <- data.frame(
  Mark = dt[Mark %in% c('CTCF', 'H3K27ac', 'H3K4me3') & BP == '500 bp' & Method == 'RSEG',Mark],
  label = paste0('RSEG (Not Shown):\n',
                 paste0(dt[Mark %in% c('CTCF', 'H3K27ac', 'H3K4me3') & BP == '250 bp' & Method == 'RSEG',formatC(Memory,digits = 2,format = 'f')],' GB\n'),
                 paste0(dt[Mark %in% c('CTCF', 'H3K27ac', 'H3K4me3') & BP == '250 bp' & Method == 'RSEG',formatC(Time/60,digits = 2,format = 'f')],' hours'),sep = ' '),
  x     = rep(1.75,3),
  y     = rep(10,3)
)

fig_memory_time <-
  ggplot(data = dt[Mark %in% c('CTCF', 'H3K27ac', 'H3K4me3') & BP == '500 bp',], aes(x = Time/60,
                                                                                     y = Memory,
                                                                                     color = Method)) +
  facet_grid(rows = 'Mark') +
  geom_point() +
  guides(color = guide_legend(nrow = 1, order = 1)) +
  theme_bw() +
  ylab('Memory (in GB)') + xlab('Time (in hours)') +
  scale_color_manual(values = colors) +
  theme(legend.direction = 'horizontal', legend.position = 'top')

save(fig_memory_time,file = './Figure_Memory_Time.RData')