# Libraries

library(data.table)
library(Polychrome)

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

files <-
  list.files('../../Public',
             'slurm*.*.csv$',
             recursive = TRUE,
             full.names = TRUE)

# Filtering

files <- files[-grep('threecells', files)]

# Looping through files

mem.dt <- rbindlist(lapply(strsplit(files, '/'), function(x) {
  # Does not run for MACS2
  if (length(grep('MACS2', x)) == 0) {
    subx <- strsplit(x[7], '_|bp')[[1]]
    
    csvx <- fread(paste(x, collapse = '/'))[2,]$MaxRSS
    
    # Do not proceed if job is not completed
    if (!is.na(csvx)) {
      # If ChIPComp or DiffBind, need to check MACS2
      if (x[4] %in% c('ChIPComp', 'DiffBind')) {
        csvmacs2 <-
          rbindlist(lapply(files[intersect(grep(paste0('MACS2/', x[5]), files), grep(paste0('Helas3|Hepg2'), files))], function(y) {
            fread(y)[JobName == 'batch',]
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

# Organizing factors

mem.dt$BP %<>% factor(levels = c(250, 500, 750, 1000)) %<>% plyr::mapvalues(from = c(250, 500, 750, 1000), to = paste0(c(250, 500, 750, 1000), ' bp'))

mem.dt$Method %<>% plyr::mapvalues(
  from = c('ChIPComp', 'DiffBind'),
  to = c('ChIPComp +\nMACS2', 'DiffBind +\nMACS2')
)

# Plot

fig <-
  ggplot(data = mem.dt[Mark %in% c('CTCF', 'H3K27ac', 'H3K4me3'), ], aes(x = Method,
                                                                           y = Memory,
                                                                           fill = Method)) +
  facet_grid(BP ~ Mark) +
  geom_col() +
  guides(x = guide_axis(angle = 45), fill = guide_legend(nrow = 1, order = 1)) +
  theme_bw() +
  ylab('Memory (in Gb)') +
  scale_fill_manual(values = colors) +
  theme(legend.direction = 'horizontal', legend.position = 'top')

# Render
ggsave(
  plot = fig,
  filename = 'Figure_Memory_Time.pdf',
  height = 7,
  width = 8,
  dpi = 'retina'
)
