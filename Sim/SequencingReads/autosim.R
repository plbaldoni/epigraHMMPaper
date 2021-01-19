# It needs is.tf and iters arguments here

# This simulation pipeline was adapted from csaw (Nucleic Acids Research 44,
# e45 (2015), http://bioinf.wehi.edu.au/csaw/).
# We adapted the codes that were made available by the authors in order to obtain
# broad peaks with diffused ChIP-seq signals. We thank the authors Aaron Lun and
# Gordon K. Smyth for providing their entire pipeline available and reproducible.
# Please, refer to their repository for the original codes and to download to the
# package 'xscss'.

# This simulation generates a mixture of complex peak structures.
# This is done by treating each peak as a mixture of three subpeaks.
# For DB peaks, one or more of those subpeaks are randomly chosen and eliminated.

# Loading libraries

library(epigraHMM)
library(SummarizedExperiment)
library(DiffBind)
library(data.table)
library(ChIPComp)
library(xscss)
library(magrittr)
library(csaw)
library(edgeR)
library(microbenchmark)
library(rhdf5)

# General parameters

fraglen <- 100
npeaks <- 4000
nde <- 1000

prior.df <- 20
dispersion <- 1 / prior.df
grouping <- c("A", "A", "B", "B")
true.width <- 10000
base.mu <- 300
dist <- c(60000, 65000)
back.mu <- c(9.99, 10.01)
back.width <- 2000
back.rlen <- 10
dist.peak <- c(0.9, 1.8)
convtime <- 6e10 # Converts nanoseconds to minutes

if (is.tf) {
    # We do not simulate TF data, we kept this for completeness
    all.fix <- paste0("tf", iters)
    radius <- fraglen
    if (is.homo) {
        prior.df <- 1e8
        all.fix <- paste0("tfx", iters)
        base.mu <- 10
        up.mu <- 20
        down.mu <- 0
    } else {
        up.mu <- 45
        down.mu <- 15
    }
    binsize = 100
} else {
    is.homo <- FALSE
    all.fix <- paste0("hist", iters)
    radius <- true.width / 2L
    binsize <- 200
}

design <- model.matrix( ~ factor(grouping))
ofile <- paste0(all.fix, "_out.tsv")
fdr.thresholds <- c(0.01, 0.05, 0.1, 0.15, 0.2)

if (iters == 1) {
    results = './results/'
    system(paste('mkdir', results))
}

result.file <- paste0(all.fix, "_result.tsv")
unlink(result.file)
dump2file <- function(id, cutoff, result) {
    # This function saves the results output from assessChIP2
    write.table(
        file = result.file,
        data.frame(
            Method = id,
            Cutoff = cutoff,
            OneMinusPrecision = 1 - result$overlap /
                result$found,
            Recall = result$recall,
            Time = result$time,
            CallSize = result$callsize,
            PeakSize = result$peaksize,
            HitsPerPeak = result$hitsperpeak,
            NCalls = result$ncalls,
            NPeaks = result$npeaks,
            TPR = result$tpr,
            FPR = result$fpr,
            PPV = result$ppv,
            Error = result$error
        ),
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        append = TRUE
    )
}

assessChIP2 = function(observed,
                       known,
                       time,
                       tol = 200,
                       checkfc = TRUE,
                       gsize = gsize,
                       is.error = F)
{
    # This function calculates the metrics (TPR, FPR, FDR, etc.) given a set
    # of filtered peak ranges (observed) and known true differential peak locations
    
    if (is.error == F) {
        gr.gensize = GRanges('chrA', IRanges(1, gsize))
        obs <- read.table(observed, header = TRUE)
        if (checkfc) {
            up.o <- obs$logFC > 0
            if (is.null(up.o)) {
                stop("need a log-FC field in the table of observed sites")
            }
        } else {
            up.o <- rep(TRUE, nrow(obs))
        }
        obranges <- GRanges(obs$chr, IRanges(obs$start, obs$end))
        if (!is.na(tol)) {
            out <- mergeWindows(obranges, sign = up.o, tol = tol)
            obranges <- out$region
            all.signs <- logical(length(obranges))
            all.signs[out$id] <- up.o
            up.o <- all.signs
        }
        kx <- read.table(known, header = TRUE)
        if (checkfc) {
            up.t <- kx$logFC > 0
            if (is.null(up.t)) {
                stop("need a log-FC field in the table of known sites")
            }
        } else {
            up.t <- rep(TRUE, nrow(kx))
        }
        if (!nrow(kx)) {
            stop("no known sites to check against")
        }
        kranges <- GRanges(kx[, 1], IRanges(kx[, 2], kx[, 3]))
        if (is.null(kx$name)) {
            kranges$name <- 1:length(kranges)
        } else {
            kranges$name <- kx$name
        }
        known.up <- kranges[up.t]
        known.down <- kranges[!up.t]
        u.olap <- findOverlaps(known.up, obranges[up.o])
        d.olap <- findOverlaps(known.down, obranges[!up.o])
        recall <-
            length(unique(known.up$name[queryHits(u.olap)])) +
            length(unique(known.down$name[queryHits(d.olap)]))
        overlapped <-
            length(unique(subjectHits(u.olap))) + length(unique(subjectHits(d.olap)))
        found <- length(obranges)
        
        # Calculating number of hits per peak/call (but first merging neighbouring windows/bins)
        all.olap <- findOverlaps(reduce(kranges), reduce(obranges))
        hitsperpeak <-
            as.data.table(all.olap)[, .N, by = 'queryHits'][, mean(N)]
        
        return(
            list(
                overlap = overlapped,
                found = found,
                recall = recall,
                time = time,
                callsize = mean(width(reduce(obranges))),
                peaksize = mean(width(reduce(kranges))),
                hitsperpeak = hitsperpeak,
                ncalls = length(reduce(obranges)),
                npeaks = length(reduce(kranges)),
                tpr = sum(width(intersect(
                    kranges, obranges
                ))) / sum(width(kranges)),
                fpr = sum(width(intersect(
                    setdiff(gr.gensize, kranges), obranges
                ))) / sum(width(setdiff(
                    gr.gensize, kranges
                ))),
                ppv = sum(width(intersect(
                    kranges, obranges
                ))) / sum(width(obranges)),
                error = 0
            )
        )
    } else{
        return(
            list(
                overlap = NA,
                found = NA,
                recall = NA,
                time = NA,
                callsize = NA,
                peaksize = NA,
                hitsperpeak = NA,
                ncalls = NA,
                npeaks = NA,
                tpr = NA,
                fpr = NA,
                ppv = NA,
                error = 1
            )
        )
    }
}

resultDump2 = function(regions,
                       tab,
                       cutoff = 0.05,
                       out = "out.tsv",
                       postprob = F,
                       is.error = F)
{
    if (is.error == F) {
        if (postprob == T) {
            # Sometimes RSEG gives a few windows with posterior probabilities >1. So, I am excluding these windows
            idx <- (tab$Postprob >= 0 & tab$Postprob <= 1)
            regions <- regions[idx]
            tab <- tab[idx, ]
            
            sig <-
                epigraHMM:::fdrControl(prob = tab$Postprob, fdr = cutoff)
        } else{
            sig <- (tab$FDR <= cutoff)
        }
        
        sig.bins <- regions[sig]
        write.table(
            data.frame(
                chr = as.character(seqnames(sig.bins)),
                start = start(sig.bins),
                end = end(sig.bins),
                tab[sig,]
            ),
            file = out,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
    } else{
        write.table(
            data.frame(
                chr = 'chrA',
                start = 0,
                end = 0
            ),
            file = out,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
    }
}

set.seed(2019 + iters)

# Simulation begins now

for (it in iters) {
    # Generating simulated data for histone mark data.
    up.pk <- 1:nde
    down.pk <- 1:nde + nde
    disp <- prior.df * dispersion / rchisq(npeaks, df = prior.df)
    
    if (!is.tf) {
        type.A.1 <- type.A.2 <- type.A.3 <- !logical(npeaks)
        chosen.drop.A <- sample(1:6, nde, replace = TRUE)
        type.A.1[up.pk] <- bitwAnd(chosen.drop.A, 0x1) > 0L
        type.A.2[up.pk] <- bitwAnd(chosen.drop.A, 0x2) > 0L
        type.A.3[up.pk] <- bitwAnd(chosen.drop.A, 0x4) > 0L
        
        type.B.1 <- type.B.2 <- type.B.3 <- !logical(npeaks)
        chosen.drop.B <- sample(1:6, nde, replace = TRUE)
        type.B.1[down.pk] <- bitwAnd(chosen.drop.B, 0x1) > 0L
        type.B.2[down.pk] <- bitwAnd(chosen.drop.B, 0x2) > 0L
        type.B.3[down.pk] <- bitwAnd(chosen.drop.B, 0x4) > 0L
    }
    
    # Single chromosome, for convenience.
    distances <- round(runif(npeaks, dist[1], dist[2]))
    pos.1 <- cumsum(distances)
    if (!is.tf) {
        pos.2 <- pos.1 + true.width * dist.peak[1]
        pos.3 <- pos.1 + true.width * dist.peak[2]
        sizes <- c(chrA = max(pos.3) + dist[2])
    } else {
        sizes <- c(chrA = max(pos.1) + dist[2])
    }
    chrs <- rep("chrA", npeaks)
    
    fnames <- list()
    cnames <- list()
    for (lib in 1:length(grouping)) {
        fname <- paste0(all.fix, "_out_", lib, ".sam")
        cname <- paste0('control', iters, "_out_", lib, ".sam")
        if (!is.tf) {
            # Simulating complex histone marks.
            if (grouping[lib] == "A") {
                drop.1 <- type.A.1
                drop.2 <- type.A.2
                drop.3 <- type.A.3
            } else {
                drop.1 <- type.B.1
                drop.2 <- type.B.2
                drop.3 <- type.B.3
            }
            
            # Generating ChIP reads
            peakFile(
                fname,
                chrs = chrs[drop.1],
                pos = pos.1[drop.1],
                mu = base.mu,
                disp = disp[drop.1],
                sizes = sizes,
                fraglen = fraglen,
                width = true.width,
                tf = FALSE
            )
            peakFile(
                fname,
                chrs = chrs[drop.2],
                pos = pos.2[drop.2],
                mu = base.mu,
                disp = disp[drop.2],
                sizes = sizes,
                fraglen = fraglen,
                width = true.width,
                tf = FALSE,
                append = TRUE
            )
            peakFile(
                fname,
                chrs = chrs[drop.3],
                pos = pos.3[drop.3],
                mu = base.mu,
                disp = disp[drop.3],
                sizes = sizes,
                fraglen = fraglen,
                width = true.width,
                tf = FALSE,
                append = TRUE
            )
            
            # Generating Control reads
            peakFile(
                cname,
                chrs = chrs[drop.1],
                pos = pos.1[drop.1],
                mu = mean(back.mu),
                disp = rep(dispersion, sum(drop.1)),
                sizes = sizes,
                fraglen = fraglen,
                width = true.width,
                tf = FALSE
            )
            peakFile(
                cname,
                chrs = chrs[drop.2],
                pos = pos.2[drop.2],
                mu = mean(back.mu),
                disp = rep(dispersion, sum(drop.2)),
                sizes = sizes,
                fraglen = fraglen,
                width = true.width,
                tf = FALSE,
                append = TRUE
            )
            peakFile(
                cname,
                chrs = chrs[drop.3],
                pos = pos.3[drop.3],
                mu = mean(back.mu),
                disp = rep(dispersion, sum(drop.3)),
                sizes = sizes,
                fraglen = fraglen,
                width = true.width,
                tf = FALSE,
                append = TRUE
            )
            
        } else {
            # Simulating simple TF changes.
            cur.mu <- rep(base.mu, npeaks)
            if (grouping[lib] == "A") {
                cur.mu[down.pk] <- down.mu
                cur.mu[up.pk] <- up.mu
            } else {
                cur.mu[up.pk] <- down.mu
                cur.mu[down.pk] <- up.mu
            }
            # Generating ChIP reads
            peakFile(
                fname,
                chrs = chrs,
                pos = pos.1,
                mu = cur.mu,
                disp = disp,
                sizes = sizes,
                fraglen = fraglen,
                width = true.width,
                tf = TRUE
            )
            # Generating Control reads
            peakFile(
                cname,
                chrs = chrs,
                pos = pos.1,
                mu = rep(mean(back.mu), length(cur.mu)),
                disp = rep(dispersion, length(disp)),
                sizes = sizes,
                fraglen = fraglen,
                width = true.width,
                tf = TRUE
            )
        }
        fnames[[lib]] <- fname
        cnames[[lib]] <- cname
    }
    
    fnames <- unlist(fnames)
    addBackground(
        fnames,
        sizes = sizes,
        width = back.width,
        rlen = back.rlen,
        back.mu  =  back.mu,
        dispersion = dispersion,
        prior.df = prior.df,
        append = TRUE
    )
    bam.files <- crunch2BAM(fnames)
    unlink(fnames)
    
    cnames <- unlist(cnames)
    addBackground(
        cnames,
        sizes = sizes,
        width = back.width,
        rlen = back.rlen,
        back.mu  =  back.mu,
        dispersion = dispersion,
        prior.df = prior.df,
        append = TRUE
    )
    bam.control.files <- crunch2BAM(cnames)
    unlink(cnames)
    
    lfile <- paste0(all.fix, "_log.txt")
    if (is.tf) {
        write.table(
            file = lfile,
            data.frame(
                chr = chrs[up.pk],
                start = pos.1[up.pk] - radius,
                end = pos.1[up.pk] + radius,
                logFC = 1
            ),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE
        )
        write.table(
            file = lfile,
            data.frame(chrs[down.pk], pos.1[down.pk] - radius, pos.1[down.pk] + radius, logFC = -1),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            append = TRUE,
            col.names = FALSE
        )
    } else {
        write.table(
            file = lfile,
            data.frame(
                chr = chrs[!type.A.1],
                start = pos.1[!type.A.1] - radius,
                end = pos.1[!type.A.1] + radius,
                logFC = -1,
                name = which(!type.A.1)
            ),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE
        )
        write.table(
            file = lfile,
            data.frame(
                chr = chrs[!type.A.2],
                start = pos.2[!type.A.2] - radius,
                end = pos.2[!type.A.2] + radius,
                logFC = -1,
                name = which(!type.A.2)
            ),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            append = TRUE,
            col.names = FALSE
        )
        write.table(
            file = lfile,
            data.frame(
                chr = chrs[!type.A.3],
                start = pos.3[!type.A.3] - radius,
                end = pos.3[!type.A.3] + radius,
                logFC = -1,
                name = which(!type.A.3)
            ),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            append = TRUE,
            col.names = FALSE
        )
        write.table(
            file = lfile,
            data.frame(
                chr = chrs[!type.B.1],
                start = pos.1[!type.B.1] - radius,
                end = pos.1[!type.B.1] + radius,
                logFC = 1,
                name = which(!type.B.1)
            ),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            append = TRUE,
            col.names = FALSE
        )
        write.table(
            file = lfile,
            data.frame(
                chr = chrs[!type.B.2],
                start = pos.2[!type.B.2] - radius,
                end = pos.2[!type.B.2] + radius,
                logFC = 1,
                name = which(!type.B.2)
            ),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            append = TRUE,
            col.names = FALSE
        )
        write.table(
            file = lfile,
            data.frame(
                chr = chrs[!type.B.3],
                start = pos.3[!type.B.3] - radius,
                end = pos.3[!type.B.3] + radius,
                logFC = 1,
                name = which(!type.B.3)
            ),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            append = TRUE,
            col.names = FALSE
        )
    }
    
    if (iters == 1) {
        system(paste('cp', lfile, results))
    }
    peakdir <- paste0(all.fix, "_peaks")
    dir.create(peakdir, showWarnings = FALSE)
    prefix <- sub("\\.bam$", "", basename(bam.files))
    prefix.control <-
        sub("\\.bam$", "", basename(bam.control.files))
    gsize <- sum(as.numeric(sizes))
    
    #############################################################
    ### Running HOMER, DiffBind, and ChIPComp
    #############################################################
    
    for (peakcaller in c("HOMER")) {
        all.peakfiles <- list()
        if (peakcaller == "HOMER") {
            timeHOMER <- microbenchmark({
                pktype <- "bed"
                hpdir <- file.path(peakdir, "homerpool")
                dir.create(hpdir)
                for (x in 1:length(bam.files)) {
                    pooldir <- file.path(hpdir, prefix[x])
                    dir.create(pooldir)
                    poolcontroldir <-
                        file.path(hpdir, prefix.control[x])
                    dir.create(poolcontroldir)
                    
                    system(
                        paste(
                            "makeTagDirectory",
                            pooldir,
                            bam.files[x],
                            "-format sam -keepAll"
                        )
                    )
                    system(
                        paste(
                            "makeTagDirectory",
                            poolcontroldir,
                            bam.control.files[x],
                            "-format sam -keepAll"
                        )
                    )
                    
                    outers <-
                        file.path(peakdir,
                                  paste0("homer_", prefix[x], ".txt"))
                    system(
                        paste(
                            "findPeaks",
                            pooldir,
                            "-i",
                            poolcontroldir,
                            "-style",
                            ifelse(is.tf, "factor", "histone"),
                            "-o",
                            outers,
                            "-gsize",
                            gsize,
                            "-fragLength",
                            fraglen,
                            "-inputFragLength",
                            fraglen,
                            "-tbp",
                            0
                        )
                    )
                    
                    # Converting to proper BED format.
                    blah <- read.table(outers)
                    write.table(
                        file = outers,
                        blah[, c(2, 3, 4, 1, 8)],
                        row.names = FALSE,
                        quote = FALSE,
                        sep = "\t",
                        col.names = FALSE
                    )
                    all.peakfiles[[x]] <- outers
                }
                unlink(hpdir, recursive = TRUE)
            }, times = 1)
        }
        # MACS was not able to call peaks from simulated data, so we are commenting
        # out the code chunks from this peak caller. Also, we did not explore SICER in this paper.
        # else if (peakcaller=="MACS") {
        #     pktype <- "macs"
        #     for (x in 1:length(bam.files)) {
        #         oprefix <- file.path(peakdir, paste0("macs_", prefix[x]))
        #         runMACS2(bam.files[x], oprefix, fraglen=fraglen, gsize=gsize, format="BAM",
        #                  extra = paste(ifelse(is.tf,'','--broad'),'--control',bam.control.files[x]))
        #         all.peakfiles[[x]] <- paste0(oprefix, "_peaks.xls")
        #     }
        #     macs.peakfiles <- all.peakfiles # For use with DBChIP.
        # } else if (peakcaller=="SICER") {
        #     if (is.tf || it > 1L) { next } # Only doing it for histone mark data, and just once, for demonstration.
        #     pktype <- "bed"
        #     scdir <- file.path(peakdir, "sicerbed")
        #     dir.create(scdir)
        #     if (!file.exists("~/tmp")) { dir.create("~/tmp") }
        #     win <- 200
        #     gap <- win*2L
        #     eval <- 1
        #
        #     all.peakfiles <- list()
        #     for (x in 1:length(bam.files)) {
        #         cur.bed <- "out.bed"
        #         convertBamToBed(bam.files[x], file.path(scdir, cur.bed))
        #         system(paste("methods/SICER_V1.1/SICER/SICER-rb.sh", scdir, cur.bed, scdir, "Simulated",
        #                      2147483647, win, fraglen, 1, gap, eval))
        #         newname <- file.path(peakdir, paste0("sicer_", x, ".bed"))
        #         blah <- read.table(file.path(scdir, sprintf("out-W%i-G%i-E%i.scoreisland", win, gap, eval)))
        #         write.table(data.frame(blah[,1:3], paste0("peak-", 1:nrow(blah)), blah[,4]), file=newname,
        #                     sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        #         all.peakfiles[[x]] <- newname
        #     }
        #     unlink("chr.list")
        #     unlink(scdir, recursive=TRUE)
        # }
        
        #############################################################
        # Running MACS with DiffBind
        #############################################################
        
        timeDiffBind <- microbenchmark({
            current <- dba(
                sampleSheet = data.frame(
                    SampleID = prefix,
                    Condition = grouping,
                    bamReads = bam.files,
                    ControlID = prefix.control,
                    bamControl = bam.control.files,
                    Peaks = unlist(all.peakfiles),
                    PeakCaller = pktype
                ),
                minOverlap = 2
            )
            current <-
                dba.count(current,
                          fragmentSize = fraglen,
                          bRemoveDuplicates = FALSE)
            current <-
                dba.contrast(
                    current,
                    group1 = current$masks$A,
                    group2 = current$masks$B,
                    name1 = "A",
                    name2 = "B",
                    minMembers = 2
                )
            current <-
                dba.analyze(DBA = current,
                            method = DBA_DESEQ2) # Avoiding tagwise.dispersion when data is homoskedastic.
            out <- dba.report(current, th = 1)
            
            curtab <- as.data.frame(elementMetadata(out))
            names(curtab)[names(curtab) == "Fold"] <- "logFC"
            
        }, times = 1)
        
        if (iters == 1) {
            write.table(
                as.data.frame(out),
                file = paste0(results, 'DiffBind_ranges.txt'),
                row.names = F,
                quote = F,
                sep = '\t'
            )
            write.table(
                as.data.frame(curtab),
                file = paste0(results, 'DiffBind_table.txt'),
                row.names = F,
                quote = F,
                sep = '\t'
            )
        }
        
        for (cutoff in fdr.thresholds) {
            resultDump2(out, curtab, cutoff, out = ofile)
            result <-
                assessChIP2(
                    ofile,
                    lfile,
                    tol = NA,
                    checkfc = FALSE,
                    gsize = gsize,
                    time = timeDiffBind$time / convtime + timeHOMER$time / convtime
                ) # Don't bother checking fold change, as it isn't well defined for complex events.
            dump2file(paste("DiffBind +", peakcaller), cutoff, result)
        }
        
        #############################################################
        ### Running ChIPComp
        ############################################################
        
        timeChIPComp <- microbenchmark({
            for (peak in unlist(all.peakfiles)) {
                write.table(
                    fread(peak)[, 1:3],
                    file = gsub('.txt', '.bed', peak),
                    col.names = F,
                    row.names = F,
                    quote = F,
                    sep = '\t'
                )
            }
            
            conf <- data.frame(
                SampleID = 1:4,
                condition = c('A', 'A', 'B', 'B'),
                factor = rep('Simulated', 4),
                ipReads = bam.files,
                ctReads = bam.control.files,
                peaks = gsub('.txt', '.bed', unlist(all.peakfiles))
            )
            
            dsg = as.data.frame(lapply(conf[, c("condition", "factor")], as.numeric)) -
                1
            dsg = as.data.frame(model.matrix( ~ condition, dsg))
            
            countSet = makeCountSet(
                conf,
                dsg,
                filetype = "bam",
                binsize = binsize,
                species = 'other'
            )
            countSet = ChIPComp(countSet)
            
            gr.ChIPComp = with(countSet$db, GRanges(chr, IRanges(start, end)))
            dt.ChIPComp = as.data.table(countSet$db)
            dt.ChIPComp[, FDR := p.adjust(pvalue.wald, method = 'BH')]
        }, times = 1)
        
        if (iters == 1) {
            write.table(
                as.data.frame(gr.ChIPComp),
                file = paste0(results, 'ChIPComp_ranges.txt'),
                row.names = F,
                quote = F,
                sep = '\t'
            )
            write.table(
                as.data.frame(dt.ChIPComp),
                file = paste0(results, 'ChIPComp_table.txt'),
                row.names = F,
                quote = F,
                sep = '\t'
            )
        }
        
        for (cutoff in fdr.thresholds) {
            resultDump2(gr.ChIPComp, dt.ChIPComp, cutoff, out = ofile)
            result <-
                assessChIP2(
                    ofile,
                    lfile,
                    tol = NA,
                    checkfc = FALSE,
                    gsize = gsize,
                    time = timeChIPComp$time / convtime + timeHOMER$time / convtime
                ) # Don't bother checking fold change, as it isn't well defined for complex events.
            dump2file(paste("ChIPComp +", peakcaller), cutoff, result)
        }
    }
    
    #############################################################
    ### Running diffReps
    ############################################################
    
    timediffReps <- microbenchmark({
        drdir <- file.path(peakdir, "diffreps")
        dir.create(drdir)
        
        for (x in 1:length(bam.files)) {
            convertBamToBed(bam.files[x], file.path(drdir, gsub(
                './', '', gsub('.bam', '.bed', bam.files[x])
            )))
        }
        for (x in 1:length(bam.control.files)) {
            convertBamToBed(bam.control.files[x], file.path(drdir, gsub(
                './', '', gsub('.bam', '.bed', bam.control.files[x])
            )))
        }
        
        write.table(
            data.frame('chrA', gsize),
            file = paste0(drdir, '/gsize.txt'),
            col.names = F,
            row.names = F,
            quote = F,
            sep = '\t'
        )
        
        bedfiles = list.files(drdir, pattern = 'hist.*.bed', full.names = T)
        bedcontrolfiles = list.files(drdir, pattern = 'control.*.bed', full.names = T)
        
        cmd = paste(
            'diffReps.pl --chrlen ',
            paste0(drdir, '/gsize.txt'),
            ' --report',
            paste0(drdir, '/diffreps.txt'),
            '--treatment',
            paste(bedfiles[1:2], collapse = ' '),
            '--control',
            paste(bedfiles[3:4], collapse = ' '),
            '--btr',
            paste(bedcontrolfiles[1:2], collapse = ' '),
            '--bco',
            paste(bedcontrolfiles[3:4], collapse = ' '),
            '--pval 1',
            '--nsd',
            ifelse(is.tf, 'sharp', 'broad'),
            '--window',
            binsize,
            '--meth nb',
            '--frag',
            fraglen,
            '--nohs --noanno'
        )
        cat('Command: ', cmd)
        system(cmd)
        
        tb.diffreps <-
            fread(paste0(drdir, '/diffreps.txt'), skip = 32)
        print(head(tb.diffreps))
        tb.diffreps[, FDR := padj]
        gr.diffreps <-
            with(tb.diffreps, GRanges(Chrom, IRanges(Start, End)))
    }, times = 1)
    
    if (iters == 1) {
        write.table(
            as.data.frame(gr.diffreps),
            file = paste0(results, 'diffReps_ranges.txt'),
            row.names = F,
            quote = F,
            sep = '\t'
        )
        write.table(
            as.data.frame(tb.diffreps),
            file = paste0(results, 'diffReps_table.txt'),
            row.names = F,
            quote = F,
            sep = '\t'
        )
    }
    
    for (cutoff in fdr.thresholds) {
        resultDump2(gr.diffreps, tb.diffreps, cutoff, out = ofile)
        result <-
            assessChIP2(
                ofile,
                lfile,
                tol = NA,
                checkfc = FALSE,
                gsize = gsize,
                time = timediffReps$time / convtime
            )
        dump2file('diffReps', cutoff, result)
    }
    unlink(drdir, recursive = TRUE)
    
    #############################################################
    ### Running RSEG
    ############################################################
    
    timeRSEG <- microbenchmark({
        drdir <- file.path(peakdir, "rseg")
        dir.create(drdir)
        
        # Pooling
        cmd = paste(
            'samtools merge',
            paste0(drdir, '/histone_A.bam'),
            paste0(bam.files[1:2], collapse = ' ')
        )
        cat(paste('Command:', cmd, '\n'))
        system(cmd)
        
        cmd = paste(
            'samtools merge',
            paste0(drdir, '/histone_B.bam'),
            paste0(bam.files[3:4], collapse = ' ')
        )
        cat(paste('Command:', cmd, '\n'))
        system(cmd)
        
        # Transforming to BED and sorting
        for (x in c(paste0(drdir, '/histone_A.bam'),
                    paste0(drdir, '/histone_B.bam'))) {
            convertBamToBed(x, gsub('.bam', '.bed', x))
            system('export LC_ALL=C')
            system(paste(
                'sort -k1,1 -k3,3n -k2,2n -k6,6r',
                gsub('.bam', '.bed', x),
                '-o',
                gsub('.bam', '.sorted.bed', x)
            ))
        }
        
        bedfiles = list.files(drdir, pattern = '*.sorted.bed', full.names = T)
        
        # Genome size and deadzones
        write.table(
            data.frame('chrA', 1, gsize),
            file = paste0(drdir, '/gsize.txt'),
            col.names = F,
            row.names = F,
            quote = F,
            sep = '\t'
        )
        write.table(
            data.frame('chrA', 1, 2),
            file = paste0(drdir, '/deadzones.txt'),
            col.names = F,
            row.names = F,
            quote = F,
            sep = '\t'
        )
        
        # Calling RSEG
        cmd = paste(
            '~/rseg/bin/rseg-diff -verbose -mode 3',
            '-out',
            paste0(drdir, '/rseg_output.txt'),
            '-score',
            paste0(drdir, '/rseg_postprob.txt'),
            '-b',
            binsize,
            '-chrom',
            paste0(drdir, '/gsize.txt'),
            '-deadzones',
            paste0(drdir, '/deadzones.txt'),
            '-duplicates',
            '-fragment_length',
            fraglen,
            bedfiles[1],
            bedfiles[2]
        )
        cat('Command: ', cmd)
        system(cmd)
    }, times = 1)
    
    if (file.exists(paste0(drdir, '/rseg_postprob.txt'))) {
        tb.rseg <- fread(paste0(drdir, '/rseg_postprob.txt'))
        tb.rseg[, Postprob := (V4 + V5)]
        gr.rseg <- with(tb.rseg, GRanges(V1, IRanges(V2, V3)))
        
        if (iters == 1) {
            write.table(
                as.data.frame(gr.rseg),
                file = paste0(results, 'RSEG_ranges.txt'),
                row.names = F,
                quote = F,
                sep = '\t'
            )
            write.table(
                as.data.frame(tb.rseg),
                file = paste0(results, 'RSEG_table.txt'),
                row.names = F,
                quote = F,
                sep = '\t'
            )
        }
        
        for (cutoff in fdr.thresholds) {
            resultDump2(gr.rseg,
                        tb.rseg,
                        cutoff,
                        out = ofile,
                        postprob = T)
            result <-
                assessChIP2(
                    ofile,
                    lfile,
                    tol = NA,
                    checkfc = FALSE,
                    gsize = gsize,
                    time = timeRSEG$time / convtime
                )
            dump2file('RSEG', cutoff, result)
        }
    } else{
        #If RSEG throws errors, flag it
        tb.rseg <- data.table('chrA', 0, 0)
        gr.rseg <- with(tb.rseg, GRanges(V1, IRanges(V2, V3)))
        
        for (cutoff in fdr.thresholds) {
            resultDump2(
                gr.rseg,
                tb.rseg,
                cutoff,
                out = ofile,
                postprob = T,
                is.error = T
            )
            result <-
                assessChIP2(
                    ofile,
                    lfile,
                    tol = NA,
                    checkfc = FALSE,
                    gsize = gsize,
                    time = timeRSEG$time / convtime,
                    is.error = T
                )
            dump2file('RSEG', cutoff, result)
        }
    }
    unlink(drdir, recursive = TRUE)
    
    #############################################################
    ### csaw, with its combined window methodology.
    ############################################################
    
    timecsaw <- microbenchmark({
        xparam <- readParam(dedup = FALSE)
        
        data <-
            windowCounts(
                bam.files,
                width = binsize,
                ext = fraglen,
                param = xparam,
                filter = 20
            )
        binned <-
            windowCounts(bam.files,
                         width = 2000,
                         bin = TRUE,
                         param = xparam)
        
        bin.ab <-
            scaledAverage(asDGEList(binned),
                          scale = median(getWidths(binned)) / median(getWidths(data)))
        threshold <- median(bin.ab) + log2(2)
        keep <- aveLogCPM(asDGEList(data)) >  threshold
        
        data <- data[keep,]
        tabres <-
            analyzeQLCounts(assay(data), design, totals = data$totals)
        merged <-
            mergeWindows(rowRanges(data),
                         tol = 100,
                         max.width = 5000)
        tabneg <- combineTests(merged$id, tabres)
        
    }, times = 1)
    
    if (iters == 1) {
        write.table(
            as.data.frame(merged$region),
            file = paste0(results, 'csaw_ranges.txt'),
            row.names = F,
            quote = F,
            sep = '\t'
        )
        write.table(
            as.data.frame(tabneg),
            file = paste0(results, 'csaw_table.txt'),
            row.names = F,
            quote = F,
            sep = '\t'
        )
    }
    
    for (cutoff in fdr.thresholds) {
        resultDump2(merged$region, tabneg, cutoff, out = ofile)
        result <-
            assessChIP2(
                ofile,
                lfile,
                tol = NA,
                checkfc = FALSE,
                gsize = gsize,
                time = timecsaw$time / convtime
            )
        dump2file("csaw", cutoff, result)
    }
    
    #############################################################
    ### Running mixNBHMM
    ############################################################
    
    timeepigraHMM <- microbenchmark({
        out.epigraHMM <-
            epigraHMM::epigraHMMDataSetFromBam(
                bamFiles = bam.files,
                colData = data.frame(
                    condition = c('1', '1', '2', '2'),
                    replicate = c(1, 2, 1, 2)
                ),
                genome = GenomicRanges::GRanges('chrA', IRanges::IRanges(1, gsize)),
                windowSize = binsize,
                gapTrack = FALSE,
                blackList = FALSE
            )
        
        control.epigraHMM <-
            epigraHMM::controlEM(
                epsilonEM = c(
                    MRCPE = 1e-03,
                    MACPE = 1e-03,
                    ARCEL = 1e-04
                ),
                tempDir = './',
                fileName = paste0('autosim', it)
            )
        
        out.epigraHMM <-
            epigraHMM::initializer(object = out.epigraHMM,
                                   control = control.epigraHMM)
        
        out.epigraHMM <-
            epigraHMM::epigraHMM(
                object = out.epigraHMM,
                control = control.epigraHMM,
                type = 'differential',
                dist = 'nb'
            )
        
        dt.epigraHMM <-
            as.data.table(SummarizedExperiment::rowRanges(out.epigraHMM))
        dt.epigraHMM[, Postprob := exp(rhdf5::h5read(
            S4Vectors::metadata(out.epigraHMM)$output,
            'logProb1'
        )[, 2])]
    }, times = 1)
    
    if (iters == 1) {
        write.table(
            as.data.frame(SummarizedExperiment::rowRanges(out.epigraHMM)),
            file = paste0(results, 'epigraHMM_ranges.txt'),
            row.names = F,
            quote = F,
            sep = '\t'
        )
        write.table(
            as.data.frame(dt.epigraHMM),
            file = paste0(results, 'epigraHMM_table.txt'),
            row.names = F,
            quote = F,
            sep = '\t'
        )
    }
    
    for (cutoff in fdr.thresholds) {
        resultDump2(
            SummarizedExperiment::rowRanges(out.epigraHMM),
            dt.epigraHMM,
            cutoff,
            out = ofile,
            postprob = T
        )
        result <-
            assessChIP2(
                ofile,
                lfile,
                tol = NA,
                checkfc = FALSE,
                gsize = gsize,
                time = timeepigraHMM$time / convtime
            )
        dump2file("epigraHMM", cutoff, result)
    }
    
    #############################################################
    ### Running THOR
    ############################################################
    
    timeTHOR <- microbenchmark({
        drdir <- file.path(peakdir, "thor")
        dir.create(drdir)
        
        write.table(
            data.frame('chrA', gsize),
            file = paste0(drdir, '/gsize.txt'),
            col.names = F,
            row.names = F,
            quote = F,
            sep = '\t'
        )
        write.table(
            rbindlist(
                list(
                    data.table(c('#rep1', bam.files[1:2])),
                    data.table(c('#rep2', bam.files[3:4])),
                    data.table(c('#inputs1', bam.control.files[1:2])),
                    data.table(c('#inputs2', bam.control.files[3:4])),
                    data.table(c(
                        '#chrom_sizes', paste0(drdir, '/gsize.txt')
                    ))
                )
            ),
            file = paste0(drdir, '/THOR.config'),
            col.names = F,
            row.names = F,
            quote = F
        )
        
        cmd = paste(
            'rgt-THOR',
            paste0(drdir, '/THOR.config'),
            '--name peaks -b',
            binsize,
            '--exts',
            paste0(rep(100, 4), collapse = ','),
            '--pvalue 1.0 --output-dir',
            drdir
        )
        cat('Command: ', cmd)
        system(cmd)
        
        dt.thor = as.data.table(read.table(
            paste0(drdir, '/peaks-diffpeaks.narrowPeak'),
            header = F
        ))
        dt.thor[, FDR := 10 ^ (-V8)]
        print(head(dt.thor))
        print(summary(dt.thor$FDR))
        gr.thor = with(dt.thor, GenomicRanges::GRanges(V1, IRanges(V2, V3)))
    }, times = 1)
    
    if (iters == 1) {
        write.table(
            as.data.frame(gr.thor),
            file = paste0(results, 'THOR_ranges.txt'),
            row.names = F,
            quote = F,
            sep = '\t'
        )
        write.table(
            as.data.frame(dt.thor),
            file = paste0(results, 'THOR_table.txt'),
            row.names = F,
            quote = F,
            sep = '\t'
        )
    }
    
    for (cutoff in fdr.thresholds) {
        resultDump2(gr.thor, dt.thor, cutoff, out = ofile)
        result <-
            assessChIP2(
                ofile,
                lfile,
                tol = NA,
                checkfc = FALSE,
                gsize = gsize,
                time = timeTHOR$time / convtime
            )
        dump2file('THOR', cutoff, result)
    }
    unlink(drdir, recursive = TRUE)
}

# Cleaning up.

if (iters > 1) {
    unlink(bam.files)
    unlink(paste0(bam.files, ".bai"))
    unlink(bam.control.files)
    unlink(paste0(bam.control.files, ".bai"))
    unlink(paste0('autosim', it, '.h5'))
}

unlink(ofile)