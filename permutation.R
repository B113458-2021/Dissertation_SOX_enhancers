#!/usr/bin/env Rscript

library(argparser, warn.conflicts=FALSE, quietly=TRUE)
argp = arg_parser("Perform circular permutation test.")
argp = add_argument(argp, "--inputA", help="Region of interest B (bed format)")
argp = add_argument(argp, "--inputB", help="Region of interest B (bed format)")
argp = add_argument(argp, "--sample", help="Sample name")
argp = add_argument(argp, "--ntimes", type="integer", help="Number of permutation tests to be performed, minimum 100")
argp = add_argument(argp, "--cores",type="integer", help="Number of cores used, min 8")
argv = parse_args(argp)

library(regioneR, warn.conflicts=FALSE, quietly=TRUE)

write(paste(Sys.time(), "Loading region of interest A = ",argv$inputA), stderr())
inputA <- sort(toGRanges(argv$inputA)) #load as a GRanges object and make sure it is sorted).

# Load region of interest to intersect
write(paste(Sys.time(), "Loading region of interest B = ",argv$inputB), stderr())
inputB <- sort(toGRanges(argv$inputB)) #load as a GRanges object and make sure it is sorted).

# Calculate overlaps
write(paste(Sys.time(), "Number of overlapping regions:"), stderr())
print(numOverlaps(inputA,inputB, count.once = TRUE))

# Perform permutation test:
write(paste(Sys.time(), "Performing circular permutation test."), stderr())

N <- argv$ntimes
print(N)
cores <- argv$cores
print(cores)


circ <- permTest(A=inputA, B=inputB, 
                    ntimes=N,
                    randomize.function=circularRandomizeRegions,
                    evaluate.function=numOverlaps,
                    count.once= TRUE,
                    genome="hg38",
                    mc.set.seed=F,
                    mc.cores=cores)
print(circ)

filename <- paste(argv$sample, "_circ_permutation.pdf", sep="")

# Plotting permutation graph
pdf(filename)
plot(circ, xlab="Number of overlaps")
dev.off()

local <- localZScore(A=inputA, B=inputB, pt=circ, window=1000, step=50, count.once=TRUE)

filename <- paste(argv$sample, "_circ_permutation.RData", sep="")
save.image(filename)
