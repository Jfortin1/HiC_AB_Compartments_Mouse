library(readr)
library(contact)
library(bsseq)
# Genome version: mm9

buildContactMatrix <- function(data, res=40*1000, chr = "chr14"){

	create.alignment.gr <- function(data, res=40*1000, chr="chr14"){
		start <- 0
		end <- (ncol(data))*res
		
		
		start.seq <- seq(start,end-res,res)
		end.seq   <- c(start.seq[-1], end)-1
		gr <- data.frame2GRanges(data.frame(chr=chr, start=start.seq, end=end.seq))
		gr
	}

	gr <- create.alignment.gr(data,res=res, chr=chr)
	contact <- Contact(matrix=data, gr=gr)
	contact
}

getCon <- function(file, chr){
	data <- read_tsv(file, col_names=FALSE)
	data <- as.matrix(data[,-ncol(data)])
	buildContactMatrix(data, chr=chr)
}

chrs <- paste0("chr", c(1:19, "X", "Y"))
data1 <- data2 <- NULL
for (jj in 1:length(chrs)){
	chr <- chrs[jj]

	file1 <- paste0("../raw/mES/uij/uij.", chr) 
	file2 <- paste0("../raw/mCO/uij/uij.", chr)
	files <- list(file1,file2)
	cons <- lapply(files, getCon, chr=chr)
	cons <- lapply(cons, removeLowCoverage)
	cons <- lapply(cons, iceNorm)
	cons <- lapply(cons, obsExpNorm)
	cons <- lapply(cons, function(x){
		x$matrix <- cor(x$matrix)
		x
	})
	
	# Common genome support for both cell types:
	gr <- contact:::intersect.granges(cons)
	cons <- lapply(cons, subsetByOverlaps, gr)

	# Obtaining the A/B compartments:
	pcs <- lapply(cons, getFirstPC)
	pcs <- lapply(pcs, meanSmoother, iter=2)
	pcs <- do.call(cbind, pcs)
	if (cor(pcs)[1,2] < 0){
		pcs[,2] <-  -pcs[,2]
	}

	# What is closed? What is open? 
	gcfile <- paste0("gc_40kb_", chr, ".rda")
	load(file.path("../misc/mm9_GC/", gcfile))
	gr <- cons[[1]]$gr
	indices <- subjectHits(findOverlaps(gr, gc))
	gc_track <- gc$gc[indices]

	for (j in 1:ncol(pcs)){
		temp <- sign(cor(pcs[,j], gc_track))
		if (temp==1){
			pcs[,j] <- -pcs[,j]
		}
	}

	# Formatting the data:
	gr <- cons[[1]]$gr
	pc <- pcs[,1]
	sign <- sign(pc)
	domain <- ifelse(sign==1, "closed", "open")
	df1 <- data.frame(chr=chr, start=start(gr), end=end(gr), eigen=pc, domain=domain)

	gr <- cons[[1]]$gr
	pc <- pcs[,2]
	sign <- sign(pc)
	domain <- ifelse(sign==1, "closed", "open")
	df2 <- data.frame(chr=chr, start=start(gr), end=end(gr), eigen=pc, domain=domain)


	print(jj)
}






# Plotting the A/B compartments:
par(mfrow=c(3,1), mar=c(2,2,2,2))
ylim=c(-0.1,0.1)
mybarplot(pcs[,1], ylim=ylim)
mybarplot(pcs[,2], ylim=ylim) 
mybarplot(gc_track-median(gc_track), ylim=c(-0.3,0.3))










