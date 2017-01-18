getGC <- function(chr="chr14", res=40*1000){
	file <- paste0(chr,"_GC_fragL200_bin200.txt")
	data <- read_tsv(file,col_types=c("dd"), col_names=FALSE)
	colnames(data) <- c("pos", "count")
	
	start <- data$pos
	end <- start[-1]-1
	n <- length(start)
	end <- c(end, end[n-1]+200)
	gr <- data.frame2GRanges(data.frame(chr=chr, start=start, end=end))

	size <- end[length(end)]+1
	start <- seq(0,size, res)
	end <- start[-1]-1
	end <- c(end, end[length(end)]+res)
	new.gr <- data.frame2GRanges(data.frame(chr=chr, start=start, end=end))

	indices <- subjectHits(findOverlaps(gr,new.gr))
	temp <- tapply(data[,2], INDEX=indices, FUN=mean)
	gc <- rep(0.5, length(new.gr))
	gc[as.numeric(names(temp))] <- temp

	new.gr$gc <- gc
	new.gr
}

chrs <- paste0("chr",c(1:19, "X", "Y", "M"))
for (i in 1:length(chrs)){
	gc <- getGC(chrs[i])
	new_file <- paste0("gc_40kb_", chrs[i], ".rda")
	save(gc, file=new_file)
	print(i)
}