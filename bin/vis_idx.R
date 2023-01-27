library(ggplot2)

fname <- 'hg38.cttaag.idx.gz'
bp_per_pixel <- 375

# read in the index
con = gzfile(fname, "r")

label_ref <- c()
label_num <- c()
label_pos <- c()
label_motif <- c()

cur_ref <- ''
idx <- 1
while (TRUE) {
    line <- readLines(con, n=1)
    if (length(line) == 0) {
      break
    }
    
    if (length(grep("#motif ", line))==1) {
        # ignore header...
    } else if (length(grep("^>", line))==1) {
        tmp <- gsub("^>", "", line)
        cur_ref <- strsplit(tmp, '\t')[[1]][1]
        ref_label_count <- as.integer(strsplit(tmp, '\t')[[1]][3])
        print(paste0(cur_ref, " (", ref_label_count, ")"))

        label_ref <- c(label_ref, vector("character", ref_label_count))
        label_num <- c(label_num, vector("character", ref_label_count))
        label_pos <- c(label_pos, vector("numeric", ref_label_count))
        label_motif <- c(label_motif, vector("numeric", ref_label_count))

    } else {
        label_ref[idx] <- cur_ref
        label_num[idx] <- as.numeric(strsplit(line, "\t")[[1]][1])
        label_pos[idx] <- as.numeric(strsplit(line, "\t")[[1]][2])
        label_motif[idx] <- as.numeric(strsplit(line, "\t")[[1]][3])
        idx <- idx + 1
    }
}
close(con)

ref.df <- data.frame(ref=label_ref, num=label_num, pos=label_pos, motif=label_motif)

ref.df$dist_right <- NA
ref.df$dist_left <- NA

ref.df[1:nrow(ref.df)-1,]$dist_right <- ref.df[2:nrow(ref.df),]$pos - ref.df[1:nrow(ref.df)-1,]$pos
ref.df[2:nrow(ref.df),]$dist_left <- ref.df[2:nrow(ref.df),]$pos - ref.df[1:nrow(ref.df)-1,]$pos


for (i in 1:nrow(ref.df)) {
    if (i < nrow(ref.df) && ref.df[i,]$ref != ref.df[i+1,]$ref) {
        ref.df[i,]$dist_right <- NA
    }
    if (i < 1 && ref.df[i,]$ref != ref.df[i-1,]$ref) {
        ref.df[i,]$dist_left <- NA
    }
}
ref.df$ref <- factor(ref.df$ref, levels=unique(label_ref))
summary(ref.df)

plot(density(log10(ref.df$dist_right), na.rm=T))

p <- ggplot(ref.df, aes(x=ref, y=log10(dist_right))) +
geom_boxplot()

p

boxplot(log10(ref.df$dist_right))


ref.df$pixel <- unlist(lapply(ref.df$pos, function(x) { floor(x / bp_per_pixel)}))

mol2648157 <- c(2251,11224,24325,26342,31524,34690,44733,47489,50589,53278,55013,56929,63792,68700,71844,80074,83351,90148,93208,107591,111844,131903,155636,165775,174338,176250)
ref_min <- signif(411639,2)*.9
ref_max <- signif(573270,2)*1.1




tmp.df <- ref.df[ref.df$ref=='chr12' & ref.df$pos >= ref_min & ref.df$pos <= ref_max,]
tmp.df$pixel_offset <- tmp.df$pixel-tmp.df[1,]$pixel

mol1_pixel <- (mol2648157/bp_per_pixel)

p <- ggplot(tmp.df, aes(x=pixel_offset), y=1) + 
    geom_point(y=1, color="#0000cc99") + 
    scale_y_continuous(limits=c(0, 2)) +
    theme_bw()

    #p <- p + geom_point(data=data.frame(pos=mol1_pixel, y=0.95), aes(pos, y=y), color='#cc000066')

for (i in 1:10) {
    mol1_pixel <- (mol1/bp_per_pixel)
    mol1_pixel <- mol1_pixel - mol1_pixel[i]
    p <- p + geom_point(data=data.frame(pos=mol1_pixel, y=1+(0.05*i)), aes(pos, y=y), color='#00cc0066')
}
p
