
# make sure files are read in properly
options(stringsAsFactors=F)

# read in tables of reads subdivided by category
fs = list.files(full.name=T)
vs <- lapply(fs,read.table, header=TRUE)

# set up how the reads are going to be binned
breaks = seq(50,600,by=5)
breaks = c(0, breaks,10000)

# read in table of all reads which contains a header "x". Remove it with vim or read table as header=TRUE
gw = read.table("../gw.fragments.txt", header=TRUE)
gw$x <- as.numeric(gw$x)
# There is one read that is much too long - discard
gw <- subset(gw, x < 10000)

# histogram of all read lengths with bins of size 5 up to 800, then all larger
gw.h = hist(gw$x,breaks=breaks,plot=F)

# using this histogram
get_hist <- function(x,breaks,gw.h) {
    # calculate frequencies of reads of each length
    x$V2 <- as.numeric(x$V2)
    # In one category, there is one read that is much too long - discard
    x <- subset(x, V2 < 10000)
    h <- hist(x$V2,breaks=breaks,plot=F)
    # total number of reads
    t.reads = length(x$V2)
    bin.counts = h$counts
    # normalise to total number of reads
    bin.counts.n = bin.counts/t.reads
    # divide by genome wide number of reads (divided by total number of reads)
    en = bin.counts.n/(gw.h$counts/sum(gw.h$counts))
    # log10
    en2 = log10(en)
    # save list of bin sizes, frequencies, total reads, normalised counts, log counts
    o = list(
        breaks = h$breaks,
        counts = h$counts,
        TReads = t.reads,
        NormC  = bin.counts.n,
        BinMid = h$mids,
        Log10En = en2
    )
    return(o)
}

# apply this function for each category of reads
hs <- lapply(vs,get_hist,breaks,gw.h)

cols = c(
    "#550000","#FF0000",
    "#FF8181","#FFFFFF",
    "#8181FF","#0000FF",
    "#000055")

# make colour scale
cols <- colorRampPalette(c("blue", "white", "red"))(n = 100)
# get the log10 count values and concatenate, keep those less than 10e4
tmp = lapply(hs,function(x) return(x$Log10En))
tmp = do.call(c,tmp)
tmp = tmp[abs(tmp) < 10e4]
# record range
rangem = range(tmp,na.rm=T)

# breaks for colours in heat map
b = seq(-max(abs(rangem)),max(abs(rangem)),length.out=length(cols)+1)

pdf("enrichments2.pdf",width=14,height=10)
# set up blank plot
plot(
    0,type="n",
    xlim=c(0,length(breaks)+1),
    ylim=c(0,length(hs)-1),
    axes=F,
    xlab='Fragment Size',
    ylab=''
)
# draw axis
axis(1,
     at=0:(length(breaks)-1),
     labels=breaks)

# legend labelling
legend.labs = c()
# for each category
for( i in 1:length(hs) ) {
# labels
    lab = fs[i]
    lab = gsub("./","",lab)
    lab = gsub(".txt","",lab)
    legend.labs[i] = lab
    # assign colours by count values
    x = hs[[i]]$Log10En
    x2 = cut(x,b)
    x3 = cols[x2]
    x3[is.na(x3)] = '#d3d3d3'
  # draw each box of the heat map and colour appropriately
    rect(
        xleft=0:(length(breaks)-1),
        xright=1:length(breaks),
        ybottom=i-1,
        ytop=i,
        col=x3,
        border=x3
    )
}
# add labels

text(
  x=rep(120,length(legend.labs)),
  y=0:length(legend.labs)+.5,
  labels=legend.labs
)


dev.off()