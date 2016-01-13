

options(stringsAsFactors=F)

fs = list.files(path='intersects',pattern='inter',full.name=T)
vs <- lapply(fs,read.table)

breaks = seq(0,800,by=5)
breaks = c(breaks,10000)

gw = read.table("gw_vals.txt")
gw.h = hist(gw$V1,breaks=breaks,plot=F)

get_hist <- function(x,breaks,gw.h) {
    h <- hist(x$V1,breaks=breaks,plot=F)
    t.reads = length(x$V1)
    bin.counts = h$counts
    bin.counts.n = bin.counts/t.reads
    en = bin.counts.n/(gw.h$counts/sum(gw.h$counts))
    en2 = log10(en)
    
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

hs <- lapply(vs,get_hist,breaks,gw.h)

cols = c(
    "#550000","#FF0000",
    "#FF8181","#FFFFFF",
    "#8181FF","#0000FF",
    "#000055")

cols <- colorRampPalette(c("blue", "white", "red"))(n = 100)

tmp = lapply(hs,function(x) return(x$Log10En))
tmp = do.call(c,tmp)
tmp = tmp[abs(tmp) < 10e4]
rangem = range(tmp,na.rm=T)

b = seq(-max(abs(rangem)),max(abs(rangem)),length.out=length(cols)+1)

png("enrichments.png",width=800,height=600)

plot(
    0,type="n",
    xlim=c(0,length(breaks)+1),
    ylim=c(0,length(hs)),
    axes=F,
    xlab='Fragment Size',
    ylab=''
)

axis(1,
     at=0:(length(breaks)-1),
     labels=breaks)

legend.labs = c()

for( i in 1:length(hs) ) {

    lab = fs[i]
    lab = gsub("inter_stats_","",lab)
    lab = gsub(".txt","",lab)
    lab = gsub("intersects/","",lab)
    legend.labs[i] = lab
    
    x = hs[[i]]$Log10En
    x2 = cut(x,b)
    x3 = cols[x2]
    x3[is.na(x3)] = '#d3d3d3'
  
    rect(
        xleft=0:(length(breaks)-1),
        xright=1:length(breaks),
        ybottom=i-1,
        ytop=i,
        col=x3,
        border=x3
    )
}

text(
    x=rep(120,length(legend.labs)),
    y=0:length(legend.labs)+.5,
    labels=legend.labs
)

dev.off()




