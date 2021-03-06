# We used this code to plot our results from analyses of 108-gene bins

# set working directory
setwd("/Your/Path/Here/")

library("ape")

# read data 
data <- read.table("Bins108_BS_support.txt", header=T, row.names=1)

data
data$MnemiopsisOtherMetazoa

row.names(data)

zz <- c(1:10)

# create cladograms to be plotted as guide/explanatory trees
CtenMet <- read.tree(text = "(Mnemiopsis,(Amphimedon,(Trichoplax,(Cnidaria,Bilateria))));")
AmphCten <- read.tree(text="(Amphimedon,(Mnemiopsis,(Trichoplax,(Cnidaria,Bilateria))));")
Coelent <- read.tree(text="(Amphimedon,(Trichoplax,((Mnemiopsis,Cnidaria),Bilateria)));")
Mandib <- read.tree(text="(Chelicerata,(Strigamia,Pancrustacea));")
Paradox <- read.tree(text="((Chelicerata,Strigamia),Pancrustacea);")
Strigam <- read.tree(text="(Strigamia,(Chelicerata,Pancrustacea));")

# open png file
png(file="108bins.png",    # create PNG for the heat map        
    width = 10*300,        # 10 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12)  

# defining the plot space,
# plotting guide trees,
# and support by rate of evolution

par(oma=c(2,2,2,2), mfrow=c(4,3))
par(mar=c(2.5,3,2.5,3))

plot(CtenMet, edge.width=4, cex=1.2, label.offset=0.2, font=1, x.lim=c(-1,8))
nodelabels(node=7, pch=21, bg="red", cex=2)
plot(AmphCten, edge.width=4, cex=1.2, label.offset=0.2, font=1, x.lim=c(-1,8))
nodelabels(node=7, pch=21, bg="red", cex=2)
plot(Coelent, edge.width=4, cex=1.2, label.offset=0.2, font=1, x.lim=c(-1,8))
nodelabels(node=9, pch=21, bg="red", cex=2)

par(mar=c(2,2,1.5,2))

plot(xlim=c(1,10),x=zz, ylim=c(0,100), y=data$MnemiopsisOtherMetazoa, "o-", mar=c(0,0,0,0), xaxt="n")
title(main="Ctenophore sister to Metazoa", cex.main=1.2, xlab='', ylab='')
axis(1, at=c(1:10))
plot(xlim=c(1,10),x=zz, ylim=c(0,100), y=data$MnemiopsisPlacozoaEumetazoa, "o-", mar=c(0,0,0,0), xaxt="n")
title(main="(Ctenophore,(Placozoa,Eumetazoa))", cex.main=1.2, xlab='', ylab='')
axis(1, at=c(1:10))
plot(xlim=c(1,10),x=zz, ylim=c(0,100), y=data$Coelenterata, "o-", mar=c(0,0,0,0), xaxt="n")
title(main="Coelenterata", cex.main=1.2, xlab='', ylab='')
axis(1, at=c(1:10))

par(mar=c(2.5,3,2.5,3))

plot(Mandib, edge.width=4, cex=1.2, label.offset=0.2, font=1, x.lim=c(-3,7))
nodelabels(node=5, pch=21, bg="red", cex=2)
plot(Paradox, edge.width=4, cex=1.2, label.offset=0.2, font=1, x.lim=c(-3,7))
nodelabels(node=5, pch=21, bg="red", cex=2)
plot(Strigam, edge.width=4, cex=1.2, label.offset=0.2, font=1, x.lim=c(-3,7))
nodelabels(node=5, pch=21, bg="red", cex=2)

par(mar=c(2,2,1.5,2))

plot(xlim=c(1,10),x=zz, ylim=c(0,100), y=data$Mandibulata, "o-", mar=c(0,0,0,0), xaxt="n")
title(main="Mandibulata", cex.main=1.2, xlab='', ylab='')
axis(1, at=c(1:10))
plot(xlim=c(1,10),x=zz, ylim=c(0,100), y=data$Paradoxopoda, "o-", mar=c(0,0,0,0), xaxt="n")
title(main="Paradoxopoda", cex.main=1.2, xlab='', ylab='')
axis(1, at=c(1:10))
plot(xlim=c(1,10),x=zz, ylim=c(0,100), y=data$StrigamiaChelicerataPancrustacea, "o-", mar=c(0,0,0,0), xaxt="n")
title(main="(Strigamia,(Chelicerata,Pancrustacea))", cex.main=1.2, xlab='', ylab='')
axis(1, at=c(1:10))

mtext("       Bootstrap support                                                   Bootstrap support", side=2,  outer=T, padj=-0.5, adj=0)
mtext("Bin number: from the slowest-evolving on the left to the fastest-evolving on the right", side=1, outer=T, padj= 0.5)
#mtext("108-gene bins", side=3, outer=T, cex=1.2)

# closing the png file
dev.off()
