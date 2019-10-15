#install.packages(c("ggplot2","gridExtra","gtable","tidyr","devtools"),dependencies=T)

# install pophelper package from GitHub
#devtools::install_github('royfrancis/pophelper')

# load library for use
library(pophelper)
library(ggplot2)
library(gridExtra)
library(gtable)
library(tidyr)
library(grid)

pdf("Faststructure_KKNPLB.pdf",width=6,height=4,paper="special",pointsize=15)
par(mar=c(4,2,2,2))
pops<-read.table("Ind_allCu.txt",stringsAsFactors=F)
pops<-pops[-37,]
data<-readQ("fS_run_K.5.meanQ")
rownames(data[[1]])<-pops$V1
#pops_V2<-data.frame(as.character(pops$V2),stringsAsFactors=FALSE)
colour_pops<-c("darkblue","seagreen3","dodgerblue","coral3","darkgreen")
Qplot<-plotQ(data,returndata=T,returnplot=T,exportplot=F,quiet=T,basesize=5,showgrplab=T,grplab=pops[,2,drop=F],selgrp="V2",clustercol=colour_pops,subsetgrp=c("Krom","Kosi","Noss","Pais","Lang","Best"),grplabsize=5,grplabpos=0.85,grplabheight=2,indlabheight=0,linecol="white",
divcol="black",divtype=1,barbordercol="black",barbordersize=0.2,linealpha=0,linepos=0,showyaxis=F,showsp=F)
y.breaks <- c(0, 0.25,0.5,0.75,1) # placeholder, since this wasn't specified in the question
y.axis.grob <- textGrob(label = y.breaks,
                        x = unit(0.75, "npc"),
                        y = unit(seq(0, 1, length.out = length(y.breaks)), "npc"),
                        just = "right")

# add text grob into main plot's gtable
Qplot$plot[[1]] <- gtable_add_grob(Qplot$plot[[1]],
                               y.axis.grob,
                               t = 1, 
                               l = 1, 
                               b = nrow(Qplot$plot[[1]]) - 1,
                               clip = "off")
grid.draw(Qplot$plot[[1]])

dev.off()