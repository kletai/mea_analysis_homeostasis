################TO READ IN DATA###########################

readInData = function(dataFile, PTXref, startTime, stopTime, stimTime){
  
  #to read in the data from a binned spike file made in MATLAB post spike-sorting (dataFile)
  #Also need a list of the units to use from python filtering in python (PTXref)
  #also input startTime of experiment, stopTime of experiment, and when stimulation was added (stimTime)
  #each dataset read in is a single treatment condition
  
  #EXAMPLE:
  #dataFile = '~/Dropbox (HMS)/for kelsey/Arc/311_spk_freq_table.csv'
  #PTXref = '~/Dropbox (HMS)/for kelsey/Arc/311_PTX_u.csv'
  #startTime = '04-Apr-2019 09:25:00'
  #stopTime = '05-Apr-2019 02:30:00'
  #stimTime = '04-Apr-2019 09:55:00'
  
  s = read.csv(file = dataFile)
  uPTX = read.csv(file=PTXref)
  uPTX = as.character(unique(uPTX$unit_name))
  
  start = which(s[,dim(s)[[2]]]==startTime)-35
  stop = which(s[,dim(s)[[2]]]==stopTime)
  
  first6h = s[(which(s[,dim(s)[[2]]]==startTime)-35):which(s[,dim(s)[[2]]]==startTime),1:(dim(s)[[2]]-1)]
  last6h = s[(which(s[,dim(s)[[2]]]==stopTime)-35):which(s[,dim(s)[[2]]]==stopTime),1:(dim(s)[[2]]-1)]
  afterStim = s[(which(s[,dim(s)[[2]]]==stimTime)):(which(s[,dim(s)[[2]]]==stimTime)+35),1:(dim(s)[[2]]-1)]
  
  first6h = first6h[,uPTX]
  last6h = last6h[,uPTX]
  afterStim = afterStim[,uPTX]

  return(list(first6h, last6h, afterStim,s[start:stop,uPTX]))
  

  
#########################TO MAKE MEDIAN TIME COURSE PLOTS####################################
  
timecourseplot = function(a, ylim = c(2/3, 3)){
  
  #for making figures from datasets (each 1 expt) generated in readInData function (e.g., Fig 1F)
  #a is a list of datasets
  #ylim controls y-axis of the plots
  
  #EXAMPLE:
  #timecourseplot(list(fl165,fl175, fl216, fl219, fl228, fl419, fl276, fl363, fl423), ylim=c(0,5))
    
    fulldata = sapply(a, '[[', 4)
    baselines = sapply(a, '[[', 1)
    
    normdata = list()
    for(i in 1:length(a)){
      x = t(t(fulldata[[i]])/apply(baselines[[i]], 2, mean))
      normdata[[i]] = x
    } #find fold change
    
    meds = lapply(normdata, apply, 1, median) #find medians
    medsM = matrix(unlist(meds), ncol = length(meds[[1]]), byrow = TRUE) 
    means = apply(medsM, 2, mean) #find mean of medians
    
    a = lowess(means, f=1/20) #smooth
    par(mar=c(4,4,2,2)) #plot
    plot(a$x*5/60 - 3, a$y, type='l', ylim = ylim, xlab = '', ylab = '',
         bty = 'l')
    mtext(text='Time after drug (h)', side=1, line=2.2)
    mtext(text='Fold change in firing rate', side=2, line=2.2)
    c = 1
    for(i in meds){
      b = lowess(i, f=1/20)
      lines(b$x*5/60 - 3, b$y, type='l', col = palette1[[c]], lwd=1)
      c = c+1
    }
    lines(a$x*5/60 - 3, a$y, type='l', ylim = c(0, 3), lwd = 1.5)
    abline(h=1, lty=2, col="grey", lwd =1.5)
}



timecourseplotSINGLE = function(expt, ylim = c(0,20)){
  
  #make a timecourse plot showing every unit in a single experiment (e.g., Figure 1D)
  #expt = experiment (from readInData)
  #ylim controls y-axis
  
  x = t(t(expt[[4]])/apply(expt[[1]], 2, mean))
  
  med = apply(x, 1, median)
  
  a = lowess(med, f=1/20)
  
  layout(matrix(1:2,ncol=1), width = c(20,20),height = c(15,3))
  par(mar=c(3,4,2,2))
  plot(a$x*5/60 - 3, a$y, type='l', ylim = ylim, xlab = '', ylab = '',
       bty = 'l')
  mtext(text='Time after drug (h)', side=1, line=2.2)
  mtext(text='Fold change in firing rate', side=2, line=2.2)
  
  for(i in 1:dim(x)[[2]]){
    b = lowess(x[,i], f=1/20)
    c = cut(log10(apply(expt[[1]], 2, mean)/300), breaks = seq(-3.1, log10(15), length.out = 100), len=100)
    colors = palette3[c]
    lines(b$x*5/60 - 3, b$y, type='l', col = colors[i])
  }
  lines(a$x*5/60 - 3, a$y, type='l', ylim = c(0, 3), lwd = 2)
  abline(h=1, lty=2, col="dark blue", lwd =1.5)
  
  legend_image <- as.raster(matrix(palette3, ncol=100))
  par(mar=c(0,4,0.5,2))
  plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
  text(y=0.4, x = seq(0.25,2.75,l=5), labels = signif(10^(seq(-3,log10(15),l=5)), digits=3))
  text(y=0.1, x = 1.5, labels='Firing Rate (Hz)')
  rasterImage(legend_image, 0.25, 0.6, 2.75,1)
}

timecourseplotTWOEXPT = function(a,b, ylim = c(2/3, 3), lab1 = 'Control', lab2 = 'KO', palette = palette2){
  
  #for making figures from two conditions on a single plot generated in readInData function (e.g., Fig 4B)
  #a is a list of datasets in condition 1 (e.g., KO)
  #b is a list of datasets in contdition 2 (e.g., het)
  #ylim controls y-axis of the plots
  
  fulldataA = sapply(a, '[[', 4)
  baselinesA = sapply(a, '[[', 1)
  
  fulldataB = sapply(b, '[[', 4)
  baselinesB = sapply(b, '[[', 1)
  
  normdataA = list()
  for(i in 1:length(a)){
    x = t(t(fulldataA[[i]])/apply(baselinesA[[i]], 2, mean))
    normdataA[[i]] = x
  }
  
  normdataB = list()
  for(i in 1:length(b)){
    x = t(t(fulldataB[[i]])/apply(baselinesB[[i]], 2, mean))
    normdataB[[i]] = x
  }
  
  medsA = lapply(normdataA, apply, 1, median)
  medsMA = matrix(unlist(medsA), ncol = length(medsA[[1]]), byrow = TRUE)
  meansA = apply(medsMA, 2, mean)
  
  medsB = lapply(normdataB, apply, 1, median)
  medsMB = matrix(unlist(medsB), ncol = length(medsB[[1]]), byrow = TRUE)
  meansB = apply(medsMB, 2, mean)
  
  a = lowess(meansA, f=1/20)
  b = lowess(meansB, f=1/20)
  
  par(mar=c(4,4,2,2))
  plot(a$x*5/60 - 3, a$y, type='l', ylim = ylim, xlab = '', ylab = '',
       bty = 'l', col = 'white')
  mtext(text='Time after drug (h)', side=1, line=2.2)
  mtext(text='Fold change in firing rate', side=2, line=2.2) 
  for(i in medsA){
    z = lowess(i, f=1/20)
    lines(z$x*5/60 - 3, z$y, type='l', col = palette[3], lwd=1)
  }
  
  
  for(i in medsB){
    z = lowess(i, f=1/20)
    lines(z$x*5/60 - 3, z$y, type='l', col = palette5[3], lwd=1)
  }
  
  lines(a$x*5/60 - 3, a$y, type='l', ylim = c(0, 3), lwd = 1.5, col = palette[5])
  lines(b$x*5/60 - 3, b$y, type='l', ylim = c(0, 3), lwd = 1.5, col = palette5[5])
  
  
  abline(h=1, lty=2, col="grey", lwd =1.5)
  
  legend(10, 3, lty=1, col=c(palette[5], palette5[5]), legend = c(lab1, lab2), bty='n')
}

#########################################TO MAKE VIOLIN PLOTS###############################
vioplot_multi_3 = function(l, i1 = 1, i2=3, i3=2, colpal='reds', lab1 = 'Baseline', lab2 = 'Post-Homeostasis'){
  
  #to make violin plots like Figure1H using 3 timepoint from the same experiments
  #l is a list of experiments from (readInData)
  #i1, i2, i3 control which time blocks will be used from experiment
  
  
  if(colpal == 'reds'){
    col1 = rgb(0.4, 0, 0.05, 0.7)
    col2 = rgb(0.984, 0.415, 0.39, 0.7)
    col3 = rgb(0.988, 0.733, 0.631, 0.7)
  }
  
  if(colpal == 'blues'){
    col1 = rgb(0.03, 0.19, 0.41, 0.7)
    col2 = rgb(0.78, 0.85, 0.94, 0.7)
  }
  
  
  library(vioplot)
  toplot1 = sapply(sapply(l, '[[', i1), function(x){apply(x, 2, mean)})
  toplot2 = sapply(sapply(l, '[[', i2), function(x){apply(x, 2, mean)})
  toplot3 = sapply(sapply(l, '[[', i3), function(x){apply(x, 2, mean)})
  toplot = list()
  for(i in 1:length(toplot1)){
    toplot[[(3*i-2)]] = toplot1[[i]]/300
    toplot[[3*i-1]] = toplot2[[i]]/300
    toplot[[3*i]] = toplot3[[i]]/300
  }
  par(mar=c(4,4,2,2), bty='l')
  plot("", xlim = c(0.5, length(l)*3+((length(l)-1))*0.5), ylim = c(log10(min(unlist(toplot))), log10(max(unlist(toplot)))+2), xaxt = "n",  bty='l',
       xlab='', ylab='')
  plot.loc = seq_along(toplot)
  for(i in 1:(length(l)-1)){
    plot.loc[(i*3+1):length(plot.loc)] = plot.loc[(i*3+1):length(plot.loc)] +0.5
  }
  lapply(seq_along(toplot), function(x)
    vioplot(log10(toplot[[x]]), col = rep(c(col1, col2, col3), length(l))[x], at = plot.loc[x],add = T)
  )
  mtext(text = "log10(Firing Rate (Hz))", side=2, line=2.2)
  labs = c()
  for(i in 1:(length(l))){
    labs = c(labs, paste('Replicate ', i))
  }
  axis(at=plot.loc[c(F,T,F)], side=1, labels = 1:length(l))
  mtext(text = "Replicate", side=1, line=2.2)
  
  p1 = c()
  for(i in 1:length(l)){
    p1=c(p1, ks.test(log10(toplot[[(i*3-2)]]), log10(toplot[[(i*3-1)]]))$p.value)
  }
  q1 = p.adjust(p1, 'fdr')
  print(q1)
  q1[which(q1<0.1)] = '*'
  q1[which(q1>0.1)] = ''
  
  p2 = c()
  for(i in 1:length(l)){
    p2=c(p2, ks.test(log10(toplot[[(i*3-2)]]), log10(toplot[[(i*3)]]))$p.value)
  }
  q2 = p.adjust(p2, 'fdr')
  print(q2)
  q2[which(q2<0.1)] = '*'
  q2[which(q2>0.1)] = ''
  
  
  for(i in 1:length(l)){
    ymax = log10(max(c(unlist(toplot[[(i*3-2)]]), unlist(toplot[[(i*3-1)]]), unlist(toplot[[(i*3)]]))))+1
    text(plot.loc[(i*3-1)], ymax-0.25, q1[[i]], cex = 1.5)
    text(plot.loc[(i*3)], ymax-0.75, q2[[i]], cex=1.5)
  }
  
  
}

vioplot_multi = function(l, i1 = 1, i2=2, colpal='reds', lab1 = 'Baseline', lab2 = 'Post-Homeostasis'){

  #to make violin plots like Figure4D using 2 timepoints from the same experiments
  #l is a list of experiments from (readInData)
  #i1, i2 control which time blocks will be used from experiment
  
  
  if(colpal == 'reds'){
    col1 = rgb(0.4, 0, 0.05, 0.7)
    col2 = rgb(0.988, 0.733, 0.631, 0.7)
  }
  
  if(colpal == 'blues'){
    col1 = rgb(0.03, 0.19, 0.41, 0.7)
    col2 = rgb(0.78, 0.85, 0.94, 0.7)
  }
  
  
  library(vioplot)
  toplot1 = sapply(sapply(l, '[[', i1), function(x){apply(x, 2, mean)})
  toplot2 = sapply(sapply(l, '[[', i2), function(x){apply(x, 2, mean)})
  toplot = list()
  for(i in 1:length(toplot1)){
    toplot[[(2*i-1)]] = toplot1[[i]]/300
    toplot[[2*i]] = toplot2[[i]]/300
  }
  plot.loc = seq_along(toplot)
  for(i in 1:(length(l)-1)){
    plot.loc[(i*2+1):length(plot.loc)] = plot.loc[(i*2+1):length(plot.loc)] +0.5
  }
  par(mar=c(4,4,2,2), bty='l')
  plot("", xlim = c(0.5, max(plot.loc)+0.5), ylim = c(log10(min(unlist(toplot))), log10(max(unlist(toplot)))+2), xaxt = "n",  bty='l',
       xlab='', ylab='')
  
  lapply(seq_along(toplot), function(x)
    vioplot(log10(toplot[[x]]), col = rep(c(col1, col2), length(l))[x], at = plot.loc[x],add = T)
  )
  mtext(text = "log10(Firing Rate (Hz))", side=2, line=2.2)
  axis(at=(plot.loc)[c(T,F)]+0.5, side=1, labels = 1:length(l))
  mtext(text = "Replicate", side=1, line=2.2)
  
  p1 = c()
  for(i in 1:length(l)){
    p1=c(p1, ks.test(log10(toplot[[(i*2-1)]]), log10(toplot[[(i*2)]]))$p.value)
  }
  q1 = p.adjust(p1, 'fdr')
  
  for(i in 1:length(l)){
    ymax = log10(max(c(unlist(toplot[[(i*2-1)]]), unlist(toplot[[(i*2)]]))))+1
    segments((plot.loc[i*2-1]), ymax-0.5, plot.loc[i*2], ymax-0.5, col='black')
    text((plot.loc[i*2]-0.5), ymax-0.25, 
         paste('q=',signif(q1[i],2)))
  }
  
  
}

vioplot_multi_2group = function(l1, l2, i1 = 1, lab1='Control', lab2='KO'){
  
  #to make violin plots like Figure4E using 1 timepoint from the 2 experimental conditions
  #l1 is a list of experiments from (readInData) for condition 1
  #l2 is a list of experiments from (readInData) for condition 2
  #i1 control which time block will be used from experiments
  
  col1 = palette2[5]
  col2 = palette5[5]
  
  par(mar=c(4,4,2,2))
  
  library(vioplot)
  toplot1 = sapply(sapply(l1, '[[', i1), function(x){apply(x, 2, mean)})
  toplot2 = sapply(sapply(l2, '[[', i1), function(x){apply(x, 2, mean)})
  toplot = list()
  for(i in 1:length(toplot1)){
    toplot[[(2*i-1)]] = toplot1[[i]]/300
    toplot[[2*i]] = toplot2[[i]]/300
  }
  plot.loc = seq_along(toplot)
  for(i in 1:(length(l1)-1)){
    plot.loc[(i*2+1):length(plot.loc)] = plot.loc[(i*2+1):length(plot.loc)] +0.5
  }
  par(mar=c(4,4,2,2), bty='l')
  plot("", xlim = c(0.5, max(plot.loc)+0.5), ylim = c(log10(min(unlist(toplot))), log10(max(unlist(toplot)))+2), xaxt = "n",  bty='l',
       xlab='', ylab='')
  #axis(1, labels = c("4cyl", "6cyl", "8cyl"), at = c(1:length(l)))
  
  lapply(seq_along(toplot), function(x)
    vioplot(log10(toplot[[x]]), col = rep(c(col1, col2), length(l))[x], at = plot.loc[x],add = T)
  )
  mtext(text = "log10(Firing Rate (Hz))", side=2, line=2.2)
  axis(at=(plot.loc)[c(T,F)]+0.5, side=1, labels = 1:length(l1))
  mtext(text = "Replicate", side=1, line=2.2)
  
  p1 = c()
  for(i in 1:length(l1)){
    p1=c(p1, ks.test(log10(toplot[[(i*2-1)]]), log10(toplot[[(i*2)]]))$p.value)
  }
  print(p1)
  q1 = p.adjust(p1, 'fdr')
  print(q1)
  
  for(i in 1:length(l1)){
    ymax = log10(max(c(unlist(toplot[[(i*2-1)]]), unlist(toplot[[(i*2)]]))))+1
    segments((plot.loc[i*2-1]), ymax-0.5, plot.loc[i*2], ymax-0.5, col='black')
    text((plot.loc[i*2]-0.5), ymax-0.25, 
         paste('q=',signif(q1[i],2)))
  }
  
  
}



########################CALCULATE TIME TO RETURN TO BASELINE#####################
return2baselineAnalysis = function(expt){
  #expt = an experiment from readInData
  
  library(TTR)
  row.names(expt[[4]]) = c(1:dim(expt[[4]])[[1]])
  m = t(t(expt[[4]])/apply(expt[[1]], 2, mean))
  meds = apply(m, 1, median)
  
  f7 <- rep(1/13,13)
  medsMA = filter(meds, f7, sides=2)
  plot(meds, col = 'grey', type = 'l')
  lines(medsMA, col = 'blue')
  medsMA[is.na(medsMA)] = 100
  names(medsMA) = 1:393
  m = min(runSD(medsMA, 72)[is.na(runSD(medsMA, 72))==F]) #running SD over 6h
  minSD = medsMA[runSD(medsMA, 72) == m][is.na(medsMA[runSD(medsMA, 72) == m])==F]*1.025 #find entry where SD is minimum (and give it 2.5% more)
  returnTime = min(as.numeric(names(medsMA[46:393][medsMA[46:393]<=minSD])))
  plot(meds, col = 'grey', type = 'l')
  abline(v=returnTime)
  abline(h=minSD)
  
  print(returnTime*5/60-3) #convert into hours after stim time, assuming 3h baseline
}
  