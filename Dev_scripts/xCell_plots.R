library(RColorBrewer)
library(gplots)
library(stats)
library(corrplot)
library(stringr)
library(psych)
library(pheatmap)
library(MCPcounter)

createMatrixNoDups = function(data) {
  M = as.matrix(data[,-1])
  samples = colnames(M)
  genes = data[,1]
  M = as.numeric(as.character(M))
  M = matrix(M,nrow=length(genes),byrow=FALSE)
  if (length(unique(genes))!=length(genes)) {
    M = aggregate(x = M, by = list(genes), FUN = mean)
    genes = M[,1]
    M = M[,-1]
  } 
  rownames(M) = genes
  colnames(M) = samples
  return (as.matrix(M))
}

mix.protocol = function(working.dir,mix_set,spill,alpha=1,do.ciber=TRUE,do.known=TRUE,do.mcp=FALSE,transpose=FALSE,ct_order=NULL,font.size=7,show.before=TRUE) {
  colors <- brewer.pal(11, "RdBu")
  m=(read.table(paste0(working.dir,'/simulations/mix_',mix_set,'_scores.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE))
  dist=(read.table(paste0(working.dir,'/simulations/mix_',mix_set,'_dist.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE))
  
  oclust='original'
  
  if(is.null(ct_order)) {
    ct_order = rownames(dist)
  }
  rxr=cor(t(m[rownames(dist),]),t(dist),method="pearson")
  pdf(paste0(working.dir,'/simulations/plots/',mix_set,'_raw.pdf'),width=4,height=4)
  rxrc = corrplot(as.matrix(rxr[ct_order,ct_order]),tl.cex = font.size/10,tl.col = "black",col=colorRampPalette(rev(colors))(200),order=oclust)
  dev.off()
  mixScatters(paste0(working.dir,'/simulations/plots/scatters/',mix_set,'_scatters_raw.pdf'),m[rownames(dist),],dist,norm=TRUE)
  
  M = transformScores(m,spill$fv)
  rxb=cor(t(M[rownames(dist),]),t(dist),method="pearson")
  pdf(paste0(working.dir,'/simulations/plots/',mix_set,'_transform.pdf'),width=4,height=4)
  rxbc = corrplot(as.matrix(rxb[ct_order,ct_order]),tl.cex = font.size/10,tl.col = "black",col=colorRampPalette(rev(colors))(200),order=oclust)
  dev.off()
  mixScatters(paste0(working.dir,'/simulations/plots/scatters/',mix_set,'_scatters_transform.pdf'),M[rownames(dist),],dist,norm=FALSE)

  M2 = spillOver(M,spill$K,alpha)
  rx=cor(t(M2[ct_order,]),t(dist[ct_order,]),method="pearson")
  pdf(paste0(working.dir,'/simulations/plots/',mix_set,'_spill.pdf'),width=4,height=4)
  corrplot(as.matrix(rx),tl.cex = font.size/10,tl.col = "black",col=colorRampPalette(rev(colors))(200),order=oclust)
  dev.off()
  A = intersect(rownames(M2),rownames(dist))
  mixScatters(paste0(working.dir,'/simulations/plots/scatters/',mix_set,'_scatters_spill.pdf'),M2[A,],dist[A,],norm=FALSE)
  mixScatters(paste0(working.dir,'/simulations/plots/scatters/spill/',mix_set,'_scatters_spill.pdf'),M2[A,],dist[A,],norm=FALSE)
  
  pdf(paste0(working.dir,'/simulations/plots/',mix_set,'_corrplots.pdf'),width=8,height=4)
  par(mfrow=c(1,2))
  rxbc=corrplot(as.matrix(rxb[ct_order,ct_order]),tl.cex = font.size/10,tl.col = "black",col=colorRampPalette(rev(colors))(200),order=oclust)
  rxba=corrplot(as.matrix(rx[ct_order,ct_order]),tl.cex = font.size/10,tl.col = "black",col=colorRampPalette(rev(colors))(200),order=oclust)
  dev.off()
  
  ddiag = sum(diag(rxba))/sum(diag(rxbc))
  drest = (sum(abs(rxba))-sum(diag(rxba)))/(sum(abs(rxbc))-sum(diag(rxbc)))
  ddist = cor(t(dist))
  
  badCor = rxba>0.25 | rxbc>0.25
  diag(badCor)<-FALSE
  diff_expectedA = rxba[badCor] - ddist[badCor]
  diff_expectedB = rxbc[badCor] - ddist[badCor]
  
  diff_expected = sum(abs(diff_expectedA))/sum(abs(diff_expectedB))
  
  badCor1 = rxba>0.25 
  badCor2 = rxbc>0.25
  diag(badCor1)<-FALSE
  diag(badCor2)<-FALSE
  
  diff_expectedA = rxba[badCor] - ddist[badCor]
  diff_expectedB = rxbc[badCor] - ddist[badCor]
  
  message(paste(mix_set,mean(diag(rxba)),ddiag,drest,sum(badCor),diff_expected,sum(badCor2),sum(badCor1),sum(badCor1 & badCor2),sep='\t'))
  
  aranbef = cbind(as.matrix(diag(rxb)),rownames(rxb),'xCell (before)')
  rownames(aranbef) = str_c(rep('xCell (before)',dim(aranbef)[1]),' ',rownames(aranbef))
  
  aran = cbind(as.matrix(diag(rx)),rownames(rx),'xCell')
  rownames(aran) = str_c(rep('xCell',dim(aran)[1]),' ',rownames(aran))
  if (show.before==TRUE) {
    R = rbind(aran,aranbef)
  } else {
    R = aran
  }
  ciber = list()
  if (do.ciber == TRUE) {
    rc = ciberAnalysis(paste0(working.dir,'/simulations/CIBERSORT.',mix_set,'.txt'),dist,paste0(working.dir,'/simulations/plots/',mix_set,'_cibersort.pdf'))
    ciber = cbind(as.matrix(diag(rc$r)),rownames(rc$r),'Newman')
    rownames(ciber) = str_c(rep('Newman',dim(ciber)[1]),' ',rownames(ciber))
    A = intersect(rownames(rc$x),rownames(dist))
    z1 = cor(as.matrix(rc$x[A,]),as.matrix(rc$y[A,]))
    z2 = cor(M2[A,],dist[A,])
    blue <- rgb(0, 0, 1, alpha=0.5)
    red <- rgb(1, 0, 0, alpha=0.5)
    pdf(paste0(working.dir,'/simulations/plots/',mix_set,'_xcell_cibersort.pdf'),width=3,height=4)
    barplot(sort(diag(z2)),border=NA,xaxt='n',ylab="R (Pearson)",col=blue,space=0)
    barplot(sort(diag(z1)),border=NA,xaxt='n',col=red,space=0,add=TRUE)
    print(sprintf('N = %d, xCell = %.2f, CIBEROSRT = %.2f',length(A),median(diag(z2)),median(diag(z1))))
    dev.off()
    R = rbind(R,ciber)
    
  } 
  if (do.known == TRUE) {
    pdf(paste0(working.dir,'/simulations/plots/',mix_set,'_known.pdf'),width=4,height=8)
    rk = knownAnalysis(dist,paste0(working.dir,'/simulations/mix_',mix_set,'_known.txt'))
    dev.off()
    R = rbind(R,rk$out)
    
  }
  if (do.mcp == TRUE) {
    rmcp = mcp.analyze('bp_n12')
    mcp = cbind(rmcp,names(rmcp),'MCP')
    rownames(mcp) = str_c(rep('MCP',length(rmcp)),' ',names(rmcp))
    R = rbind(R,mcp)
  }
  
  if (do.ciber==TRUE && do.known==TRUE) {
    colnames(rc$r2) = paste('CIBERSORT',colnames(rc$r2))
    known_ciber_mat = cbind(rk$mat,rc$r2[rownames(rk$mat),])
    pdf(paste0(working.dir,'/simulations/plots/',mix_set,'_known_ciber.pdf'),width=4,height=8)
    mat = corrplot(as.matrix(known_ciber_mat),tl.cex = 0.7,tl.col = "black",col=colorRampPalette(rev(colors))(200))
    dev.off()
  } else {
    known_ciber_mat = rk$mat
  }
  
  if (do.ciber==TRUE || do.known==TRUE) {
    mat = compareCorrelations(R,c("xCell","xCell (before)","Bindea","Charoentong","Palmer","Rooney","Tirosh","Newman",'MCP'),NULL,paste0(working.dir,'/simulations/plots/',mix_set,'_compare.pdf'),8,3,drop.only1 = TRUE,transpose,ct_order=ct_order,font.size=font.size)
    
    
    pdf(paste0(working.dir,'/simulations/plots/page/',mix_set,'_page.pdf'),width=8,height=8)
    par(mfrow=c(2,2))
    corrplot(as.matrix(rxb[ct_order,ct_order]),tl.cex = font.size/10,tl.col = "black",col=colorRampPalette(rev(colors))(200),order=oclust)
    corrplot(as.matrix(rx[ct_order,ct_order]),tl.cex = font.size/10,tl.col = "black",col=colorRampPalette(rev(colors))(200),order=oclust)
    tl.col = c(rep('black',dim(rk$mat)[2]),rep('blue',dim(rc$r2)[2]))
    corrplot(as.matrix(known_ciber_mat),tl.cex = 0.7,tl.col = tl.col,col=colorRampPalette(rev(colors))(200))
    corrplot(as.matrix(mat),method="pie",na.label="-",tl.cex = font.size/10,cl.cex=font.size/10,number.cex=(font.size+3)/10,tl.col = "black",col=colorRampPalette(rev(colors))(200))
    dev.off()
    
    R 
  } else {
    rx
  }
  
}

mcp.analyze = function(mix_set) {
  ex=read.table(paste0(working.dir,'/simulations/mix_',mix_set,'_expr.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE)
  dist=read.table(paste0(working.dir,'/simulations/mix_',mix_set,'_dist.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE)
  mcp = MCPcounter.estimate(ex,'HUGO_symbols')
  rownames(mcp) = c('T-cells','CD8+ T-cells','Cytotoxic','NK cells','B-cells','Monocytes','DC','Neutrophils','Endothelial cells','Fibroblasts')
  A = intersect(rownames(mcp),rownames(dist))
  r = corr.test(t(mcp[A,]),t(dist[A,]))
  diag(r$r)
}

knownAnalysis = function(dist,fn) {
  colors <- brewer.pal(11, "RdBu")
  
  mk=read.table(fn,sep="\t",header=TRUE,row.names=1, as.is=TRUE)
  source = unlist(lapply(rownames(mk),function(x) {unlist(strsplit(x,'[_%]'))[1]}))
  known = unlist(lapply(rownames(mk),function(x) {unlist(strsplit(x,'[_%]'))[2]}))
  rownames(mk) = str_c(source," ",known)
  id = sort(known,index.return=TRUE)
  known = known[id$ix]
  source = source[id$ix]
  mk = mk[id$ix,]
  rows = known %in% rownames(dist)
  dist = dist[sort(rownames(dist)),]
  r=cor(t(mk),t(dist),method='pearson')
  out = data.frame()
  A = known %in% rownames(dist)[1]
  out = cbind(r[A,1,drop=FALSE],known[A],source[A])
  for (i in 2:length(rownames(dist))) {
    A = known %in% rownames(dist)[i]
    out = rbind(out,cbind(r[known %in% rownames(dist)[i],i,drop=FALSE],known[A],source[A]))
  }
  colnames(out) = c('r','cell_type','source')
  r2 = r[rows,]
  colors <- brewer.pal(11, "RdBu")
  
  mat = corrplot(as.matrix(t(r2)),tl.cex = 0.7,tl.col = "black",col=colorRampPalette(rev(colors))(200))
  res = list()
  res$mat = mat
  res$out = out
  res
}

fixScales = function(X,Y) {
  z = matrix(0,dim(X)[1])
  for (i in 1:dim(X)[1]) {
    x = as.matrix(X[i,])
    y = as.matrix(t(Y[i,]))
    fit <- lm(x ~ y)
    z[i] = coef(fit)[2]
  }
  rownames(z) = rownames(X)
  z
  
}

compareCorrelations = function(br,usr,uct,fn,w,h,drop.only1=TRUE,transpose=FALSE,font.size=7,ct_order=NULL) {
  colors <- brewer.pal(11, "RdBu")
  
  abr = aggregate(as.numeric(br[,1])~rownames(br),FUN=mean)
  nbr = sort(unique(str_c(br[,3],'_',br[,2])))
  ct = unlist(lapply(nbr,function(x) {unlist(strsplit(x,'[_%]'))[2]}))
  sr = tolower(unlist(lapply(nbr,function(x) {unlist(strsplit(x,'[_%]'))[1]})))
  #usr = unique(sr)
  if (is.null(uct)) {
    uct = sort(unique(ct))
  }
  mat = matrix(nrow=length(usr),ncol=length(uct))
  for (i in 1:length(usr)) {
    for (j in 1:length(uct)) {
      x = abr[sr==tolower(usr[i]) & ct==uct[j],2]
      if(length(x)==1)
        mat[i,j] = x
    }
  }
  colnames(mat) = uct
  rownames(mat) = usr
  n = apply(mat,2,function(x) sum(!is.na(x)))
  if(drop.only1==TRUE)
    mat = mat[,n>1]
  mat = mat[!apply(mat, 1, function(x){sum(!is.na(x))==0}),]
  pdf(fn,width=w,height=h)
  if (is.null(ct_order))
    ct_order = colnames(mat)
  if (transpose==TRUE)
    mat = t(mat[,ct_order])
  corrplot(as.matrix(mat),method="pie",na.label="-",tl.cex = font.size/10,cl.cex=font.size/10,number.cex=(font.size+3)/10,tl.col = "black",col=colorRampPalette(rev(colors))(200))
  dev.off()  
  mat
}

ciberAnalysis = function(fn,dist,fn_pdf,remove.dups=FALSE) {
  m=read.table(fn,sep="\t",header=TRUE,row.names=NULL, as.is=TRUE)
  if (remove.dups==TRUE) {
    m = createMatrixNoDups(m)
    B = intersect(rownames(m),colnames(dist))
    m = m[B,]
    dist = dist[,B]
  } else {
    m = m[,-1]
  }
  m = t(m[,1:22])
  rownames(m) = c("naive B-cells","Memory B-cells","Plasma cells","CD8+ T-cells","CD4+ naive T-cells","CD4+ memory T-cells rest","CD4+ memory T-cells act","Tfh","Tregs","Tgd","NK cells rest","NK cells act","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","DC rest","DC act","Mast cells rest","Mast cells act","Eosinophils","Neutrophils")
  cd4 = colSums(m[c("CD4+ naive T-cells","CD4+ memory T-cells rest","CD4+ memory T-cells act","Tregs"),])
  cd4mem = colSums(m[c("CD4+ memory T-cells rest","CD4+ memory T-cells act"),])
  nk = colSums(m[c("NK cells rest","NK cells act"),])
  dc = colSums(m[c("DC rest","DC act"),])
  bcell = colSums(m[c("naive B-cells","Memory B-cells","Plasma cells"),])
  mast = colSums(m[c("Mast cells rest","Mast cells act"),])
  macs = colSums(m[c("Macrophages M0","Macrophages M1","Macrophages M2"),])
  m = rbind(m,cd4,cd4mem,cd4mem,nk,dc,bcell,mast,macs)
  rownames(m)[23:30] = c("CD4+ T-cells","CD4+ Tcm","CD4+ Tem","NK cells","DC","B-cells","Mast cells","Macrophages")
  A = intersect(rownames(m),rownames(dist))
  
  M = m[A,]
  d = dist[A,]
  r=corr.test(t(M),t(d),method="pearson")
  r$r[is.na(r$r)] = 0
  pdf(fn_pdf,width=4,height=4)
  colors <- brewer.pal(11, "RdBu")
  corrplot(as.matrix(r$r),tl.cex = 0.7,tl.col = "black",col=colorRampPalette(rev(colors))(200))
  dev.off()
  r2=cor(M,d)
  barplot(sort(diag(r2)))
  #print(median(diag(r2)))
  r$x = M
  r$y = d
  
  r2 = cor(t(M),t(dist),method="pearson")
  r$r2 = t(r2)
  r
}


mixScatters = function(fn,m,dist,norm=FALSE,point.cex=0.25) {
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  corr_eqn <- function(x,y, digits = 2) {
    corr_coef <- round(cor(x, y,method="pearson"), digits = digits)
    paste("R = ", corr_coef)
  }
  
  names = rownames(dist)
  n = dim(dist)[1]
  nr = round(n/2+0.05)
  pdf(fn,width=nr,height=2)
  layout(matrix(c(1:(nr*2)),2,nr))
  par( mai = c(0, 0, 0, 0))
  for (i in 1:n) {
    y = t(as.matrix(dist[names[i],]))
    x = as.matrix(m[names[i],])
    
    if (norm==FALSE) {
      maxxy1 = max(max(x),max(y))
      maxxy2 = max(max(x),max(y))
    } else {
      x = x-min(x)
      y = y-min(y)
      maxxy1 = max(x)*1.1
      maxxy2 = max(y)*1.1
    }
    plot(x,y,pch=16,cex=point.cex,col="darkblue",yaxt='n',xaxt='n',xlim=c(0,maxxy1),ylim=c(0,maxxy2))
    #lines(lowess(x,y), col="red")
    abline(lm(as.vector(y)~as.vector(x)),col='red')
    r=sprintf('%.2f',cor(as.vector(x),as.vector(y),method="pearson"))
    #text(0.01, 0.95, labels=names[i], adj=c(0, .5),cex=0.65)
    text(0.005, maxxy2*0.95, labels=names[i], adj=c(0, .5),cex=0.85)
    text(maxxy1*0.6, maxxy2*0.02, labels=bquote(rho == ~ .(r)), adj=c(0, .5),cex=0.85)
    
  }
  dev.off()
}


library(psych)

analyze.facs = function(fcs,scores,sets,known=FALSE,spill=NULL,alpha=0.5,fn=NULL) {
  a = intersect(colnames(scores),colnames(fcs))
  rownames(fcs) = sets
  
  if(known==TRUE) {
    source = unlist(lapply(rownames(scores),function(x) {unlist(strsplit(x,'[_%]'))[1]}))
    ct = unlist(lapply(rownames(scores),function(x) {unlist(strsplit(x,'[_%]'))[2]}))
    #rownames(scores) = str_c(source," ",ct)
    y = scores[ct %in% sets,a]
    ct2 = ct[ct %in% sets]
    x = sapply(ct2,function(x) fcs[x,a])
    x <- matrix(unlist(x), ncol = dim(x)[2], byrow = FALSE)
    colnames(x) = rownames(y)
  } else {
    ct = intersect(sets,rownames(scores))
    y = scores[ct,a]
    x = t(fcs[ct,a])
  }
  
  if (is.null(spill)) { 
    s = as.matrix(y)
  } else {
    #spill$fv[,2]=1
    z = transformScores(as.matrix(y),spill$fv)
    s = spillOver(z,spill$K,alpha)
    #s = y
    A = intersect(colnames(x),rownames(s))
    x = x[,A]
    s=s[A,]
  }
  r = corr.test(x,t(s),method="pearson",adjust='none')
  dr = diag(r$r)
  names(dr) = colnames(r$r)
  print(as.matrix(dr))
  colors <- brewer.pal(11, "RdBu")
  corrplot(r$r,tl.cex=0.5,tl.col = "black",col=colorRampPalette(rev(colors))(200))
  
  if (!is.null(fn)) {
    pdf(paste0(fn,'.pdf'),width=3.2,height=4.7)
    corrplot(r$r,tl.cex=0.5,tl.col = "black",col=colorRampPalette(rev(colors))(200))
    dev.off()
    thres= sort(diag(r$r),decreasing = TRUE)[7]
    print(thres)
    A = diag(r$r)>thres
    mixScatters(paste0(fn,'_mix_scatters.pdf'),s[A,],as.matrix(t(x[,A])),norm=TRUE,point.cex=0.5)
    
  }
  r2 = cor(t(x),s,method="spearman")
  barplot(sort(diag(r2)))
  print(median(diag(r2)))
  r$x = t(x)
  r$y = s
  r
  
}

compare.facs = function(r,k,ciber,fn,drop.only1=TRUE) {
  aran = cbind(diag(r$r),colnames(r$r),'xCell',diag(r$p))
  rownames(aran) = str_c(rep('xCell',dim(aran)[1]),' ',aran[,2])
  cbr = cbind(diag(ciber$r),colnames(ciber$r),'Newman',diag(ciber$p))
  rownames(cbr) = str_c(rep('Newman',dim(cbr)[1]),' ',cbr[,2])
  source = unlist(lapply(colnames(k$r),function(x) {unlist(strsplit(x,'[_%]'))[1]}))
  ct = unlist(lapply(colnames(k$r),function(x) {unlist(strsplit(x,'[_%]'))[2]}))
  known = cbind(diag(k$r),ct,source,diag(k$p))
  rownames(known) = str_c(source," ",ct)
  abr = rbind(aran,known,cbr)
  
  usr = c('xCell','Bindea','Charoentong','Palmer','Rooney','Tirosh','Newman')
  ct = abr[,2]
  uct = unique(abr[,2])
  sr = tolower(abr[,3])
  
  mat.r= matrix(nrow=length(usr),ncol=length(uct))
  mat.p = mat.r
  for (i in 1:length(usr)) {
    for (j in 1:length(uct)) {
      x = abr[which(sr==tolower(usr[i]) & ct==uct[j]),]
      mat.r[i,j] = as.numeric(x[1])
      mat.p[i,j] = as.numeric(x[4])
      
    }
  }
  colnames(mat.r) = uct
  rownames(mat.r) = usr
  n = apply(mat.r,2,function(x) sum(!is.na(x)))
  if(drop.only1==TRUE) {
    mat.r = mat.r[,n>1]
    mat.r = mat.r[,!is.na(mat.r[1,])]
    mat.p = mat.p[,n>1]
    mat.p = mat.p[,!is.na(mat.r[1,])]
    avg.fcs=apply(fcs,1,function(x) mean(x,na.rm=TRUE))
    avg.fcs2 = which(colnames(mat.r) %in% rownames(fcs)[avg.fcs>0.01])
    mat.r = mat.r[,avg.fcs2]
    mat.p = mat.p[,avg.fcs2]
  }
  pdf(fn,width=5,height=3)
  corrplot(as.matrix(mat.r),p.mat=mat.p,sig.level=0.05,pch.col = "gray", pch.cex = 1,method="pie",na.label="-",tl.cex = 0.7,tl.col = "black",col=colorRampPalette(rev(colors))(200))
  dev.off()  
  mat.r  
  
}

fit.scores.plots = function(fv,families,working.dir,ctrl1,ctrl2,ctrl2_type=NULL,platform) {
  files1 = list.files(path=paste0(working.dir,'/mixtures/scores/'),pattern=paste0('_',ctrl1,'_'))
  files2 = list.files(path=paste0(working.dir,'/mixtures/scores/'),pattern=paste0('_',ctrl2,'_'))
  if (!is.null(ctrl2_type)) {
    for (i in 1:length(files1)) {
      ct = gsub('(.txt)','',strsplit(files1[i],'_')[[1]][3]);
      if(ctrl2_typefamilies[ct,1]==ctrl2_type) 
        files1[i] = 'NA'
    }
    for (i in 1:length(files2)) {
      ct = gsub('(.txt)','',strsplit(files2[i],'_')[[1]][3]);
      if(families[ct,1]!=ctrl2_type) 
        files2[i] = 'NA'
    }
  }
  files = c(files1,files2)
  files = files[files!="NA"]
  n = 32
  for (j in 1:3) {
    p = list()
    N = min(n,length(files)-n*(j-1))
    for (k in 1:N) {
      i = k+(j-1)*n
      ct = gsub('(.txt)','',strsplit(files[i],'_')[[1]][3]);
      sr = gsub('(.txt)','',strsplit(files[i],'_')[[1]][1]);
      mix <- as.matrix(read.table(paste0(working.dir,'/mixtures/scores/',files[i]), header=TRUE, sep="\t", row.names=1, as.is=TRUE))
      d = list()
      d$x = (mix[ct,seq(2,33)]-mix[ct,2])/5000
      d$y=seq(0.008,0.256,0.008)
      
      d$z = fv[ct,1]*d$x^fv[ct,2]
      d = as.data.frame(d)
      
      #message(sprintf('%d\t%s\t%s\t%f\t%f\t%f',i,sr,ct,mix[ct,2],coef(z)[1],coef(z)[2]))
      
      if (sr == 'blueprint') {
        sr = 'B'
      } else if (sr == 'encode') {
        sr 
      }
      if (ct == 'Class-switched memory B-cells') {
        ct = 'CS memory B-cells'
      }
      p[[k]] = ggplot(d, aes(y=y)) + geom_line(aes(x=x), colour="red") + geom_line(aes(x=z), colour="blue")+theme_bw()+geom_point(aes(y=y,x=x),size=0.25,col='black')+geom_point(aes(y=y,x=z),size=0.25,col='black')+annotate("text", hjust=0, x = 0.02, y = 0.205, size=2.5,label = sprintf('y==%.2f*x^%.2f',fv[ct,1],fv[ct,2]),parse=T)+annotate("text", hjust=0,x = 0.02, y = 0.24,size=2.5,label = paste0(ct,' (',toupper(substr(sr,1,1)),')'), fontface =2)
      #p[[k]] = ggplot(d, aes(y=y)) + geom_line(aes(x=x), colour="red") + geom_line(aes(x=z), colour="blue")+theme_bw()+geom_point(aes(y=y,x=x),size=0.25,col='black')+geom_point(aes(y=y,x=z),size=0.25,col='black')+annotate("text", hjust=0, x = 0.02, y = 0.205, size=2.5,label = sprintf('y==%.2f*x^%.2f', coef(z)[1],coef(z)[2]),parse=T)+annotate("text", hjust=0,x = 0.02, y = 0.24,size=2.5,label = paste0(ct,' (',toupper(sr),')'), fontface =2)
      p[[k]] = p[[k]] + theme(plot.margin=unit(c(0,0,0,0), "cm"),axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.background=element_blank())
    }
    pdf(paste0(working.dir,'/plots/supp3_',platform,'_',j,'.pdf'),width=6,height=ceil(N/4))
    do.call(grid.arrange,c(p,ncol=4))
    dev.off()
  }
}


test.mixture = function(mix.file,spill,alpha,cibersort=TRUE) {
  scores = as.matrix(read.table(paste0(mix.file,'_scores.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE))
  dist = as.matrix(read.table(paste0(mix.file,'_dist.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE))
  
  A = intersect(rownames(spill$fv),rownames(scores))
  scores = scores[A,]
  a = scores-apply(scores,1,min)
  calib = as.vector(spill$fv[rownames(scores),'calib'])
  fit_power = as.vector(spill$fv[rownames(scores),'fit'])
  transformed = ((a/16000)^fit_power)/(2*calib)
  transformed=t(t(transformed)/apply(transformed,2,sum))
  A = intersect(rownames(spill$K),rownames(transformed))
  K = spill$K[A,A]*alpha
  diag(K)<-1
  scores <- apply(transformed[A, ], 2, function(x) lsqlincon(K,x, lb = 0))
  rownames(scores) = A
  
  r = cor(t(scores),t(dist))
  if (cibersort==TRUE) {
    ciber = as.matrix(read.table(paste0(mix.file,'_CIBERSORT.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE))
    ciber = ciber[,1:(dim(ciber)[2]-3)]
    A = intersect(colnames(ciber),rownames(dist))
    rc  = cor(ciber,t(dist))
    rc[is.na(rc)] = 0
    A = intersect(rownames(r),rownames(rc))
    B = intersect(colnames(r),colnames(rc))
    res = list(dist=dist,xCell=scores,CIBERSORT=ciber,RX=r[A,B],RC=rc[A,B])
  } else {
    res = list(dist=dist,xCell=scores,RX=r)
  }
  
  colors <- brewer.pal(11, "RdBu")
  col=colorRampPalette(rev(colors))(200)
  
  pdf(paste0(mix.file,'.pdf'))
  A = intersect(colnames(res$RX),rownames(res$RX))
  corrplot(res$RX[A,A],title='xCell',col=col)
  if (cibersort==TRUE) {
    A = intersect(colnames(res$RC),rownames(res$RC))
    corrplot(res$RC[A,A],title='CIBERSORT',col=col)
    b = cbind(diag(res$RX[A,A]),diag(res$RC[A,A]))
    colnames(b) = c('xCell','CIBERSORT')
    corrplot(b,method='number',col=col)
  }
  
  corrplot(res$RX,title='xCell - All',col=col)
  
  if (cibersort==TRUE) {
    corrplot(res$RC,title='CIBERSORT - All',col=col)
  }
  mixScatters(res$xCell,dist)
  mtext("xCell", outer=TRUE,  cex=1, line=-0.5)
  if (cibersort==TRUE) {
    mixScatters(t(res$CIBERSORT),dist)
    mtext("CIBERSORT", outer=TRUE,  cex=1, line=-0.5)
  }
  dev.off()
  res
}

choose.types.to.use = function(samples,n,dependencies) {
  A = sample(length(samples),length(samples))
  types.to.use = c()
  i=1
  while (i <= length(A)) {
    id = which(dependencies$types==samples[A[i]])
    if (length(id)==1) {
      A = A[!(samples[A] %in% dependencies$dep[[id]])]
    }
    i=i+1
  }
  samples[A[1:min(n,length(A))]]
}