library(GSVA)
library(GSEABase)
library(pracma)
library(RColorBrewer)
library(pheatmap)

create.train.test.ref.sets = function(working.dir,data_sets,nsamples) {
  
  ref = lapply(data_sets, function (x) {
    message(x)
    ref.samples.file = file.path(working.dir,'reference_sets',paste0(x,'_samples.txt'))
    ref.expr.file = file.path(working.dir,'reference_sets',paste0(x,'_expr.txt'))
    
    expr = read.table(ref.expr.file,sep="\t",header=TRUE,row.names=1, as.is=TRUE)
    samples = read.table(ref.samples.file,sep="\t",header=FALSE,row.names=NULL, as.is=TRUE)
    A = !(samples[,2]=="NaN")
    expr = expr[,A]
    samples = samples[A,]
    
    types = sort(unique(samples[,2]))
    A = table(samples[,2])>nsamples
    train = c()
    test = c()
    for (i in types[A]) {
      ids = which(samples[,2]==i)
      sid = sample(ids,1)
      test = c(test,sid)
      train = c(train,setdiff(ids,sid))
    }
    for (i in types[!A]) {
      ids = which(samples[,2]==i)
      train = c(train,ids)
    }
    
    write.table(expr[,train],file.path(working.dir,'reference_sets',paste0(x,'_expr_train.txt')),sep="\t",row.names=TRUE,quote =FALSE,col.names = TRUE)
    write.table(expr[,test],file.path(working.dir,'reference_sets',paste0(x,'_expr_test.txt')),sep="\t",row.names=TRUE,quote =FALSE,col.names = TRUE)
    write.table(samples[train,],file.path(working.dir,'reference_sets',paste0(x,'_samples_train.txt')),sep="\t",row.names=FALSE,quote =FALSE,col.names = FALSE)
    write.table(samples[test,],file.path(working.dir,'reference_sets',paste0(x,'_samples_test.txt')),sep="\t",row.names=FALSE,quote =FALSE,col.names = FALSE)
    
  })
}

read.data.sets = function(working.dir,data_sets,train_test) {
  
  ref = lapply(data_sets, function (x) {
    ref.samples.file = file.path(working.dir,'reference_sets',paste0(x,'_samples_',train_test,'.txt'))
    ref.expr.file = file.path(working.dir,'reference_sets',paste0(x,'_expr_',train_test,'.txt'))
    
    expr = read.table(ref.expr.file,sep="\t",header=TRUE,row.names=1, as.is=TRUE)
    samples = read.table(ref.samples.file,sep="\t",header=FALSE,row.names=NULL, as.is=TRUE)
    ref = list(expr = expr, samples = samples, name=x)
  })
}


read.types.dependencies = function(working.dir,file.name) {
  con  <- file(file.path(working.dir,file.name), open = "r")
  out <- list()
  i=1
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    vec <- unlist(strsplit(oneLine, "\t"))
    out$types[i] = vec[1]
    n = max(which(!(vec=="")))
    if (n<2)
      out$dep[i] = list(types="")
    else
      out$dep[i] = list(types=vec[2:n])
    i=i+1
  } 
  close(con)
  out
}

create.signatures = function(ref,dependencies,genes.use,other.gene.set,no.filter.types,working.dir) {
  
  dir.create(file.path(working.dir, 'signatures'), showWarnings = FALSE)
  
  temp=lapply(ref,function (x) {
    message(paste('Dataset -',x$name))
    expr = x$expr[genes.use,] 
    if(x$type=='seq') {
      expr = log2(expr)
      expr[expr<0] = 0
    }
    
    A = !is.na(x$samples[,2])
    expr = expr[,A]
    samples = x$samples[A,2]
    probs = c(.1,.25,.33333333,.5,.6666666,.75,.9)
    diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2,3,4,5)
    
    types = unique(samples)
    
    message('get quantiles...')
    q = lapply(types,function(y) {
      A = samples==y
      if (sum(A)==1) {
        ex = cbind(expr[,samples==y],expr[,samples==y])
      } else {
        ex = expr[,A]
      }
      
      q = apply(ex,1,function(z) quantile(z,probs,na.rm=TRUE))
      
      # -- remove over expressed genes, but not for epithelial cells
      if (!(y %in% no.filter.types))
        q[,!(rownames(expr) %in% other.gene.set)]=0
      q
    })
    
    message('create all signatures...')
    
    for (i in 1:length(types)) {
      message(types[i])
      id = which(dependencies$types==types[i])
      if (length(id)==1) {
        A = !(types %in% dependencies$dep[[id]])
      } else {
        A = rep(TRUE,length(types))
      }
      ntypes = sum(A)
      
      for (diff in diff_vals) {
        for (j in 1:round(length(probs)/2+0.5)) {
          diff_genes = lapply(q[A],function(x) q[[i]][j,]>x[length(probs)-j,]+diff)
          output <- matrix(unlist(diff_genes), nrow = ntypes, byrow = TRUE)
          for (n in (ntypes-1):(ntypes-3)) {
            g = colSums(output)>=n
            if (sum(g)>7 & sum(g)<201)
              write.table(rownames(expr)[g],file=file.path(working.dir,'signatures',sprintf("%s%%%s%%j%g%%d%g%%n%g.txt",types[i],toupper(x$name),round(probs[j]*100),diff,ntypes-n)),sep="\t",row.names=FALSE,quote =FALSE,col.names = FALSE)
          }
        }
      }
    }
  })
  return(NULL)
}

score.ref.signatures = function(ref,genes.use,working.dir) {
  dir.create(file.path(working.dir, 'scores'), showWarnings = FALSE)
  egc = read.signatures.dir(file.path(working.dir,'signatures'))
  lapply(ref, function(x) {
    scores = gsva(as.matrix(x$expr[genes.use,]),egc,method='ssgsea',ssgsea.norm=FALSE)
    write.table(scores,file=file.path(working.dir,"scores",paste(x$name,"_ssgsea.txt")),sep="\t",row.names=TRUE,quote =FALSE,col.names = NA)
  })
}

choose.best.signatures = function(ref,nsigs,dependencies,top.diff.percent=1,working.dir) {
  dir.create(file.path(working.dir, 'best_signatures'), showWarnings = FALSE)
  lapply(ref, function(x) {
    
    samples = x$samples[,2]
    types = unique(samples)
    message(x$name)
    scores = read.table(file.path(working.dir,"scores",paste(x$name,"_ssgsea.txt")),sep="\t",header=TRUE,row.names=1, as.is=TRUE)
    
    signatures_split = unlist(strsplit(rownames(scores),'%'))
    cell_types = signatures_split[seq(1,length(signatures_split),5)]
    sources = signatures_split[seq(2,length(signatures_split),5)]
    
    message(ref$name)
    for (i in types) {
      id = which(dependencies$types==i)
      if (length(id)==1)
        opt_in = types[!(types %in% dependencies$dep[[id]])]
      else
        opt_in = types
      
      message(i)
      sig.use = cell_types==i
      if (length(unique(sources[sig.use]))>1) {
        sig.use = cell_types==i & !(sources==toupper(x$name))
      }
      if (sum(sig.use)>0) {
        ttest = apply(scores[sig.use,],1,function(y) {
          others = samples %in% setdiff(opt_in,i)
          n = round(sum(others)*top.diff.percent)
          s = sort(y[others],decreasing = TRUE)
          A = samples==i
          if (sum(A)==1) {
            B = others & as.vector(y>as.numeric(s[n]))
            z = (y-min(y[B]))/(max(y[B])-min(y[B]))
            as.numeric(z[A]-max(z[B]))
          } else {
            t.test(y[A],y[others & as.vector(y>as.numeric(s[n]))],alternative='greater')$statistic
          }
        })
        ttest.ordered = sort(ttest,decreasing = TRUE)
        for (j in 1:nsigs){
          fn = paste0(i,'%',toupper(x$name),'%',j,'.txt')
          file.copy(from=file.path(working.dir, 'signatures',names(ttest.ordered[j])), to=file.path(working.dir, 'best_signatures',fn))
        }
      }
      
    }
  })
}  

read.signatures.dir = function(path) {
  sig = list()
  files = list.files(path=path)
  for (i in 1:length(files)) {
    sig[i] = GeneSet(scan(paste0(path,'/',files[i]), what="", sep="\n"), setName=files[i])
  }
  signatures <- GeneSetCollection(sig)
  toGmt(signatures,'~/Documents/signatures/Mouse/signatures.txt')
  signatures
}

create.mix = function(expr,types,types.to.use,signatures,genes.use,nrep,noise.add=0,out.name) {
  types.to.use = intersect(types.to.use,types)
  n.oftypes = length(types.to.use)
  dist = rand(n.oftypes,nrep)
  dist = as.matrix(apply(dist,2,function(x) x/sum(x)))
  mix = matrix(0,nrep,dim(expr)[1])
  for (i in 1:n.oftypes) {
    message(types.to.use[i])
    A = which(types==types.to.use[i])
    if(length(A)>1) {
      ex = expr[,sample(A,nrep,replace=TRUE)]
    } else {
      ex = expr[,rep(A,nrep)]
    }
    noise = runif(dim(ex)[1],-noise.add,noise.add)
    ex = ex+ex*noise
    mix = mix + t(ex)*dist[i,]
    
  }
  rownames(dist) = types.to.use
  mix = t(mix)
  colnames(mix) = colnames(dist)
  write.table(mix, file = paste0(out.name,'_n',n.oftypes,'_expr.txt'), sep = "\t", col.names=NA,quote = FALSE)
  write.table(dist, file = paste0(out.name,'_n',n.oftypes,'_dist.txt'), sep = "\t", col.names=NA,quote = FALSE)
  
  scores = rawEnrichmentAnalysis(as.matrix(mix),signatures,genes.use,paste0(out.name,'_n',n.oftypes,'_scores.txt'))
  
  res = list()
  res$mix = mix
  res$dist = dist
  res$scores = scores
  res
}

create.onecell.mix = function(ref,genes.use,control.type,nrep,working.dir) {
  dir.create(file.path(working.dir, 'mixtures'), showWarnings = FALSE)
  dir.create(file.path(working.dir, 'mixtures/expr'), showWarnings = FALSE)
  
  lapply(ref, function(x) {
    
    samples = x$samples[,2]
    types = unique(samples)
    CTRL =  as.matrix(apply(as.matrix(x$expr[genes.use,samples==control.type]),1,mean))
    for (i in 1:length(types)) {
      ct = types[i]
      message(ct)
      d = t(as.matrix(seq(0,1,by=1/nrep)))
      mix = matrix(0,length(d),dim(x$expr)[1])
      A = which(samples==ct)
      if (length(A)==1)
        expr = as.matrix(x$expr[genes.use,A])
      else
        expr = as.matrix(apply(x$expr[genes.use,A],1,mean))
      mix = expr%*%d+CTRL%*%(1-d)
      write.table(mix, file = file.path(working.dir,'/mixtures/expr/',paste0(x$name,'_',control.type,'_',ct,'.txt')), sep = "\t", col.names=NA,quote = FALSE)
    }
  })
}

score.onecell.mix = function(signatures,genes.use,working.dir,rescore.all=TRUE) {
  dir.create(file.path(working.dir, 'mixtures/scores'), showWarnings = FALSE)
  
  files = list.files(path=paste0(working.dir,'/mixtures/expr/'),pattern='txt')
  for (i in files) {
    message(i)
    out_file = file.path(working.dir,'/mixtures/scores/',i)
    if(rescore.all==FALSE) {
      if (file.exists(out_file))
        next;
    } 
    mix=read.table(file.path(working.dir,'mixtures/expr',i),sep="\t",header=TRUE,row.names=1, as.is=TRUE)
    scores = rawEnrichmentAnalysis(as.matrix(mix[,1:35]),signatures,genes.use,out_file)
  }
}

create.ref.mix = function(ref,genes.use,control.type,signatures,percent,working.dir) {
  dir.create(file.path(working.dir, 'mixtures'), showWarnings = FALSE)
  
  lapply(ref, function(x) {
    samples = x$samples[,2]
    types = unique(samples)
    
    CTRL =  as.matrix(apply(as.matrix(x$expr[genes.use,which(samples==control.type)]),1,mean))
    
    expr = apply(as.matrix(types),1,function(y) {apply(as.matrix(x$expr[genes.use,samples==y]),1,mean)})
    colnames(expr) = types
    rownames(expr) = genes.use
    mix = expr*percent + (CTRL%*%matrix(1,1,length(types)))*(1-percent)
    
    fn = file.path(working.dir,'mixtures',paste0(x$name,'_',control.type,'_',round(100*percent),'p.txt'))
    
    scores = rawEnrichmentAnalysis(as.matrix(mix),signatures,genes.use,fn)
    as.matrix(scores)
  })
}

learn.transform.parameters = function(data_set,scores_path,ctrl_type,nrep,percent) {
  files = c()
  for (i in data_sets)
    files = c(files,list.files(path=scores_path,pattern=paste0(i,'_',ctrl_type,'_')))
  
  fit_table = list()
  
  d = t(as.matrix(seq(0,1,by=1/nrep)))
  pid = max.col(1-abs((d-percent)))+1
  
  for (i in 1:length(files)) {
    ct = gsub('(.txt)','',strsplit(files[i],'_')[[1]][3]);
    sr = gsub('(.txt)','',strsplit(files[i],'_')[[1]][1]);
    if (!(ct==ctrl_type)) {
      mix = as.matrix(read.table(file.path(scores_path,files[i]),sep="\t",header=TRUE,row.names=1, as.is=TRUE))
      df = data.frame(x=(mix[ct,seq(2,pid)]-mix[ct,2])/5000,y=d[2:pid])
      
      if (tail(df$x,1)<df$x[1]) {
        fit_table$value1[i] = 1
        fit_table$value2[i] = 1
      } else {
        z = try(nls(y ~ a * x^b, data = df, start = list(a=1, b=1),control = list(maxiter = 500)))
        fit_table$name[i] = ct
        if (class(z) == "try-error") {
          fit_table$value1[i] = NA
          fit_table$value2[i] = NA
          fit_table$calib[i] = NA
        } else {
          fit_table$value1[i] = coef(z)[1]
          fit_table$value2[i] = coef(z)[2]
          fit_table$calib[i] = coef(lm(df$x^coef(z)[2]~df$y))[2]
          message(sprintf('%d\t%s\t%s\t%f\t%f\t%f\t%f',i,sr,ct,mix[ct,2],coef(z)[1],coef(z)[2],fit_table$calib[i]))
          
        }
      }
      
    }
  }
  
  fv = aggregate(cbind(fit_table$value1,fit_table$value2,fit_table$calib),by=list(fit_table$name),FUN=mean)
  rownames(fv) = fv[,1]
  fv = fv[,-1]
  fv
}

transform.and.scale = function(mat,fv) {
  md <- apply(mat, 1, median)
  A <- rownames(mat)
  mat <- (mat - md)/5000
  mat[mat < 0] <- 0
  mat <- (mat[A,]^fv[A,2])/fv[A,3]
}

combine.ref.mix = function(fv,ref1,ref2,families,ref2_types) {
  temp = t(cbind(transform.and.scale(ref1[[1]],fv),transform.and.scale(ref1[[2]],fv),transform.and.scale(ref1[[3]],fv)))
  temp=aggregate(temp,by=list(rownames(temp)),FUN=mean)
  rownames(temp) = temp[,1]
  temp = temp[,-1]
  ref1 = t(temp)
  
  temp = t(cbind(transform.and.scale(ref2[[1]],fv),transform.and.scale(ref2[[2]],fv),transform.and.scale(ref2[[3]],fv)))
  temp=aggregate(temp,by=list(rownames(temp)),FUN=mean)
  rownames(temp) = temp[,1]
  temp = temp[,-1]
  ref2 = t(temp)
  
  ref.mix = ref1
  
  temp = intersect(rownames(ref.mix),rownames(families)[families[,1] %in% ref2_types])
  ref.mix[temp,] = ref2[temp,]
  ref.mix[,temp] = ref2[,temp]
  ref.mix
}

combine.ref.mix.array = function(fv.array,ref1,ref2,ref2_types) {
  temp = t(cbind(transform.and.scale(ref1[[1]],fv.array),transform.and.scale(ref1[[2]],fv.array)))
  temp=aggregate(temp,by=list(rownames(temp)),FUN=mean)
  rownames(temp) = temp[,1]
  temp = temp[,-1]
  ref1 = t(temp)
  
  temp = t(cbind(transform.and.scale(ref2[[1]],fv.array),transform.and.scale(ref2[[2]],fv.array),transform.and.scale(ref2[[3]],fv.array)))
  temp=aggregate(temp,by=list(rownames(temp)),FUN=mean)
  rownames(temp) = temp[,1]
  temp = temp[,-1]
  ref2 = t(temp)
  
  ref.array.mix = rbind(ref1[!(rownames(ref1) %in% ref2_types),],ref2[ref2_types,colnames(ref1)])
  ref.array.mix = cbind(ref.array.mix[,!(colnames(ref1) %in% ref2_types)],ref2[rownames(ref.array.mix),ref2_types])
  ref.array.mix
}

create.spillover.matrix = function(ref.mix,fv,families,working.dir) {
  K = ref.mix/diag(ref.mix)
  fv[rownames(K),3] = 2*diag(ref.mix)*fv[rownames(K),3]
  K[K>0.5] = 0.5
  families = families[rownames(K),]
  for (i in 1:dim(K)[1]) {
    id = which(dependencies$types==rownames(K)[i])
    if (length(id)==1) {
      A = rownames(K) %in% dependencies$dep[[id]]
      K[i,A] = 0
      K[A,i] = 0
    }
   if (families[i,2]=="Parent") {
      K[i,which(families[,2]!="Parent")]=0
    } else {
      K[i,rownames(K) == families[rownames(K)[i],2]]=0
      K[i,intersect(which(!(families[,2] %in% families[rownames(K)[i],2])),which(families[,2]!="Parent"))]=0
    }
  }

  diag(K) <- 1
  
  pheatmap(K,cluster_rows=FALSE,cluster_cols=FALSE)
  
  spill = list(K = K,fv = fv)
}
