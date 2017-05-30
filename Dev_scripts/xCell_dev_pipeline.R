working.dir = "~/Documents/xCell"
source(paste0(working.dir,'/xCell.R'))
source(paste0(working.dir,'/xCell_dev_functions.R'))

# -- read reference data sets -- #
data_sets = c('fantom','blueprint','encode','iris','novershtern','hpca')
data_types = c('seq','seq','seq','array','array','array')

families=read.table(file.path(working.dir,'cells_families.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE)
dependencies = read.types.dependencies(working.dir,'types_dependecies.txt')

# --- create train and test sets --- #
create.train.test.ref.sets(working.dir,data_sets,nsamples=4)

# read reference data sets
ref = read.data.sets(working.dir,data_sets,'train')
for(i in 1:length(ref)) {
  ref[[i]]$type=  data_types[i]
}

# get common genes
genes.use = rownames(ref[[1]]$expr)
for (i in 2:length(ref)) {
  genes.use = intersect(genes.use,rownames(ref[[i]]$expr))
}
length(genes.use)

# -- create all signatures and choose best signatures -- #
not_over_expressed = unlist(read.delim(file.path(working.dir,'not_overexpressed_genes.txt'),sep='\n'))
create.signatures(ref,dependencies,genes.use,not_over_expressed,rownames(families)[families[,3]=="Epithelial"],working.dir)
score.ref.signatures(ref,genes.use,working.dir)
choose.best.signatures(ref,3,dependencies,1,working.dir)

# save signatures files as Gmt format
signatures = read.signatures.dir(file.path(working.dir,'best_signatures'))
toGmt(signatures,file.path(working.dir,'xCell_signatures.txt'))
signatures = getGmt(file.path(working.dir,'xCell_signatures.txt'))

# -- trasnformation RNA-seq -- #
nrep = 125

# -- create synthetic mixtures for learning transformation parameters and score them -- #
nrep=125
create.onecell.mix(ref[data_types=="seq"],genes.use,'MPP',nrep,working.dir)
create.onecell.mix(ref[data_types=="seq"],genes.use,'Endothelial cells',nrep,working.dir)
create.onecell.mix(ref[data_types=="array"],genes.use,'Erythrocytes',nrep,working.dir)
create.onecell.mix(ref[data_types=="array"],genes.use,'Monocytes',nrep,working.dir)
score.onecell.mix(signatures,genes.use,working.dir,FALSE)

# -- learn transformation parameteres - seq -- #
fv.mpp = learn.transform.parameters(data_sets[data_types=="seq"],paste0(working.dir,'/mixtures/scores'),'MPP',125,0.25)
fv.endo = learn.transform.parameters(data_sets[data_types=="seq"],paste0(working.dir,'/mixtures/scores'),'Endothelial cells',125,0.25)

fv = fv.mpp[intersect(rownames(fv.mpp),rownames(families)[!(families[,1]=="HSC")]),]
fv = rbind(fv,fv.endo[intersect(rownames(fv.endo),rownames(families)[families[,1]=="HSC"]),])
temp = setdiff(rownames(families),rownames(fv))
x = t(matrix(c(mean(fv[,1]),mean(fv[,2]),mean(fv[,3])),3,length(temp)))
rownames(x) = temp
fv = rbind(fv,x)

# -- learn transformation parameteres - array -- #
fv.eryt = learn.transform.parameters(data_sets[5:6],paste0(working.dir,'/mixtures/scores'),'Erythrocytes',125,0.25)
fv.mono = learn.transform.parameters(data_sets[4:6],paste0(working.dir,'/mixtures/scores'),'Monocytes',125,0.25)

fv.array = rbind(fv.eryt,fv.mono)
temp = setdiff(rownames(families),rownames(fv.array))
x = t(matrix(c(mean(fv.array[,1]),mean(fv.array[,2]),mean(fv.array[,3])),3,length(temp)))
rownames(x) = temp
fv.array = rbind(fv.array,x)

# -- create synthetic references (25%) -- #
ref.mpp = create.ref.mix(ref[data_types=="seq"],genes.use,'MPP',signatures,0.25,working.dir)
ref.endo = create.ref.mix(ref[data_types=="seq"],genes.use,'Endothelial cells',signatures,0.25,working.dir)
ref.eryt = create.ref.mix(ref[5:6],genes.use,'Erythrocytes',signatures,0.25,working.dir)
ref.mono = create.ref.mix(ref[data_types=="array"],genes.use,'Monocytes',signatures,0.25,working.dir)

# -- combine reference matrices -- #
ref.mix = combine.ref.mix(fv,ref.mpp,ref.endo,families,c('HSC'))
ref.array.mix = combine.ref.mix.array(fv.array,ref.eryt,ref.mono,c('Erythrocytes','Th1 cells','Th2 cells','CD4+ memory T-cells','MEP','Platelets'))

# add missing cell types
temp = ref.array.mix[rownames(ref.mix),setdiff(colnames(ref.array.mix),colnames(ref.mix))]
ref.mix = cbind(ref.mix,temp)
ref.mix = ref.mix[rownames(ref.mix),row.names(ref.mix)]

temp = ref.mix[rownames(ref.array.mix),setdiff(colnames(ref.mix),colnames(ref.array.mix))]
ref.array.mix = cbind(ref.array.mix,temp)
ref.array.mix = ref.array.mix[rownames(ref.array.mix),row.names(ref.array.mix)]

# -- create spillover matrices and save data -- #
spill = create.spillover.matrix(ref.mix,fv,families,working.dir)
spill.array = create.spillover.matrix(ref.array.mix,fv.array,families,working.dir)

xCell.data = list(spill=spill,spill.array=spill.array,signatures=signatures,genes=genes.use)
save(xCell.data,file=file.path(working.dir,'xCell.data.Rdata'))
save(ref.mpp,ref.endo,ref.mix,fv.mpp,fv.endo,fv,ref.eryt,ref.mono,ref.array.mix,fv.eryt,fv.mono,fv.array,file=file.path(working.dir,'dev.data.Rdata'))
