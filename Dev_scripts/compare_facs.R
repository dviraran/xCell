library(psych)
#spill = createSpilloverMatrix()

expr = read.table("~/Documents/Immport/GE_SDY311-2010.txt",sep="\t",header=TRUE,row.names=1, as.is=TRUE)
scores = rawEnrichmentAnalysis(as.matrix(expr),signatures,genes.use,paste0(working.dir,'/Immport/sdy311_scores.txt'))

expr = read.table("~/Documents/Immport/GE_SDY420.filtered.txt",sep="\t",header=TRUE,row.names=1, as.is=TRUE)
scores = rawEnrichmentAnalysis(as.matrix(expr),signatures,genes.use,paste0(working.dir,'/Immport/sdy420_scores.txt'))

sets = c('B-cells',	'CD16- monocytes',	'CD16+ monocytes',	'CD4+ T-cells',	'CD8+ T-cells',	'CD4+ Tcm',	'CD8+ Tcm',	'effector CD4+ T cells',	'effector CD8+ T cells',	'CD4+ Tem',	'CD8+ Tem',	'Tgd cells',	'lymphocytes',	'Memory B-cells',	'Monocytes',	'naive B-cells',	'CD4+ naive T-cells',	'CD8+ naive T-cells',	'NK cells',	'NKT',	'Plasma cells',	'T cells',	'pro B-cells',	'Tregs')

fcs = read.table("~/Documents/Immport/FCS_SDY311-2010.txt",sep="\t",header=TRUE,row.names=1, as.is=TRUE)
rownames(fcs) = sets
colors <- brewer.pal(11, "RdBu")

scores = read.table(paste0(working.dir,'/Immport/sdy311_scores.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE)
colnames(scores) = gsub("\\.1","",colnames(scores))
scores = aggregate(t(scores)~colnames(scores),FUN=mean)
rownames(scores) = scores[,1]
scores = scores[,-1]
scores311 =t(scores)

fcs= fcs[,-which(colnames(fcs) %in% c("SUB134240","SUB134283"))]
fcs311=fcs
r311=analyze.facs(fcs311,scores311,sets,FALSE,spill.array,alpha,paste0(working.dir,'/simulations/plots/sdy311'))

kscores = read.table("~/Documents/signatures/scores/sdy311_known.txt",sep="\t",header=TRUE,row.names=1, as.is=TRUE)
colnames(kscores) = gsub("\\.1","",colnames(kscores))
kscores = aggregate(t(kscores)~colnames(kscores),FUN=mean)
rownames(kscores) = kscores[,1]
kscores = kscores[,-1]
kscores311 =t(kscores)

k311=analyze.facs(fcs311,kscores311,sets,TRUE)
c311 = ciberAnalysis("~/Documents/signatures/cibersort/CIBERSORT.SDY311.txt",fcs311,paste0(working.dir,'/simulations/plots/sdy311_cibersort.pdf'),remove.dups=TRUE)

fcs = read.table("~/Documents/Immport/FCS_SDY420.filtered.txt",sep="\t",header=TRUE,row.names=1, as.is=TRUE)
rownames(fcs) = sets
fcs420=fcs

scores420 = read.table(paste0(working.dir,'/Immport/sdy420_scores.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE)
r420=analyze.facs(fcs420,scores420,sets,FALSE,spill.array,alpha,paste0(working.dir,'/simulations/plots/sdy420'))
kscores420 = read.table("~/Documents/signatures/scores/sdy420_known.txt",sep="\t",header=TRUE,row.names=1, as.is=TRUE)
k420=analyze.facs(fcs420,kscores420,sets,TRUE)

c420 = ciberAnalysis("~/Documents/signatures/cibersort/CIBERSORT.SDY420.txt",fcs420,paste0(working.dir,'/simulations/plots/sdy420_cibersort.pdf'))

m311 = compare.facs(r311,k311,c311,paste0(working.dir,'/simulations/plots/sdy311_compare.pdf'))
m420 = compare.facs(r420,k420,c420,paste0(working.dir,'/simulations/plots/sdy420_compare.pdf'))
