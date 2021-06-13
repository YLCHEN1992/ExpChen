# TCGA simple base codes
# Version 2021.6.13
# Including  

# need following packages
library(pheatmap)
library(R.utils)
library(ggplot2) 
library(PerformanceAnalytics)

# public functions 

V=function(x){
if(length(which(x==""))>0){
x=x[-which(x=="")]}
x}

hotmap=function(ssdata,rowss,signs="热图",Larg=600){
adr=as.character(getwd())
cat("开始绘制热图\n")
library(pheatmap)
library(ggplot2) 
ann_colors = list(SampleType= c(Primary_Tumor="#D95F02", Solid_Tissue_Normal="#1B9E77"))
rrPlot=pheatmap(ssdata,color = colorRampPalette(c("green", "grey2", "red"))(500),scale = "row",
cluster_row = FALSE,annotation_row=rowss,border=FALSE,annotation_colors = ann_colors,show_rownames=F )
if (file.exists("./Graph_TCGA")==TRUE){cat("阁下目标文件夹 Graph_TCGA 已存在\n")}else{
dir.create("./Graph_TCGA", recursive=TRUE)
cat("目标文件夹 Graph_TCGA 已为阁下创建\n")}
setwd("./Graph_TCGA")
gNAME=paste("阁下 Graph_TCGA",signs,"已绘制完成",gsub(":","_",Sys.time()),".tiff")
ggsave(filename=gNAME,rrPlot,dpi=Larg,width=8,height=8)
cat("阁下 Graph_TCGA",signs,"已绘制完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(adr)}

volcane=function(ssdata,signs="火山图",Larg=600){
adr=as.character(getwd())
cat("开始绘制火山图\n")
plotdata=subset(ssdata,as.character(p_value)!="NaN" & as.character(LOG2)!="Inf")
Significant=c()
for(o in 1:nrow(plotdata)){
slog=0
spv=0
slog=as.numeric(plotdata$LOG2[o])
spv=as.numeric(plotdata$p_value[o])
if(slog>=1 & spv<=0.05){
Significant[o]="Up"}else if(
slog<=-1 & spv<=0.05){
Significant[o]="Down"}else{
Significant[o]="No"}}
plotdata$Significant=Significant
library(ggplot2) 
Plot=ggplot(plotdata,aes(LOG2,-1*log10(p_value)))
addcolor=Plot + geom_point(aes(color =Significant))
rPlot=addcolor+ labs(title="Volcanoplot",x=expression(Log2),y=expression(-log10(p_value)))
rrPlot=rPlot+geom_hline(yintercept=1.3)+geom_vline(xintercept=c(-1,1))
if (file.exists("./Graph_TCGA")==TRUE){cat("阁下目标文件夹 Graph_TCGA 已存在\n")}else{
dir.create("./Graph_TCGA", recursive=TRUE)
cat("目标文件夹 Graph_TCGA 已为阁下创建\n")}
setwd("./Graph_TCGA")
gNAME=paste("阁下 Graph_TCGA",signs,"已绘制完成",gsub(":","_",Sys.time()),".tiff")
ggsave(filename=gNAME,rrPlot,dpi=Larg,width=8,height=8)
cat("阁下 Graph_TCGA",signs,"已绘制完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(adr)}

ExpChen=function(gene,mir,ns,control="OOO" ,g="O",m="O",Larg=600){
cat("欢迎使用TCGA表达分析程序\n★ExpChen★\n")
cat("Current Version: 2021 \n")
cat("程序开始运行\n")
adr=as.character(getwd())
ns=read.csv(deparse(substitute(ns)))
if(as.character(g)==""){cat("Gene不导入\n")}else{
cat("开始写入gene表单\n")
xsheet=read.csv(deparse(substitute(gene)))
xfiles=as.character(xsheet$File.ID)
xfilenames=as.character(xsheet$File.Name)
xproject=as.character(xsheet$Project.ID)
xcase=as.character(xsheet$Case.ID)
xsample=as.character(xsheet$Sample.ID)
xsampletype=as.character(xsheet$Sample.Type)
xdatatype=as.character(xsheet$Data.Type)
xnumsheet=nrow(xsheet)
cat("正在导入gene数据\n")
library(R.utils)
for(i in 1:xnumsheet){
T=Sys.time()
novel=""
novel=paste("./gene_gdc_download/",xfiles[i],sep="")
setwd(novel)
assign(paste("D",i,sep=""),read.table(gzfile(xfilenames[i]),
header=FALSE,col.names=c("ENSEMBL","FPKM")))
setwd(adr)
TP=Sys.time()
DT=as.numeric(TP-T)
cat("Gene已读取",round(i*100/xnumsheet,3),
"%\n 当前速度为",as.character(round(as.numeric(DT),3)),"秒/个\n",
"大约还需要",as.character(round(as.numeric(DT),3)*(xnumsheet-i)),"秒\n")} #over for i
if(substring(control,1,1)=="O"|substring(control,2,2)=="O"){
cat("Gene进行Primary Tumor和Solid Tissue Normal分类\n")
gtdatas=matrix(c(1:(nrow(subset(xsheet,Sample.Type=="Primary Tumor"))+1)),ncol=1)
gndatas=matrix(c(1:(nrow(subset(xsheet,Sample.Type=="Solid Tissue Normal"))+1)),ncol=1)
DDD=as.character(ns$ENSEMBL)
DDD=V(DDD)
for(j in 1:length(DDD)){
T=Sys.time()
nid=""
nid=as.character(DDD[j])
totalt=c()
totaln=c()
for(z in 1:xnumsheet){
if(xsheet$Sample.Type[z]=="Primary Tumor"){
csf=0
csf=subset(get(paste("D",z,sep="")),substring(ENSEMBL,1,15)==nid)$FPKM
totalt=c(totalt,csf)}else if(xsheet$Sample.Type[z]=="Solid Tissue Normal"){
csf=0
csf=subset(get(paste("D",z,sep="")),substring(ENSEMBL,1,15)==nid)$FPKM
totaln=c(totaln,csf)}}#over for z & if
tdata=matrix()
tdata=matrix(c(nid,totalt),ncol=1)
gtdatas=cbind(gtdatas,tdata)
ndata=matrix()
ndata=matrix(c(nid,totaln),ncol=1)
gndatas=cbind(gndatas,ndata)
TP=Sys.time()
DT=as.numeric(TP-T)
cat(" Gene【分类】完成",round(j*100/length(DDD),3),
"%\n 当前速度为",as.character(round(as.numeric(DT),3)),"秒/个\n",
"大约还需要",as.character(round(as.numeric(DT),3)*(length(DDD)-j)),"秒\n")}#over for j
gtdatas=gtdatas[,-1]
gndatas=gndatas[,-1]
cat("分类完成\n开始统计数据\n")
mirnaid=c()
pvalues=c()
TMFC=c()
NMFC=c()
TSDFC=c()
NSDFC=c()
for(i in 1:ncol(gtdatas)){
T=Sys.time()
mirnaid=c(mirnaid,as.character(gtdatas[1,i]))
tcvs=c()
ncvs=c()
tcvs=as.numeric(gtdatas[2:nrow(gtdatas),i])
ncvs=as.numeric(gndatas[2:nrow(gndatas),which(gndatas[1,]==mirnaid[i])])
TMFC=c(TMFC,mean(tcvs))
TSDFC=c(TSDFC,sd(tcvs))
NMFC=c(NMFC,mean(ncvs))
NSDFC=c(NSDFC,sd(ncvs))
pp=t.test(tcvs,ncvs)$p.value
pvalues=c(pvalues,pp)
TP=Sys.time()
DT=as.numeric(TP-T)
cat(" Gene【统计】完成",as.character(round(i*100/ncol(gtdatas),3)),
"%\n 当前速度为",as.character(round(as.numeric(DT),3)),"秒/个\n",
"大约还需要",as.character(round(as.numeric(DT),3)*(ncol(gtdatas)-i)),"秒\n")}#over for i
cat("统计完成\n开始制表\n")
gssdata=data.frame(Gene_ID=mirnaid,
TumorMFC=as.numeric(TMFC),
NormalMFC=as.numeric(NMFC),
LOG2=log(TMFC/NMFC,2),
p_value=pvalues,
TumorSDFC=as.numeric(TSDFC),
NormalSDFC=as.numeric(NSDFC))
NAME=paste("阁下 TCGA-Gene 分析结果数据",gsub(":","_",Sys.time()),".csv")
if (file.exists("./TCGA")==TRUE){cat("阁下目标文件夹 TCGA 已存在\n")}else{
dir.create("./TCGA", recursive=TRUE)
cat("目标文件夹 TCGA 已为阁下创建\n")}
setwd("./TCGA")
write.csv(gssdata,NAME,row.names=FALSE)
cat(NAME,"已统计完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(adr)}#over for 1O2O
if(substring(control,1,1)=="O"){
volcane(gssdata,signs="Gene火山图",Larg=Larg)}else{cat("火山图不绘制\n")}
if(substring(control,2,2)=="O"){
GEMB=as.character(ns$ENSEMBL)
GEMB=V(GEMB)
xxrownames=c()
for(i in 1:length(GEMB)){
xxrownames=c(xxrownames,as.character(ns$SYMBOL[which(ns$ENSEMBL==as.character(GEMB[i]))]))}
colnames(gtdatas)=as.character(gtdatas[1,])
colnames(gndatas)=as.character(gndatas[1,])
gtdatas=gtdatas[-1,]
gndatas=gndatas[-1,]
gheatmap=data.frame(rbind(gtdatas,gndatas),stringsAsFactors = F)
gheatmap=as.data.frame(lapply(gheatmap,as.numeric))
rownames(gheatmap)=c(paste("id",c(1:nrow(gheatmap)),sep=""))
trow_ann=data.frame(SampleType=c(rep("Primary_Tumor",nrow(gtdatas))))
nrow_ann=data.frame(SampleType=c(rep("Solid_Tissue_Normal",nrow(gndatas))))
grow_annotation=data.frame(rbind(trow_ann,nrow_ann),stringsAsFactors = F)
rownames(grow_annotation)=c(paste("id",c(1:nrow(grow_annotation)),sep=""))
colnames(gheatmap)=xxrownames
for(n in 1:ncol(gheatmap)){
gheatmap[,n]=gheatmap[,n]/mean(as.numeric(gndatas[,n]))}
hotmap(gheatmap,grow_annotation,signs="Gene热图",Larg=Larg)}else{cat("热图不绘制\n")}
}#over if gene
if(as.character(m)==""){cat("miRNA不导入\n")}else{
cat("开始写入miRNA表单\n")
ysheet=read.csv(deparse(substitute(mir)))
yfiles=as.character(ysheet$File.ID)
yfilenames=as.character(ysheet$File.Name)
yproject=as.character(ysheet$Project.ID)
ycase=as.character(ysheet$Case.ID)
ysample=as.character(ysheet$Sample.ID)
ysampletype=as.character(ysheet$Sample.Type)
ydatatype=as.character(ysheet$Data.Type)
ynumsheet=nrow(ysheet)
cat("正在导入miRNA数据\n")
for(j in 1:ynumsheet){
T=Sys.time()
novel=""
novel=paste("./mirna_gdc_download/",yfiles[j],sep="")
setwd(novel)
assign(paste("S",j,sep=""),read.table(yfilenames[j],header=TRUE))
setwd(adr)
TP=Sys.time()
DT=as.numeric(TP-T)
cat("miRNA已读取",round(j*100/ynumsheet,3),
"%\n 当前速度为",as.character(round(as.numeric(DT),3)),"秒/个\n",
"大约还需要",as.character(round(as.numeric(DT),3)*(ynumsheet-j)),"秒\n")}#over j
if(substring(control,1,1)=="O"|substring(control,2,2)=="O"){
cat("miRNA进行Primary Tumor和Solid Tissue Normal分类\n")
tdatas=matrix(c(1:(nrow(subset(ysheet,Sample.Type=="Primary Tumor"))+1)),ncol=1)
ndatas=matrix(c(1:(nrow(subset(ysheet,Sample.Type=="Solid Tissue Normal"))+1)),ncol=1)
SSS=as.character(ns$IDmirna)
SSS=V(SSS)
for(j in 1:length(SSS)){
T=Sys.time()
nid=""
nid=as.character(SSS[j])
totalt=c()
totaln=c()
for(z in 1:ynumsheet){
if(ysheet$Sample.Type[z]=="Primary Tumor"){
csf=0
csf=subset(get(paste("S",z,sep="")),miRNA_ID==nid)$reads_per_million_miRNA_mapped
totalt=c(totalt,csf)}else if(ysheet$Sample.Type[z]=="Solid Tissue Normal"){
csf=0
csf=subset(get(paste("S",z,sep="")),miRNA_ID==nid)$reads_per_million_miRNA_mapped
totaln=c(totaln,csf)}}#over if for z
tdata=matrix()
tdata=matrix(c(nid,totalt),ncol=1)
tdatas=cbind(tdatas,tdata)
ndata=matrix()
ndata=matrix(c(nid,totaln),ncol=1)
ndatas=cbind(ndatas,ndata)
TP=Sys.time()
DT=as.numeric(TP-T)
cat(" miRNA【分类】完成",as.character(round(j*100/length(SSS),3)),
"%\n 当前速度为",as.character(round(as.numeric(DT),3)),"秒/个\n",
"大约还需要",as.character(round(as.numeric(DT),3)*(length(SSS)-j)),"秒\n")}#over for j
tdatas=tdatas[,-1]
ndatas=ndatas[,-1]
cat("分类完成\n开始统计数据\n")
mirnaid=c()
pvalues=c()
TMFC=c()
NMFC=c()
TSDFC=c()
NSDFC=c()
for(i in 1:ncol(tdatas)){
T=Sys.time()
mirnaid=c(mirnaid,as.character(tdatas[1,i]))
tcvs=c()
ncvs=c()
tcvs=as.numeric(tdatas[2:nrow(tdatas),i])
ncvs=as.numeric(ndatas[2:nrow(ndatas),which(ndatas[1,]==mirnaid[i])])
TMFC=c(TMFC,mean(tcvs))
TSDFC=c(TSDFC,sd(tcvs))
NMFC=c(NMFC,mean(ncvs))
NSDFC=c(NSDFC,sd(ncvs))
pp=t.test(tcvs,ncvs)$p.value
pvalues=c(pvalues,pp)
TP=Sys.time()
DT=as.numeric(TP-T)
cat(" miRNA【统计】完成",as.character(round(i*100/ncol(tdatas),3)),
"%\n 当前速度为",as.character(round(as.numeric(DT),3)),"秒/个\n",
"大约还需要",as.character(round(as.numeric(DT),3)*(ncol(tdatas)-i)),"秒\n")}#over for i
cat("开始制表\n")
ssdata=data.frame(miRNA_ID=mirnaid,
TumorMFC=as.numeric(TMFC),
NormalMFC=as.numeric(NMFC),
LOG2=log(TMFC/NMFC,2),
p_value=pvalues,
TumorSDFC=as.numeric(TSDFC),
NormalSDFC=as.numeric(NSDFC))
NAME=paste("阁下 TCGA-miRNA 分析结果数据",gsub(":","_",Sys.time()),".csv")
if (file.exists("./TCGA")==TRUE){cat("阁下目标文件夹 TCGA 已存在\n")}else{
dir.create("./TCGA", recursive=TRUE)
cat("目标文件夹 TCGA 已为阁下创建\n")}
setwd("./TCGA")
write.csv(ssdata,NAME,row.names=FALSE)
cat(NAME,"已统计完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(adr)}#over if 1O2O
if(substring(control,1,1)=="O"){
volcane(ssdata,signs="miRNA火山图",Larg=Larg)}else{cat("火山图不绘制\n")}
if(substring(control,2,2)=="O"){
colnames(tdatas)=as.character(tdatas[1,])
colnames(ndatas)=as.character(ndatas[1,])
tdatas=tdatas[-1,]
ndatas=ndatas[-1,]
miheatmap=data.frame(rbind(tdatas,ndatas),stringsAsFactors = F)
miheatmap=as.data.frame(lapply(miheatmap,as.numeric))
rownames(miheatmap)=c(paste("id",c(1:nrow(miheatmap)),sep=""))
trow_ann=data.frame(SampleType=c(rep("Primary_Tumor",nrow(tdatas))))
nrow_ann=data.frame(SampleType=c(rep("Solid_Tissue_Normal",nrow(ndatas))))
mirow_annotation=data.frame(rbind(trow_ann,nrow_ann),stringsAsFactors = F)
rownames(mirow_annotation)=c(paste("id",c(1:nrow(mirow_annotation)),sep=""))
for(n in 1:ncol(miheatmap)){
miheatmap[,n]=miheatmap[,n]/mean(as.numeric(ndatas[,n]))}
hotmap(miheatmap,mirow_annotation,signs="miRNA热图",Larg=Larg)}else{cat("热图不绘制\n")}
}#over if mirna
gene_id=as.character(ns$ENSEMBL)
gene_id=V(gene_id)
mirna_id=as.character(ns$IDmirna)
mirna_id=V(mirna_id)
if(substring(control,3,3)!="O"){cat("相关性不绘制\n")}else{
cat("绘制关联图\n")
cat("获取位置\n")
intercaseID=intersect(xsample,ysample)
xsites=c()
ysites=c()
for(z in 1:length(intercaseID)){
xv=which(xsample==intercaseID[z])
yv=which(ysample==intercaseID[z])
xsites=c(xsites,xv)
ysites=c(ysites,yv)}
cat("提取所需数据\n")
cat(gene_id,"\n")
cat(mirna_id,"\n")
xFPKM=c()
yMiread=c()
for(c in 1:length(xsites)){
FPKM=c()
Miread=c()
for(i in 1:length(gene_id)){
PKM=as.numeric(subset(get(paste("D",xsites[c],sep="")),as.character(substring(ENSEMBL,1,15))==gene_id [i])$FPKM)
FPKM=c(FPKM,PKM)}
for(j in 1:length(mirna_id)){
iread=as.numeric(subset(get(paste("S",ysites[c],sep="")),as.character(miRNA_ID)==mirna_id [j])$reads_per_million_miRNA_mapped)
Miread=c(Miread,iread)}
xFPKM=c(xFPKM,FPKM)
yMiread=c(yMiread,Miread)}
xmat=matrix(xFPKM,byrow=TRUE,ncol=length(gene_id))
ymat=matrix(yMiread,byrow=TRUE,ncol=length(mirna_id))
xymat=cbind(xmat,ymat)
xxrownames=c()
for(i in 1:length(gene_id)){
xxrownames=c(xxrownames,as.character(ns$SYMBOL[which(ns$ENSEMBL==as.character(gene_id[i]))]))}
colnames(xymat)=c(xxrownames,mirna_id)
cat("开始绘制chart.Correlation\n")
library(ggplot2) 
library(PerformanceAnalytics)
chart.Correlation(xymat)
if (file.exists("./Graph_TCGA")==TRUE){cat("阁下目标文件夹 Graph_TCGA 已存在\n")}else{
dir.create("./Graph_TCGA", recursive=TRUE)
cat("目标文件夹 Graph_TCGA 已为阁下创建\n")}
setwd("./Graph_TCGA")
gNAME=paste("阁下 Graph_TCGA 已绘制完成",gsub(":","_",Sys.time()),".tiff")
ggsave(filename=gNAME,chart.Correlation(xymat),dpi=Larg,width=12,height=8)
cat("阁下 Graph_TCGA 已绘制完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(adr)
XXXNAME=paste("阁下 TCGA 相关性分析前体结果数据",gsub(":","_",Sys.time()),".csv")
if (file.exists("./TCGA")==TRUE){cat("阁下目标文件夹 TCGA 已存在\n")}else{
dir.create("./TCGA", recursive=TRUE)
cat("目标文件夹 TCGA 已为阁下创建\n")}
setwd("./TCGA")
write.csv(xymat,XXXNAME,row.names=FALSE)
cat(XXXNAME,"已统计完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(adr)
if(nrow(gheatmap)==nrow(miheatmap) &
substring(control,1,2)=="OO"){
Bighot=cbind(gheatmap,miheatmap)
hotmap(Bighot,grow_annotation,signs="Gene+miRNA热图",Larg=Larg)}
chart.Correlation(xymat)
}#over if 3O
if(substring(control,3,3)=="O" & g=="O"&m=="O"){
cat("统计如下：\n基因数目为：",length(gene_id),
"\nmiRNA数目为：",length(mirna_id),
"\n肿瘤组织样本数为：", length(which(xsampletype[xsites]=="Primary Tumor")),
"\n正常组织样本数为：", length(which(xsampletype[xsites]=="Solid Tissue Normal")),"\n")}
if(g=="O"){
cat("Gene统计如下：\n基因数目为：",length(gene_id),
"\nGene肿瘤组织样本数为：", length(which(xsampletype=="Primary Tumor")),
"\nGene正常组织样本数为：", length(which(xsampletype=="Solid Tissue Normal")),"\n")}
if(m=="O"){
cat("miRNA统计如下：\nmiRNA数目为：",length(mirna_id),
"\nmiRNA肿瘤组织样本数为：", length(which(ysampletype=="Primary Tumor")),
"\nmiRNA正常组织样本数为：", length(which(ysampletype=="Solid Tissue Normal")),"\n")}
}