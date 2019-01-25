library(readr)
fileName = "simul_15_0.9_2"
load(paste(fileName,".RData",sep=""))
greyscale=gray((0:32)/32)

################## HESS
## gamma
gamma_hess = as.matrix( read_table(file=paste(sep = "","results/",fileName,"_HESS_gamma_out.txt"),col_names=FALSE))
# image(gamma_hess,col=greyscale,main="hess")

## G
# image(diag(s)[,s:1])

################## SSUR
## gamma
gamma_ssur = as.matrix( read_table(file=paste(sep = "","results/",fileName,"_SSUR_gamma_out.txt"),col_names=FALSE))
# image(gamma_ssur,col=greyscale)

## G
g_ssur = as.matrix( read_table(file=paste(sep = "","results/",fileName,"_SSUR_G_out.txt"),col_names=FALSE))
# image(g_ssur[,s:1],col=greyscale)

################## dSUR
## gamma
gamma_dsur = as.matrix( read_table(file=paste(sep = "","results/",fileName,"_dSUR_gamma_out.txt"),col_names=FALSE))
# image(gamma_ssur,col=greyscale)

######
par(mfrow=c(2,3))
image(gamma_ssur,col=greyscale,main="SSUR")
image(gamma_hess,col=greyscale,main="HESS")
image(gamma_dsur,col=greyscale,main="dSUR")

image(gamma[-1,],col=greyscale,mai="true")
image(as.matrix(G)[,s:1],col=greyscale,main="true G")
image ((g_ssur+diag(s))[,s:1],col=greyscale,main="SSUR G")
par(mfrow=c(1,1))


