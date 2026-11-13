run_scode <- function(fdata, ftime, dir, tfnum, pnum, cnum, maxite) {
maxB <- 2.0
minB <- -10.0

dir.create(dir, recursive = TRUE, showWarnings = FALSE)

X <- as.matrix(read.table(fdata, sep="\t"))[1:tfnum,1:cnum]
W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
WZ <- matrix(nrow=tfnum, ncol=cnum)

pseudotime <- read.table(ftime, sep="\t")[1:cnum,2]
pseudotime <- pseudotime/max(pseudotime)

new_B <- rep(0, pnum)
old_B <- rep(0, pnum)

RSS <- Inf
for(i in 1:pnum){
    new_B[i] <- runif(1, min=minB, max=maxB)
    old_B[i] <- new_B[i]
}

sample_Z <- function(){
    for(i in 1:pnum){
        for(j in 1:cnum){
            Z[i,j] <<- exp(new_B[i]*pseudotime[j]) + runif(1, min=-0.001, max=0.001)
        }
    }
}

for(ite in 1:maxite){
    target <- floor(runif(1, min=1, max=pnum+1))
    new_B[target] <- runif(1, min=minB, max=maxB)

    if(ite == maxite){
        for(i in 1:pnum){
            new_B[i] <- old_B[i]
        }
    }

    sample_Z()

    for(i in 1:tfnum){
        X.lm <- lm(X[i,] ~ t(Z)-1)
        for(j in 1:pnum){
            W[i,j] <- X.lm$coefficients[j]
        }
        WZ[i,] <- W[i,] %*% Z
    }

    tmp_RSS <- sum((X-WZ)**2)
    if(tmp_RSS < RSS){
        RSS <- tmp_RSS
        old_B[target] <- new_B[target]
    }
    else{
        new_B[target] <- old_B[target]
    }
}

write.table(RSS, paste(dir,"/RSS.txt",sep=""), row.names=F, col.names=F, sep="\t")

write.table(W, paste(dir,"/W.txt",sep=""), row.names=F, col.names=F, sep="\t")

B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
for(i in 1:pnum){
    B[i,i] <- new_B[i]
}
invW <- MASS::ginv(W)
A <- W %*% B %*% invW

write.table(A, paste(dir,"/A.txt",sep=""), row.names=F, col.names=F, sep="\t")
write.table(B, paste(dir,"/B.txt",sep=""), row.names=F, col.names=F, sep="\t")

}
