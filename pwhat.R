install.packages("vsm_0.1.tar.gz", repos = NULL, type="source")
load("CCFA.RData")

require(vsm)

require(mgcv)
require(fda)
require(vows)


set.seed(789)
methods <- c("OLS", "GLS", "Post-smoothed")

# Fit the models
obj = vector("list", 3)
names(obj) = methods

cat("\nDoing OLS \n")
print(system.time(obj[[1]] <- tps.fd(cc, age, argvals=arclength, method="P-OLS")))
  
cat("\nDoing GLS \n")
print(system.time(obj[[2]] <- tps.fd(cc, age, argvals=arclength, method="P-GLS")))

cat("\nDoing Post-smoothed \n")
lsp.f. = -20:10
print(system.time({cv.smooth <- cv.vsm2s(cc, age, argvals=arclength, lsp.f = lsp.f.)
                   idx.cv = which.min(apply(cv.smooth$cvmat, 2,sum))
                   obj[[3]] <- vsm2s(cc, age, argvals=arclength, lsp.f = lsp.f.[idx.cv], varest=TRUE)}))
cat("lsp.f chosen is", lsp.f.[idx.cv],"\n")


############################################################################
# MAKE FIGURES
############################################################################

# Rainbow plot of raw data
quartz("Figure 1 right", 10)
rainbow.plot(cc, age, ylim=range(cc), argvals=arclength, subset=unique(c(which.min(age), which.max(age), sample(length(age), 30))), xlab="Arc length (mm)", ylab="Fractional anisotropy", posy=c(.09,.44), zval=10*1:4, digit=0, cex.lab=1.5, cex.axis=1.5, legend=FALSE)

##############################################################
##############################################################

# Pointwise RLRT
rlrt.tract = rlrt.mp(cc, age, loginvsp=-30:5) 

quartz("Figure 3", 12)
cexlab=1.4
cexmain=1.8
layout(matrix(c(1,2,1,3,1,4),2))

# Plot RLRT p-values
plot(arclength, rlrt.tract$pvalue, main=" ", xlab="Arc length (mm)", ylab="p-value", log="y", cex.lab=cexlab, cex.main=cexmain)  
abline(h=.05, col=2)

# Scatterplots for three selected voxels
whichvox = c(44,60,97)
for (vv in whichvox) abline(v=arclength[vv], lty=2, col='blue')
region.names = c("Sensorimotor", "Posterior prefrontal", "Anterior prefrontal")
for (j in 1:3) {
    plot(age, cc[ , whichvox[j]], main=region.names[j], xlab="Age (years)", ylab=if (j==1) "Fractional anisotropy" else "", cex.main=cexmain, cex.lab=cexlab)
    mod = gam(cc[ , whichvox[j]] ~ s(age), method="REML")
    agegrid=7:48
    pred = predict(mod, newdata=data.frame(age=agegrid), se.fit=TRUE)
    matlines(agegrid, cbind(pred$fit, pred$fit-2*pred$se, pred$fit+2*pred$se), col="blue", lty=c(1,2,2))
}


##############################################################
##############################################################

# Figure 4

# Top: Rainbow plots for fitted values
quartz("Figure 4 top", 18)
par(mfrow=c(1,3), mar=c(4.5, 4.5, 3, 2)) 
ageseq = 7:48
mains=c("Ordinary least squares", "Generalized least squares", "Post-smoothed")
yrange = NULL
for (i in 1:3) yrange = range(yrange, c(predict(obj[[i]], newx=ageseq)))
for (i in 1:3) {
    plot(obj[[i]], newx=ageseq, xlab="Arc length (mm)", ylab=if (i==1) "Fractional anisotropy" else "", main=mains[i], cex.main=3, cex.lab=1.8, ylim=yrange)
    for (vv in c(44,60,97)) abline(v=arclength[vv], lty=2, col='blue')
}

##############################################################

# Bottom: Pointwise DF and R^2 vs. CV 

# Obtain functional R^2 for each method
# Trapezoidal approximation
trap = function(vec, argvals=1:length(vec))
              sum(diff(argvals) * (head(vec,-1)+tail(vec,-1))) / 2

R2vsm = function(Y,mod) {
    fitted = if (length(dim(mod$fitted))==3) mod$fitted[,,1] else mod$fitted
    Ybarmat = outer(rep(1,nrow(Y)), colMeans(Y))
    1 - sum(apply((Y-fitted)^2, 1, trap, argvals=arclength)) / sum(apply((Y-Ybarmat)^2, 1, trap, argvals=arclength))	
}

R2vec = c()
for (mm in 1:3) R2vec[mm] = R2vsm(cc, obj[[mm]])  

####################

# Obtain CV score for each method, based on repeated multifold CV (takes a while!)
n.subject <- length(age)
n.folds <- 5
n.splits <- 10  

obj2 = vector("list", 3)
names(obj2) = methods

cv.array <- array(NA, c(3, n.folds, n.splits))

set.seed(1567)  
for (nn in 1:n.splits) {
    cat("\n******* Split", nn, "*******\n")
    
    all.folds <- split(sample(1:n.subject), rep(1:n.folds, length=n.subject))
    cv.tab <- matrix(NA, 3, n.folds)
    rownames(cv.tab) = methods
    
    for (kk in 1:n.folds) {
        cat("\n******* Fold", kk, "*******\n")
        IDX.cv <- all.folds[[kk]]
        IDX.train <- (1:n.subject)[-IDX.cv]

        cat("\nDoing OLS \n")
        obj2[[1]] = tps.fd(cc[IDX.train,], age[IDX.train], argvals=arclength, method="P-OLS")          
        # y.pred <- predict.tps.fd(obj2[[1]], age[IDX.cv])
        y.pred <- predict(obj2[[1]], age[IDX.cv])
        cv.tab[1,kk] <- sum((y.pred - cc[IDX.cv,])^2)

        cat("\nDoing GLS \n")
        obj2[[2]] = tps.fd(cc[IDX.train,], age[IDX.train], argvals=arclength, method="P-GLS")          
        # y.pred <- predict.tps.fd(obj2[[2]], age[IDX.cv])
        y.pred <- predict(obj2[[2]], age[IDX.cv])
        cv.tab[2,kk] <- sum((y.pred - cc[IDX.cv,])^2)

        cat("\nDoing Post-smoothed \n")
        cv.smooth = cv.vsm2s(cc[IDX.train,], age[IDX.train], argvals=arclength, lsp.f = lsp.f.)
        idx.cv = which.min(apply(cv.smooth$cvmat, 2,sum))
        obj2[[3]] = vsm2s(cc[IDX.train,], age[IDX.train], argvals=arclength, lsp.f = lsp.f.[idx.cv], varest=FALSE)
        cat("lsp.f chosen is", lsp.f.[idx.cv],"\n")          
        # y.pred <- predict.vsm2s(obj2[[3]], age[IDX.cv])
        y.pred <- predict(obj2[[3]], age[IDX.cv])
        cv.tab[3,kk] <- sum((y.pred - cc[IDX.cv,])^2)   
    }
    cv.array[,,nn] <- cv.tab
}
cv.mean = apply(cv.array, 1, sum) / (n.splits*length(age))

quartz("Figure 4 bottom", 16)
layout(mat=matrix(c(1,2),nrow=1), widths=c(1.8,1.2))
colvec = c(2,4,1)
ltyvec = c(3,5,1)
all.pwdf = c()
for (k in 1:3) all.pwdf = c(all.pwdf, obj[[k]]$pwdf)
with(obj[["Post-smoothed"]], plot(argvals, pwdf.raw, ylim=range(all.pwdf), xlab="Arc length (mm)", ylab="Pointwise degrees of freedom", log="y", cex.lab=1.3, main="", lty=1, cex.main=1.6))
for (k in 1:3) with(obj[[k]], lines(argvals, pwdf, col=colvec[k], lty=ltyvec[k]))
for (vv in c(44,60,97)) abline(v=arclength[vv], lty=2, col='grey')

plot(R2vec, cv.mean, xlab=expression(R^2), ylab="Cross-validation score", log="x", col=colvec, cex.lab=1.3, main="", cex.main=1.6)
text(R2vec+c(-.025,.0006,0), cv.mean+c(-.001,.001,.001), methods, col=colvec, cex=1)


##############################################################
##############################################################

# Figure 5: FA vs. age in eight prefrontal cortex voxels
vox = 91:98
colo = rainbow(length(vox), start=.55, end=.75)
agegrid = seq(7,48,,600)
fdn = list("Age", "voxel", "FA")
xy.peak = function(fdobj) {
    fd.eval = eval.fd(agegrid, fdobj)
    xvec = agegrid[apply(fd.eval, 2, which.max)]
    yvec = apply(fd.eval, 2, max)
    cbind(xvec, yvec)
}

fdlist = vector("list",4)
fdlist[[1]] = with(obj[[1]], fd(basis=basis.x, coef=tcrossprod(coef, B.f[vox,]), fdnames=fdn))
fdlist[[2]] = with(obj[[2]], fd(basis=basis.x, coef=tcrossprod(coef, B.f[vox,]), fdnames=fdn))
fdlist[[3]] = with(obj[[3]], fd(basis=basis.x, coef=xcoef[,vox,1], fdnames=fdn))
fdlist[[4]] = with(obj[[3]], fd(basis=basis.x, coef=coef.raw[,vox], fdnames=fdn))

yrange = range(eval.fd(agegrid, fdlist[[1]]))
for (i in 2:4) yrange = range(c(yrange, range(eval.fd(agegrid, fdlist[[i]]))))

peakplot = function(fdobj, ...) {
    plot(fdobj, col=colo, lty=1, cex.main=1.8, cex.lab=1.4, ylim=yrange, ...)
    points(xy.peak(fdobj), col=colo, pch=16)
}

quartz("Figure 5", 13)
layout(mat=matrix(c(0,1,1,0,2,3,4,5), nrow=2, byrow = TRUE))

# Brain map with voxels of interest highlighted
par(mar = c(1,5,1,5))  
xsub=104:166
ysub=72:98
othervox = matrix(NA, nrow(anat), ncol(anat))
for (k in (1:107)[-vox]) othervox[vc[k,1],vc[k,2]] = 1
image((-125:92)[xsub], (-71:110)[ysub], anat[xsub,ysub], col=grey(0:40/40), axes=FALSE, xlab="", ylab="")  
image((-125:92)[xsub], (-71:110)[ysub], othervox[xsub,ysub], add=TRUE, col="red")
for (k in 1:length(vox)) {
    thisvox = matrix(NA, nrow(anat), ncol(anat))
    thisvox[vc[vox[k],1],vc[vox[k],2]] = 1
    image((-125:92)[xsub], (-71:110)[ysub], thisvox[xsub,ysub], add=TRUE, col=colo[k])
}

# Fitted curves for those voxels
par(mar = c(4.5, 4.5, 3, 2))
peakplot(fdlist[[4]], main="Separate smooths", xlab="Age (years)", ylab="Fractional anisotropy")
for (i in 1:3) peakplot(fdlist[[i]], main=mains[i],xlab="Age (years)", ylab="")
