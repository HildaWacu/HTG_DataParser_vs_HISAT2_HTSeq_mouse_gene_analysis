# Preparation x window
col.vec <- c("black","red")
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2), cex.main=0.7, mgp=c(2,1,0), mai=c(0.5,0.5,0.5,0.5))

# Genes as points and samples as dimensions
a <- c(1, 4)
b <- c(2, 3)
m1 <- rbind(a, b)
print(m)

colnames(m) <- c("s1", "s2")
plot(m , pch=16, xlim=c(0,3), 
     ylim=c(0,5), main="Genes as points and samples as dimensions",
     col=col.vec,
     panel.first=grid(lty=1))
suppressWarnings(arrows(a[1], a[2], b[1], b[2], angle=0))
suppressWarnings(arrows(a[1], a[2], b[1], a[2], angle=0, lty=3))
suppressWarnings(arrows(b[1], a[2], b[1], b[2], angle=0, lty=3))
text(m, lab=rownames(m), pos = 2, offset = 1)
text(7, 4,label="dist(a,b)")


# Sample as points and genes as dimensions
plot(t(m) , pch=16, 
     ylim=c(0,5), xlim=c(0,5), 
     main="Sample as points and genes as dimensions",
     panel.first=grid(lty=1))
text(t(m), lab=colnames(m), pos = 2, offset = 1)


# x axis correspond to samples and the y axis represent the intensities

matplot(t(m) , 
        ylim=c(0,5), xlim=c(0,3), 
        main="x axis for samples (n=2) and y axis for intensities", 
        ylab="Intensities",
        xlab="samples",
        xaxt = "n",
        type="n")
grid(lty=1)

axis(1, 0:4, c("", "s1", "s2", "", ""))
suppressWarnings(arrows(1, a[1], 2, a[2], angle=0))
suppressWarnings(arrows(1, b[1], 2, b[2], angle=0))

matpoints(t(m) , pch=16, 
          col=col.vec)


a <- c(7, 7, 7, 7, 6, 6, 10, 10)
b <- c(8, 2, 5,  6, 1, 6, 1, 4)
m <- rbind(a, b)
matplot(t(m),
        xlim=c(0,10), ylim=c(0,12),
        main="x axis for samples (n=8) and y axis for intensities",
        pch=16,
        col="black",
        ylab="Intensities",
        xlab="samples",
        lty=1,
        type="n")
grid(lty=1)

for(i in 1:length(a)){
  suppressWarnings(arrows(i, a[i] , i, b[i], angle=0, col=col.vec[i], lty=3))
}

matpoints(t(m),
          pch=16,
          type="b",
          lty=1)

points(a, type="p", col=col.vec[1],  pch=16)
points(b, type="p", col=col.vec[2],  pch=16)

