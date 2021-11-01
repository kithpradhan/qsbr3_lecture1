library(GEOquery) 


# load data

dat <- getGEO("GSE7148")		# does this work for you?

# if not 
# load("dat_gse7148.RData")

es <- dat[[1]]
edat <- exprs(es)
pdat <- pData(es)

# sanity checking the data...

range(edat, na.rm=T)

boxplot(edat)
summary(edat)

# missing data?

f <- apply(edat, 1, function(x){ sum(is.na(x)) })

lonely.index <- pdat[,10]
table(lonely.index)

# run T-test

t.test(edat[1,lonely.index=="HighLonely"], edat[1,lonely.index=="LowLonely"])

runTP <- function(x,y){ t.test(x[y=="HighLonely"], x[y=="LowLonely"])$p.value}

tpvals <- apply(edat, 1, runTP, y=lonely.index)
hist(tpvals)
summary(tpvals)

# adjust for multiple testing

tapvals <- p.adjust(tpvals, "BH")
summary(tapvals)

# make volcano plot

avg.lonely <- apply(edat[,lonely.index=="HighLonely"], 1, mean)
avg.social <- apply(edat[,lonely.index=="LowLonely"], 1, mean)
fch <- avg.lonely - avg.social
logp <- -log(tpvals, 10)

plot(fch, logp, pch=1, xlab="Log2(L/S)", ylab="-Log10(Pval)", main="Volcano Plot")

abline(v=c(-1,1), col="red", lwd=3, lty=2)


# run limma analysis

library(limma)

lonely.index <- as.factor(lonely.index)

design <- model.matrix(~0+lonely.index)
colnames(design) <- levels(lonely.index)
fit <- lmFit(es, design)

cont.diff <- makeContrasts(HighLonely-LowLonely, levels=levels(lonely.index))
fit2 <- contrasts.fit(fit, cont.diff)
fit2 <- eBayes(fit2) 

tt <- topTable(fit2, n=nrow(exprs(es)), adjust="fdr")
sum(tt$adj.P.Val < .05)

#test out github
2+2

3+3

5+5
