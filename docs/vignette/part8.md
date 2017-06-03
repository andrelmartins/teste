---
layout: docs
title: Characterizing the enzymatic clean up and ligation sequence bias
---

# Characterizing the enzymatic clean up and ligation sequence bias

We had previously assumed that the three bases upstream and downstream of a DNase-nick site are equally likely to have adapters ligated and to be sequenced (Figure 2). However, we find that this is not the case and there is a bias in which 3-mer is ultimately detected by sequencing. Note that there is no inbalance for reverse palindromic 6-mers, for example: GCATGC

## Plotting post-nicking enzymatic sequence biases

For each DNA nick we tally the number of times a plus strand and minus strand read detects the nick event. For simplicity, we assume that the mappability of plus and minus strand-aligned reads are the same.
The analysis below is exclusively for the plus strand analysis, but a minus or combined strand analysis gives the same results.

```r
source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
counts.table = read.table('~/DNase_ENCODE/hg38_36.6.3.3.IMR90_Naked_DNase.txt')

ligation = cbind(counts.table, substring(counts.table[,2],1,3)) 
ligation[,8] = apply(ligation, 1, function(row) revcomp(row[7])) 
ligation[,9] = substring(counts.table[,2],4,6)
colnames(ligation) = c(colnames(counts.table), 'V7', 'rcKmerUp','KmerDown')

mat = data.frame(matrix(nrow=64, ncol= 64)) 
count = 0
for (mer in unique(ligation$KmerDown)) {
    count = count + 1
    temp = ligation[ligation$KmerDown == mer,] rto = temp[,6]/temp[,5]
    mat[,count] = rto
}
colnames(mat) = unique(ligation$KmerDown)
mat = do.call(data.frame,lapply(mat, function(x) replace(x, is.infinite(x),NA))) 
rownames(mat) = unique(temp[,8])
mat = mat[order(rownames(mat)) , order(colnames(mat))]
mat = as.matrix(mat)

pdf('lig_bias_matrix.pdf', width=10.4, height=9.5)
heatmap.2(log(mat, base=10), col=colorpanel(30, "blue", "white","red"),
    symbreaks=T,scale="none",na.rm=TRUE,dendrogram = 'none', symm =TRUE, 
    density.info=c("none"), key.xlab = expression('log'[10]*' ratio of bias (x-axis/y-axis)'), 
    key.title= '',trace=c("none"),Rowv = FALSE, lhei=c(0.75,4), lwid = c(1.2, 4))
dev.off()

pdf('lig_bias_matrix_row.pdf', width=10.4, height=9.5) 
heatmap.2(log(mat, base=10), col=colorpanel(30, "blue", "white","red"),
    symbreaks=T,scale="none",na.rm=TRUE,dendrogram = 'row', symm =TRUE,
    density.info=c("none"), key.xlab = expression('log'[10]*' ratio of bias (x-axis/y-axis)'), 
    key.title= '',trace=c("none"), Rowv = TRUE, lhei=c(0.75,4), lwid = c(1.2, 4))
dev.off()

n.diag = mat
diag(n.diag) = NA
df = data.frame(x=c(n.diag, diag(mat)), group=factor(c(rep("non-Palindromic", length(n.diag)),
    rep("Palindromic", length(diag(mat))))))

pdf('lig_bias_bwplot.pdf', width=3, height=4) 
trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=2),
                box.rectangle = list(col = 'black', lwd=1.6),plot.symbol = list(col='black', lwd=1.6, pch ='.')) 
bwplot(log(x, base = 2) ~ group, df,
       scales=list(x=list(relation = "free", rot = 45)),
       ylab = expression('log'[2]*' ratio of detection bias'),
       pch = '|', 
       col= 'black' )
dev.off()

pdf('avg_ligation_bias.pdf', width=12, height=4) 
print(barchart((colMeans(mat, na.rm =TRUE)~colnames(mat)),
                col='grey85',
                ylim = c(0, max(colMeans(mat, na.rm =TRUE))+ 0.01 * max(colMeans(mat, na.rm =TRUE))), 
                ylab=paste('relative ligation preference of each 3-mer', sep = ' '),
                xlab = '3-mer',
                origin = 0,
                scales=list(x=list(rot=45)),
                panel=function(...) {
                    panel.barchart(...)
                    panel.abline(h=1, lty= 2, col = 'grey40') 
                }
            ))
dev.off()

#for (mer in unique(ligation$KmerDown)) {
for (mer in c('AAA', 'ATA', 'ATC', 'GAT','GAC')) { 
    plot.barchart.lig(ligation[ligation$KmerDown == mer,],
                    filename=paste('ligation_bias_',mer,'.pdf', sep = ''), w = 12, h = 4) 
}

ligation$ratio = ligation[,6]/ligation[,5]
x.ligation = ligation[with(ligation, order(ratio)), ]

pswm.func(head(x.ligation[,2], 205), out = 'low_205.txt', positions = 6) pswm.func(tail(x.ligation[,2], 205), out = 'high_205.txt', positions = 6)
```

<img src="{{site.url}}/{{site.baseurl}}/assets/images/ligation_bias_GAC.jpg" style="width:30%;cursor:zoom-in" onclick="document.getElementById('modal15').style.display='block'">
<div id="modal15" class="w3-modal" onclick="this.style.display='none'">
    <span class="w3-button w3-hover-red w3-xlarge w3-display-topright">&times;</span>
    <div class="w3-modal-content w3-animate-zoom">
        <img src="{{site.url}}/{{site.baseurl}}/assets/images/ligation_bias_GAC.jpg" style="width:100%">
        <div class="w3-modal-caption">
            Figure 15: For all sequence-detected DNase-nicked 6-mers that end in ’GAC’ we compare the ratio of sequence reads that start with ’GAC’ (’<strong>GAC</strong>CAGATGACA’ in Figure 2) to the oppositely oriented 3-mer (’<strong>ATC</strong>ATATCCCGT’ in Figure 2).
        </div>
    </div>
</div>

figure 16 ...