---
layout: docs
title: Processing of DNase-seq data from ENCODE
---

<div id="sidebar">
<a href="#">Processing of DNase-seq data from ENCODE</a>
<a href="part2.html">`seqOutBias` to generate scaled bigWig files</a>
</div>

# `seqOutBias` to generate scaled bigWig files

The software `seqOutBias` will scale aligned bam read counts by the ratio of genome-wide observed
read counts to the sequence based counts for each k-mer. The k-mer counts take into account the
mappability at a given read length. The `seqOutBias` program allows for flexibility in specifying k-mer
size, strand-specific offsets, and spaced k-mers.

## Using `seqOutBias` to scale DNase-seq files by 6-mer nick preference

The specificity of DNase is strongly influenced by the three bases that flank each side of the DNase cut
site (Figure 1) (He *et al.*, 2014; Yardımcı *et al.*, 2014).

<img src="assets/images/DNase_cut_preference.jpg" style="width:30%;cursor:zoom-in" onclick="document.getElementById('modal01').style.display='block'">

  <div id="modal01" class="w3-modal" onclick="this.style.display='none'">
    <span class="w3-button w3-hover-red w3-xlarge w3-display-topright">&times;</span>
    <div class="w3-modal-content w3-animate-zoom">
      <img src="img_fjords.jpg" style="width:100%">
    </div>
  </div>

ACGGGATATGATGACCAGATGACA
TGCCCTATACTACTGGTCTACTGT
DNase nicking
Figure 1: The six base pair window centered on the DNase nick dictates cleavage preference. (He *et al.*,
2014)

`seqOutBias` will calculate the genome-wide occurences of each specified k-mer centered on the DNase
nick site accounting for the mappability of the specified read length (note the default is –read-size=36).
For each case below the offsets are half the value of the kmer-size parameter, which is the sequence
length (k-mer) that surronds the nick-site and influences specificity, therefore the program will calculate
the frequency of k-mers centered on the nick-site. Experimentally, we assume that we are equally
likely to sequence either end of a DNase nick site, so the –shift-counts parameter is used to shift
the Crick strand alignments in line with the Watson strand alignments (Figure 2). DNase nicks can be
offset or in line, as shown. Note that generating the mappability files for a given genome and read
length is time-consuming, but once these files are made, `seqOutBias` will recognize the existence of
these files and avoid timely recomputing and regeneration of these files.

```bash
bam=UW_MCF7_both.bam
`seqOutBias` hg38.fa $bam --no-scale --bw=MCF7_0-mer.bigWig --shift-counts --skip-bed
`seqOutBias` hg38.fa $bam --kmer-size=6 --bw=MCF7_6-mer.bigWig --plus-offset=3 --minus-offset=3 --shift-counts --skip-bed
`seqOutBias` hg38.fa $bam --kmer-size=10 --bw=MCF7_10-mer.bigWig --plus-offset=5 --minus-offset=5 --shift-counts --skip-bed
bam=IMR90_Naked_DNase.bam
`seqOutBias` hg38.fa $bam --no-scale --bw=Naked_0-mer.bigWig --shift-counts --skip-bed
`seqOutBias` hg38.fa $bam --kmer-size=6 --bw=Naked_6-mer.bigWig --plus-offset=3 --minus-offset=3 --shift-counts --skip-bed
`seqOutBias` hg38.fa $bam --kmer-size=10 --bw=Naked_10-mer.bigWig --plus-offset=5 --minus-offset=5 --shift-counts --skip-bed
```
