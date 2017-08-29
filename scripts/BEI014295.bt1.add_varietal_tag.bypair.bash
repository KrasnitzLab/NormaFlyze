python /mnt/wigclust5/data/safe/kendall/nla3_88/add.varietal.tag01.bt1.bypair.py <(gunzip -c /mnt/wigclust5/data/safe/kendall/nla3_88/b2/s_1_1_sequence.b1.txt.gz) <(gunzip -c /mnt/wigclust5/data/safe/kendall/nla3_88/b2/s_1_2_sequence.b1.txt.gz) <(/mnt/wigclust5/data/safe/kendall/samtools-0.1.19/samtools view -h /mnt/wigclust5/data/safe/kendall/nla3_88/b2/BEI014295.r1.ucsc.hg19dm6.bt1.sorted.bam) 1 | /mnt/wigclust5/data/safe/kendall/samtools-0.1.19/samtools view -Sbh - > /mnt/wigclust5/data/safe/kendall/nla3_88/b2/BEI014295.r1.ucsc.hg19dm6.bt1.sorted.vt.bypair.bam

/mnt/wigclust5/data/safe/kendall/samtools-0.1.19/samtools index /mnt/wigclust5/data/safe/kendall/nla3_88/b2/BEI014295.r1.ucsc.hg19dm6.bt1.sorted.vt.bypair.bam

python /mnt/wigclust5/data/safe/kendall/nla3_88/add.varietal.tag01.bt1.bypair.py <(gunzip -c /mnt/wigclust5/data/safe/kendall/nla3_88/b2/s_1_1_sequence.b1.txt.gz) <(gunzip -c /mnt/wigclust5/data/safe/kendall/nla3_88/b2/s_1_2_sequence.b1.txt.gz) <(/mnt/wigclust5/data/safe/kendall/samtools-0.1.19/samtools view -h /mnt/wigclust5/data/safe/kendall/nla3_88/b2/BEI014295.r2.ucsc.hg19dm6.bt1.sorted.bam) 2 | /mnt/wigclust5/data/safe/kendall/samtools-0.1.19/samtools view -Sbh - > /mnt/wigclust5/data/safe/kendall/nla3_88/b2/BEI014295.r2.ucsc.hg19dm6.bt1.sorted.vt.bypair.bam

/mnt/wigclust5/data/safe/kendall/samtools-0.1.19/samtools index /mnt/wigclust5/data/safe/kendall/nla3_88/b2/BEI014295.r2.ucsc.hg19dm6.bt1.sorted.vt.bypair.bam