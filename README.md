# BAM STATS

Statistics for bam file. Including mapping summary(mapping ratio, identity) and
mapped reads fragment information(position on reference, position on query, strand, mapping quality, mapping identity)

## Dependency

* [htslib](https://github.com/samtools/htslib)

## Install

Build with cmake: 

```shell
git clone --recursive git@github.com:Archieyoung/bam_stats.git
cd bam_stats
mkdir build && cd build
cmake ..
make
```
an executable binary file bam_stats will appear in build/bin.


