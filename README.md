# BAM STATS

Statistics for bam file. Including mapping quality, aligned lenght, splited reads number, splited read fragments number, splited reads fragment lenght and so on.

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


