CXXFLAGS = -g -Wall --std=c++11

hts_base=./htslib
bam_stats: from_cigar.o bam_stats.o
	g++ $(CXXFLAGS) -I $(hts_base) -L $(hts_base) -lhts -o bam_stats from_cigar.o bam_stats.o

from_cigar.o: from_cigar.cpp from_cigar.h
	g++ $(CXXFLAGS) -c -I $(hts_base) -o $@ from_cigar.cpp

bam_stats.o: bam_stats.cpp bam_stats.h
	g++ $(CXXFLAGS) -c -I $(hts_base) -o $@ bam_stats.cpp

.PHONY: clean

clean:
	-rm -f *.o
