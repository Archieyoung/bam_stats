#include <iostream>
#include <algorithm> 
#include <numeric>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <getopt.h>


#include "bam_stats.h"

inline bool is_unmapped(const uint16_t &_flag) {
    return _flag & BAM_FUNMAP;
}

inline int error_prob_to_phred(double ep) {
    if (ep == 0) {
        return 41; // max phred score
    }
    double phred_d = -10*std::log10(ep); 
    if ((phred_d) > 41) {
        return 41;
    } else {
        return round(phred_d);
    }
}

inline double phred_to_error_prob(int phred) {
    return std::pow(10.0, -phred/10.0);
}

int mean_phred(std::vector<int> &q_qual) {
    double sum_ep = 0;
    int32_t n = 0;
    auto beg = q_qual.cbegin();
    auto end = q_qual.cend();
    while (beg != end) {
        sum_ep += phred_to_error_prob(*beg);
        ++beg;
        n += 1;
    }
    if (n != 0) {
        return error_prob_to_phred(sum_ep/n);
    } else {
        return -1; // error, quality string is empty
    }
}

inline float mean_identity(std::vector<float> &identity_vec) {
    double sum = std::accumulate(identity_vec.cbegin(), identity_vec.cend(), 0.0);
    return sum/identity_vec.size();
}

// uint64_t total_reads(std::vector<std::string> &qnames) {
//     std::sort(qnames.begin(), qnames.end()); // sort inplace
//     auto end_uniq = std::unique(qnames.begin(), qnames.end());
//     qnames.erase(end_uniq, qnames.end());
//     return qnames.size();
// }

int get_fragment_qual(bam1_t *b, const uint32_t &l_query,
    const uint32_t &query_start, const uint32_t &query_end) {
    /*get mapping fragment quality string, and then calculate mean sequencing
    quality of it*/
    // get quality string
    // bam_get_qual function get the phred score, no nead to minus 33
    uint8_t *bam_quality_string = bam_get_qual(b);

    std::vector<int> fragment_qual;
    // fragment mean quality
    if (bam_is_rev(b)) {
        for (uint32_t i = l_query-query_start-1; i <= l_query-query_end-1; ++i) {
            fragment_qual.push_back(int(*(bam_quality_string+i)));
        }
    } else {
        for (uint32_t i = query_start; i <= query_end; ++i) {
                fragment_qual.push_back(int(*(bam_quality_string+i)));
        }
    }
    // std::cout << fragment_qual << std::endl;
    int fragment_mean_qual = mean_phred(fragment_qual);
    return fragment_mean_qual;
    
}

int bam_stats(const char *input_bam, const std::string prefix, bool get_qual)
{
    // open bam
    samFile *fp = sam_open(input_bam, "r");
    bam_hdr_t *h = sam_hdr_read(fp);
    bam1_t *b;

    b = bam_init1();

    std::ofstream out_hd;
    out_hd.open(prefix+".mapping.fragment.stats.txt");

    if (get_qual) {
        // output table header
        out_hd << "#CHROM\tREF_START\tREF_END\tQUERY_NAME\tQUERY_POS1\tQUERY_POS2\t"
"MAPPING_QAULITY\tFRAGMENT_MEAN_Q\tFRAGMENT_IDENTITY" << std::endl;
    } else {
        // output table header
        out_hd << "#CHROM\tREF_START\tREF_END\tQUERY_NAME\tQUERY_POS1\tQUERY_POS2\t"
"MAPPING_QAULITY\tFRAGMENT_IDENTITY" << std::endl;
    }
    
    // key qnames, value query length
    std::map<std::string, uint32_t> qlen;
    // unmapped reads number
    uint64_t unmapped_num = 0;

    // total bases
    uint64_t total_bases = 0;
    // mapped bases
    uint64_t mapped_bases = 0;

    // identity vector
    std::vector<float> identity_vec;
    int ret;
    // read record one by one
    while ((ret = sam_read1(fp, h, b) >= 0)) {
        // query name
        const char *query_name = bam_get_qname(b);
        const uint32_t l_query = b->core.l_qseq;

        if (qlen.find(query_name) == qlen.end()) {
            qlen[query_name] = l_query;
            total_bases += l_query;
        }

        const uint16_t bam_flag = b->core.flag;
        // skip unmapped record
        if (is_unmapped(bam_flag)) {
            // std::cout << "Unmapped." << std::endl;
            ++unmapped_num;
            continue;
        }

        const uint32_t _cigar_array_len = b->core.n_cigar;
        
        const uint32_t pos = b->core.pos;
        const char *ref_name = h->target_name[b->core.tid];
        const uint16_t mapping_qual = b->core.qual;

        // std::cout << _cigar_array_len << std::endl;
        // from_cigar constructor
        from_cigar _cp(bam_get_cigar(b), _cigar_array_len);        
        const uint32_t query_start = _cp.get_query_start(bam_flag, l_query);
        const uint32_t query_end = _cp.get_query_end(bam_flag, l_query);
        // mapped bases is the span of the segment, cliping not included
        uint32_t qspan = abs(query_end - query_start) + 1;
        mapped_bases += qspan;

        // "NM" tag
        uint8_t *NM_tag;
        // check tag existence
        if (bam_aux_get(b, "NM") != nullptr) {
            NM_tag = bam_aux_get(b, "NM");     
        } else {
            std::cerr << "Can't find \"NM\" tag" << std::endl;
            continue;
        }

        int64_t NM_tag_value;
        // check value get
        NM_tag_value = bam_aux2i(NM_tag);

        // mapping identity
        float percent_identity = 100*(1 - float(NM_tag_value) / float(qspan));
        identity_vec.push_back(percent_identity);

        int fragment_mean_qual = 0;

        if (get_qual) {
            fragment_mean_qual = get_fragment_qual(b, l_query, query_start, query_end);
            out_hd << ref_name
            << "\t"
            << _cp.get_ref_start(pos)
            << "\t"
            << _cp.get_ref_end(pos)
            << "\t"
            << query_name
            << "\t"
            << query_start
            << "\t" 
            << query_end
            << "\t"
            << mapping_qual
            << "\t"
            << fragment_mean_qual
            << "\t"
            << percent_identity
            << std::endl;
        } else {
            out_hd << ref_name
            << "\t"
            << _cp.get_ref_start(pos)
            << "\t"
            << _cp.get_ref_end(pos)
            << "\t"
            << query_name
            << "\t"
            << query_start
            << "\t" 
            << query_end
            << "\t"
            << mapping_qual
            << "\t"
            << percent_identity
            << std::endl;
        }
    }

    out_hd.close();

    std::ofstream out_hd1;

    out_hd1.open(prefix+".mapping.summary.txt");

    out_hd1 << "Total_Reads\tMapped_Reads\t"
"Mapped_Reads_Rate\tTotal_Bases\tMapped_Bases\tMapped_Bases_Rate\t"
"Mean_Percent_of_Identity" << std::endl;

    uint64_t total_reads_num = qlen.size();
    uint64_t mapped_reads_num = total_reads_num - unmapped_num;
    out_hd1 << total_reads_num << "\t"
            << mapped_reads_num << "\t"
            << float(mapped_reads_num)/total_reads_num << "\t"
            << total_bases << "\t"
            << mapped_bases << "\t"
            << float(mapped_bases)/total_bases << "\t"
            << mean_identity(identity_vec) << "\t"
            << std::endl;
   
    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fp);
    return 0;
} // need simplify


void usage() {
    std::cerr << "Long Reads mapping Fragments stats" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "-h, --help                print this message and exit\n"
        << "-b, --bam, FILE           input bam file [default: None]\n"
        << "-p, --prefix, FILE        output file prefix [default: None]\n"
        << "-q, --get_qual            get fragment sequencing quality [default FALSE]"
        << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc <= 1) {
        usage();
        return 0;
    }

    // long_option array
    static const struct option long_options[] = {
        {"bam", required_argument, 0, 'b'},
        {"out", required_argument, 0, 'o'},
        {"get_qual", no_argument, 0, 'q'},
        {"help", no_argument, 0, 'h'}
    };

    int c = 100, long_idx;
    const char *opt_str = "b:p:hq";

    const char *_input_bam;
    const char *_prefix;
    bool _get_qual = 0; // default false

    // int getopt_long (int argc, char *const *argv, const char *shortopts, const struct option *longopts, int *indexptr)
    while ((c = getopt_long(argc, argv, opt_str, long_options, &long_idx)) != -1) {
        // cout << "test" << endl;
        // cout << c << endl;
        if (c == 'b') {
            _input_bam = optarg;
        } else if (c == 'p') {
            _prefix = optarg;
        } else if (c == 'h') {
            usage();
            return 0;
        } else if (c == 'q') {
            _get_qual = 1;
        }
        else {
            std::cerr << "Invalid Arguments.";
        }
    }

    bam_stats(_input_bam, _prefix, _get_qual);

    return 0;
}

