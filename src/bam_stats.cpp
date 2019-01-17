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

int get_fragment_qual(bam1_t *b, int32_t &l_query,
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

bool is_hard_clip(const uint32_t *cigar_array_prt, const uint32_t cigar_array_len, const uint16_t bam_flag) {
    if (is_unmapped(bam_flag)) {
        return false;
    }
    const int c_op_s = bam_cigar_op(*cigar_array_prt);
    const int c_op_e = bam_cigar_op(*(cigar_array_prt + cigar_array_len - 1));
    if (c_op_s == BAM_CHARD_CLIP || c_op_e == BAM_CHARD_CLIP) {
        return true;
    } else {
        return false;
    }
}

/*
compute gap compressed identity. https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
perl -ane 'if(/NM:i:(\d+)/){$n=$1;$l=0;$l+=$1 while/(\d+)[MID]/g;print(($l-$n)/$l,"\n")}'
identity = (NM - D - I + O)/(M + O), remove gap extend but keep gap open, when calculate identity.
NM is NM tag in sam file, which represent not-match bases. O is number of gap open.
*/
float gap_compressed_identity(const int32_t &nm, const uint32_t *ca_prt,
    const uint32_t &ca_len)
{
    uint32_t l_m = 0;
    uint32_t l_d = 0;
    uint32_t l_i = 0;
    uint32_t l_o = 0;
    for (uint32_t i = 0; i < ca_len; ++i) {
        const int c_op = bam_cigar_op(*(ca_prt+i));
        const int c_oplen = bam_cigar_oplen(*(ca_prt+i));
        if (c_op == BAM_CMATCH) {
            l_m += c_oplen;
        } else if (c_op == BAM_CDEL) {
            l_d += c_oplen;
            ++l_o;
        } else if (c_op == BAM_CINS) {
            l_i += c_oplen;
            ++l_o;
        }
    }
    float _identity = 1 - float(nm - l_d - l_i + l_o)/(l_m + l_o);
    return _identity;
}

int bam_stats(const char *input_bam, const std::string prefix,
    int identity_type, uint8_t min_mq, bool get_qual)
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
        out_hd << "#CHROM\tREF_START\tREF_END\tQUERY_NAME\tQUERY_POS1\tQUERY_POS2\tQUERY_LEN\t"
"MAPPING_QAULITY\tFRAGMENT_MEAN_Q\tFRAGMENT_IDENTITY" << std::endl;
    } else {
        // output table header
        out_hd << "#CHROM\tREF_START\tREF_END\tQUERY_NAME\tQUERY_POS1\tQUERY_POS2\tQUERY_LEN\t"
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
        const uint16_t bam_flag = b->core.flag;
        const uint32_t *_cigar_array_prt = bam_get_cigar(b);
        const uint32_t _cigar_array_len = b->core.n_cigar;
        
        const uint32_t pos = b->core.pos;
        const char *ref_name = h->target_name[b->core.tid];
        const uint16_t mapping_qual = b->core.qual;

        // query name
        const char *query_name = bam_get_qname(b);
        int32_t l_query = b->core.l_qseq;
        
        
        // std::cout << _cigar_array_len << std::endl;
        // from_cigar constructor
        from_cigar _cp(_cigar_array_prt, _cigar_array_len);        


        /*query len from CIGAR*/
        if (l_query <= 0 || is_hard_clip(_cigar_array_prt, _cigar_array_len, bam_flag)) {
            l_query = _cp.get_query_length_cigar();
        }

        if (qlen.find(query_name) == qlen.end()) {
            qlen[query_name] = l_query;
            total_bases += l_query;
        }

        // skip segment with mapping qual less than min_mq
        if (b->core.qual < min_mq) {
            continue;
        }
        
        // skip unmapped record
        if (is_unmapped(bam_flag)) {
            ++unmapped_num;
            continue;
        }

        const uint32_t query_start = _cp.get_query_start(bam_flag, l_query);
        const uint32_t query_end = _cp.get_query_end(bam_flag, l_query);
        // mapped bases is the span of the segment, cliping not included
        uint32_t qspan = abs(query_end - query_start) + 1;
        if (!(bam_flag&BAM_FSECONDARY)) {
            mapped_bases += qspan;
        }

    
        // "NM" tag
        uint8_t *NM_tag;
        // check tag existence
        if (bam_aux_get(b, "NM") != nullptr) {
            NM_tag = bam_aux_get(b, "NM");     
        } else {
            std::cerr << "Can't find \"NM\" tag" << std::endl;
            continue;
        }

        int32_t NM_tag_value;
        // check value get
        NM_tag_value = bam_aux2i(NM_tag);

        uint32_t aligned_len = _cp.get_aligned_len();

        // mapping identity
        float percent_identity = -1;
        if (identity_type == 0) {
            percent_identity = 100*gap_compressed_identity(NM_tag_value, _cigar_array_prt, _cigar_array_len);
        } else if (identity_type == 1) {
            percent_identity = 100*(1 - float(NM_tag_value)/float(aligned_len)); 
        } else {

        }
        if (percent_identity < 0 || percent_identity > 100) {
            std::cerr << "Error, percent of identity out of range for query "
                << query_name << std::endl;
            continue;
        }

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
            << l_query
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
            << l_query
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
    std::cerr << "-h, --help                   print this message and exit\n"
        << "-b, --bam, FILE              input bam file [default: None]\n"
        << "-p, --prefix, FILE           output file prefix [default: None]\n"
        << "-t, --type_of_identity INT   type of identity, 0 for gap_compressed identity and 1 for blast identity [default: 0]\n"
        << "-m, --min_mq                 minimum mapping quality [default: 0]\n"
        << "-q, --get_qual               get fragment sequencing quality [default FALSE]\n"
        << "-V, --version                print version"
        << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc <= 1) {
        usage();
        return 0;
    }

    std::string __version__ = "v0.0.1";

    // long_option array
    static const struct option long_options[] = {
        {"bam", required_argument, 0, 'b'},
        {"out", required_argument, 0, 'o'},
        {"type_of_identity", required_argument, 0, 't'},
        {"get_qual", no_argument, 0, 'q'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'V'}
    };

    int c = 100, long_idx;
    const char *opt_str = "b:p:t:m:hqV";

    const char *_input_bam;
    const char *_prefix;
    bool _get_qual = 0; // default false
    int identity_type = 0; // default gap compressed identity
    uint8_t min_mq = 0;

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
        } else if (c == 't') {
            identity_type = atoi(optarg);
        } else if (c == 'm') {
            min_mq = atoi(optarg);
        } else if (c == 'V') {
            std::cout << "bam_stats version: " << __version__ << std::endl;
            return 0;
        }
        else {
            std::cerr << "Invalid Arguments.";
            std::exit(1);
        }
    }

    bam_stats(_input_bam, _prefix, identity_type, min_mq, _get_qual);

    return 0;
}
