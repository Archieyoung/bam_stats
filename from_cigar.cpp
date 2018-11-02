#include <iostream>

#include "from_cigar.h"


int32_t from_cigar::get_query_length_cigar() {
    // local copy of cigar_array pointer
    const uint32_t *_ca_prt = ca_prt;
    
    int32_t query_l = 0;
    // loop through cigar array
    for (uint32_t i = 0; i < ca_len; ++i) {
        // bit 2 set if the cigar operation consumes the reference
        const int c_op = bam_cigar_op(*(_ca_prt+i));
        const int c_type = bam_cigar_type(c_op);
        // std::cout << "test: " << bam_cigar_opchr(*(_ca_prt+i)) << std::endl;
        if (c_type&1 || c_op == BAM_CHARD_CLIP) {
            const int c_oplen = bam_cigar_oplen(*(_ca_prt+i));
            query_l += c_oplen;
        }
    }
    // std::cerr << query_l << std::endl;
    return query_l;
}


uint32_t from_cigar::get_query_start(const uint16_t &_flag,
    int32_t &_l_query) {
    const int c_op = bam_cigar_op(*ca_prt);
    const int c_oplen = bam_cigar_oplen(*ca_prt);
    // if the read is reversed
    if (_flag&BAM_FREVERSE) {
        // is the start seq of the query clipped?
        if (c_op == BAM_CSOFT_CLIP || c_op == BAM_CHARD_CLIP) {        
            return _l_query - 1 - c_oplen;
        } else {
            return _l_query - 1;
        }
    } else {
        // is the start seq of the query clipped?
        if (c_op == BAM_CSOFT_CLIP || c_op == BAM_CHARD_CLIP) {
            return c_oplen;
        } else {
            return 0;
        }
    }
}

uint32_t from_cigar::get_query_end(const uint16_t &_flag,
    int32_t &_l_query) {
    const int c_op = bam_cigar_op(*(ca_prt+ca_len-1));
    const int c_oplen = bam_cigar_oplen(*(ca_prt+ca_len-1));
    // if the read is reversed
    if (_flag&BAM_FREVERSE) {
        // is the end seq of the query clipped?
        if (c_op == BAM_CSOFT_CLIP || c_op == BAM_CHARD_CLIP) {        
            return c_oplen;
        } else {
            return 0;
        }
    } else {
        // is the end seq of the query clipped?
        if (c_op == BAM_CSOFT_CLIP || c_op == BAM_CHARD_CLIP) {
            return _l_query - 1 - c_oplen;
        } else {
            return _l_query - 1;
        }
    }
}

uint32_t from_cigar::get_ref_span() {
    // local copy of cigar_array pointer
    const uint32_t *_ca_prt = ca_prt;
    uint32_t ref_span = 0;
    // loop through cigar array
    for (uint32_t i = 0; i < ca_len; ++i) {
        // bit 2 set if the cigar operation consumes the reference
        const int c_op = bam_cigar_op(*(_ca_prt+i));
        const int c_type = bam_cigar_type(c_op);
        // std::cout << "test: " << bam_cigar_opchr(*(_ca_prt+i)) << std::endl;
        if (c_type&2) {
            const int c_oplen = bam_cigar_oplen(*(_ca_prt+i));
            ref_span += c_oplen;
        }
    }
    return ref_span;
}

uint32_t from_cigar::get_ref_start(const uint32_t &_pos) {
    return _pos; // 0-based leftmost coordinate
}

uint32_t from_cigar::get_ref_end(const uint32_t &_pos) {
    return _pos + get_ref_span() - 1; // 0-based rightmost coordinate
}

/* test
int main(int argc, char const *argv[])
{
    // argument check
    if (argc < 2) {
        std::cerr << "input error!" << std::endl;
        return -1;
    }

    // open bam
    samFile *fp = sam_open(argv[1], "r");
    bam_hdr_t *h = sam_hdr_read(fp);
    bam1_t *b;

    b = bam_init1();


    int ret;
    // read record one by one
    while ((ret = sam_read1(fp, h, b) > 0)) {
        const uint32_t _cigar_array_len = b->core.n_cigar;
        const uint16_t bam_flag = b->core.flag;
        const uint32_t l_query = b->core.l_qseq;
        const uint32_t pos = b->core.pos;
        const char *ref_name = h->target_name[b->core.tid];
        // std::cout << _cigar_array_len << std::endl;
        // from_cigar constructor
        from_cigar _cp(bam_get_cigar(b), _cigar_array_len);
        std::cout << ref_name
            << "\t"
            << _cp.get_ref_start(pos)
            << "\t"
            << _cp.get_ref_end(pos)
            << "\t"
            << _cp.get_query_start(bam_flag, l_query)
            << "\t" 
            <<  _cp.get_query_end(bam_flag, l_query)
            << std::endl;
    }
    
    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fp);
    return 0;
}
*/
