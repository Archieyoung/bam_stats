#ifndef FROM_CIGAR_H
#define FROM_CIGAR_H

#include "htslib/sam.h"

class from_cigar {
    public:
    from_cigar() = default;
    from_cigar(const uint32_t *cigar_array_prt, const uint32_t cigar_array_len):
        ca_prt(cigar_array_prt), ca_len(cigar_array_len) {};

    /*
     * return query lenght
    */
    int32_t get_query_length_cigar();
    /*
    parameter: bam flag, length of the query
    return: start position on query, position is "0" based
    */
    uint32_t get_query_start(const uint16_t &_flag, int32_t &_l_query);

    /*
    parameter: bam flag, length of the query
    return: end position on query, position is "0" based
    */
    uint32_t get_query_end(const uint16_t &_flag, int32_t &_l_query);

    /*
    parameter: start position on reference
    return: start position on reference, position is "0" based
    */
    uint32_t get_ref_start(const uint32_t &_pos); // start from "0"

    /*
    return span on reference
    */
    uint32_t get_ref_span();

    /*
    parameter: bam flag
    return: end position on reference, position is "0" based
    */
    uint32_t get_ref_end(const uint32_t &_pos); // start from "0"

    private:
    const uint32_t *ca_prt;
    const uint32_t ca_len;
};

#endif
