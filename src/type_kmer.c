#include "dna.h"
#include <string.h>
#include <ctype.h>

/*
 * K-mer type input/output functions
 */

/*
 * K-mer input function
 * Format: "ACGT"
 */
PG_FUNCTION_INFO_V1(kmer_in);
Datum
kmer_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    int len = strlen(str);
    kmer *result;
    int i;
    
    if (len == 0)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("k-mer cannot be empty")));
    
    /* Validate input */
    for (i = 0; i < len; i++)
    {
        char c = toupper(str[i]);
        if (!is_valid_nucleotide(c))
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("invalid nucleotide character in k-mer: %c", str[i])));
    }
    
    /* Allocate and populate result */
    result = (kmer *) palloc(VARHDRSZ + sizeof(int32) + len);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + len);
    result->k = len;
    
    for (i = 0; i < len; i++)
    {
        result->data[i] = toupper(str[i]);
    }
    
    PG_RETURN_KMER_P(result);
}

/*
 * K-mer output function
 */
PG_FUNCTION_INFO_V1(kmer_out);
Datum
kmer_out(PG_FUNCTION_ARGS)
{
    kmer *k = PG_GETARG_KMER_P(0);
    char *result = palloc(k->k + 1);
    
    memcpy(result, k->data, k->k);
    result[k->k] = '\0';
    
    PG_RETURN_CSTRING(result);
}

/*
 * K-mer binary receive function
 */
PG_FUNCTION_INFO_V1(kmer_recv);
Datum
kmer_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
    kmer *result;
    int32 k;
    
    k = pq_getmsgint(buf, 4);
    if (k <= 0)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_BINARY_REPRESENTATION),
                 errmsg("invalid k-mer length in binary representation")));
    
    result = (kmer *) palloc(VARHDRSZ + sizeof(int32) + k);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + k);
    result->k = k;
    
    pq_copymsgbytes(buf, result->data, k);
    
    PG_RETURN_KMER_P(result);
}

/*
 * K-mer binary send function
 */
PG_FUNCTION_INFO_V1(kmer_send);
Datum
kmer_send(PG_FUNCTION_ARGS)
{
    kmer *k = PG_GETARG_KMER_P(0);
    StringInfoData buf;
    
    pq_begintypsend(&buf);
    pq_sendint32(&buf, k->k);
    pq_sendbytes(&buf, k->data, k->k);
    
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

/*
 * Get k-mer as C string
 */
char *
kmer_get_str(const kmer *k)
{
    char *result = palloc(k->k + 1);
    
    memcpy(result, k->data, k->k);
    result[k->k] = '\0';
    
    return result;
}

/*
 * Get k-mer length
 */
int
kmer_get_k(const kmer *k)
{
    return k->k;
}

/*
 * K-mer comparison function
 */
int
kmer_compare_internal(const kmer *a, const kmer *b)
{
    int result;
    
    /* First compare the k values */
    if (a->k < b->k)
        return -1;
    else if (a->k > b->k)
        return 1;
    
    /* If k values are equal, compare sequences */
    result = memcmp(a->data, b->data, a->k);
    
    return result;
}