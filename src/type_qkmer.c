#include "dna.h"
#include <string.h>
#include <ctype.h>

/*
 * Quality K-mer type input/output functions
 * QKmer stores a k-mer with associated quality scores
 */

/*
 * QKmer input function
 * Format: "ACGT:!\"#$" (sequence:quality)
 */
PG_FUNCTION_INFO_V1(qkmer_in);
Datum
qkmer_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    char *colon_pos = strchr(str, ':');
    qkmer *result;
    int seq_len;
    int i;
    
    if (!colon_pos)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("qkmer format must be sequence:quality")));
    
    seq_len = colon_pos - str;
    
    if (seq_len == 0)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("qkmer sequence cannot be empty")));
    
    if (strlen(colon_pos + 1) != seq_len)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("quality string length must match sequence length")));
    
    /* Validate sequence */
    for (i = 0; i < seq_len; i++)
    {
        char c = toupper(str[i]);
        if (!is_valid_nucleotide(c))
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("invalid nucleotide character in qkmer: %c", str[i])));
    }
    
    /* Allocate and populate result */
    result = (qkmer *) palloc(VARHDRSZ + sizeof(int32) + seq_len * 2);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + seq_len * 2);
    result->k = seq_len;
    
    /* Copy sequence */
    for (i = 0; i < seq_len; i++)
    {
        result->sequence[i] = toupper(str[i]);
    }
    
    /* Copy quality scores */
    memcpy(result->sequence + seq_len, colon_pos + 1, seq_len);
    
    PG_RETURN_QKMER_P(result);
}

/*
 * QKmer output function
 */
PG_FUNCTION_INFO_V1(qkmer_out);
Datum
qkmer_out(PG_FUNCTION_ARGS)
{
    qkmer *qk = PG_GETARG_QKMER_P(0);
    char *result = palloc(qk->k * 2 + 2); /* sequence + ':' + quality + '\0' */
    
    memcpy(result, qk->sequence, qk->k);
    result[qk->k] = ':';
    memcpy(result + qk->k + 1, qk->sequence + qk->k, qk->k);
    result[qk->k * 2 + 1] = '\0';
    
    PG_RETURN_CSTRING(result);
}

/*
 * QKmer binary receive function
 */
PG_FUNCTION_INFO_V1(qkmer_recv);
Datum
qkmer_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
    qkmer *result;
    int32 k;
    
    k = pq_getmsgint(buf, 4);
    if (k <= 0)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_BINARY_REPRESENTATION),
                 errmsg("invalid qkmer length in binary representation")));
    
    result = (qkmer *) palloc(VARHDRSZ + sizeof(int32) + k * 2);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + k * 2);
    result->k = k;
    
    pq_copymsgbytes(buf, result->sequence, k * 2);
    
    PG_RETURN_QKMER_P(result);
}

/*
 * QKmer binary send function
 */
PG_FUNCTION_INFO_V1(qkmer_send);
Datum
qkmer_send(PG_FUNCTION_ARGS)
{
    qkmer *qk = PG_GETARG_QKMER_P(0);
    StringInfoData buf;
    
    pq_begintypsend(&buf);
    pq_sendint32(&buf, qk->k);
    pq_sendbytes(&buf, qk->sequence, qk->k * 2);
    
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

/*
 * Get average quality score for a qkmer
 */
PG_FUNCTION_INFO_V1(qkmer_avg_quality);
Datum
qkmer_avg_quality(PG_FUNCTION_ARGS)
{
    qkmer *qk = PG_GETARG_QKMER_P(0);
    double sum = 0.0;
    int i;
    
    for (i = 0; i < qk->k; i++)
    {
        /* Convert Phred quality to numeric (assuming Phred+33 encoding) */
        sum += (double)(qk->sequence[qk->k + i] - 33);
    }
    
    PG_RETURN_FLOAT8(sum / qk->k);
}

/*
 * Get minimum quality score for a qkmer
 */
PG_FUNCTION_INFO_V1(qkmer_min_quality);
Datum
qkmer_min_quality(PG_FUNCTION_ARGS)
{
    qkmer *qk = PG_GETARG_QKMER_P(0);
    int min_qual = qk->sequence[qk->k] - 33;
    int i;
    
    for (i = 1; i < qk->k; i++)
    {
        int qual = qk->sequence[qk->k + i] - 33;
        if (qual < min_qual)
            min_qual = qual;
    }
    
    PG_RETURN_INT32(min_qual);
}

/*
 * Filter qkmer by minimum quality threshold
 */
PG_FUNCTION_INFO_V1(qkmer_filter_quality);
Datum
qkmer_filter_quality(PG_FUNCTION_ARGS)
{
    qkmer *qk = PG_GETARG_QKMER_P(0);
    int min_threshold = PG_GETARG_INT32(1);
    int i;
    
    for (i = 0; i < qk->k; i++)
    {
        int qual = qk->sequence[qk->k + i] - 33;
        if (qual < min_threshold)
            PG_RETURN_BOOL(false);
    }
    
    PG_RETURN_BOOL(true);
}