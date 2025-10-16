#include "dna.h"
#include <string.h>
#include <ctype.h>

/*
 * DNA type input/output functions
 */

/*
 * DNA input function
 */
PG_FUNCTION_INFO_V1(dna_in);
Datum
dna_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    int len = strlen(str);
    dna *result;
    int i;
    
    /* Validate input */
    for (i = 0; i < len; i++)
    {
        char c = toupper(str[i]);
        if (!is_valid_nucleotide(c))
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("invalid nucleotide character: %c", str[i])));
    }
    
    /* Allocate and populate result */
    result = (dna *) palloc(VARHDRSZ + len);
    SET_VARSIZE(result, VARHDRSZ + len);
    
    for (i = 0; i < len; i++)
    {
        result->data[i] = toupper(str[i]);
    }
    
    PG_RETURN_DNA_P(result);
}

/*
 * DNA output function
 */
PG_FUNCTION_INFO_V1(dna_out);
Datum
dna_out(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int len = VARSIZE_ANY_EXHDR(d);
    char *result = palloc(len + 1);
    
    memcpy(result, d->data, len);
    result[len] = '\0';
    
    PG_RETURN_CSTRING(result);
}

/*
 * DNA binary receive function
 */
PG_FUNCTION_INFO_V1(dna_recv);
Datum
dna_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
    dna *result;
    int len;
    
    len = pq_getmsgint(buf, 4);
    if (len < 0)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_BINARY_REPRESENTATION),
                 errmsg("invalid length in DNA binary representation")));
    
    result = (dna *) palloc(VARHDRSZ + len);
    SET_VARSIZE(result, VARHDRSZ + len);
    
    pq_copymsgbytes(buf, result->data, len);
    
    PG_RETURN_DNA_P(result);
}

/*
 * DNA binary send function
 */
PG_FUNCTION_INFO_V1(dna_send);
Datum
dna_send(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    StringInfoData buf;
    int len = VARSIZE_ANY_EXHDR(d);
    
    pq_begintypsend(&buf);
    pq_sendint32(&buf, len);
    pq_sendbytes(&buf, d->data, len);
    
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

/*
 * DNA comparison function
 */
int
dna_compare_internal(const dna *a, const dna *b)
{
    int len_a = VARSIZE_ANY_EXHDR(a);
    int len_b = VARSIZE_ANY_EXHDR(b);
    int result;
    
    result = memcmp(a->data, b->data, Min(len_a, len_b));
    
    if (result != 0)
        return result;
    
    /* If prefixes are equal, shorter string is less */
    if (len_a < len_b)
        return -1;
    else if (len_a > len_b)
        return 1;
    else
        return 0;
}