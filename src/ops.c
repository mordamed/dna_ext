#include "dna.h"

/*
 * DNA and K-mer operators
 * Implements comparison and containment operators
 */

/*
 * DNA equality operator
 */
PG_FUNCTION_INFO_V1(dna_eq);
Datum
dna_eq(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_BOOL(dna_compare_internal(a, b) == 0);
}

/*
 * DNA inequality operator
 */
PG_FUNCTION_INFO_V1(dna_ne);
Datum
dna_ne(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_BOOL(dna_compare_internal(a, b) != 0);
}

/*
 * DNA less than operator
 */
PG_FUNCTION_INFO_V1(dna_lt);
Datum
dna_lt(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_BOOL(dna_compare_internal(a, b) < 0);
}

/*
 * DNA less than or equal operator
 */
PG_FUNCTION_INFO_V1(dna_le);
Datum
dna_le(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_BOOL(dna_compare_internal(a, b) <= 0);
}

/*
 * DNA greater than operator
 */
PG_FUNCTION_INFO_V1(dna_gt);
Datum
dna_gt(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_BOOL(dna_compare_internal(a, b) > 0);
}

/*
 * DNA greater than or equal operator
 */
PG_FUNCTION_INFO_V1(dna_ge);
Datum
dna_ge(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_BOOL(dna_compare_internal(a, b) >= 0);
}

/*
 * DNA comparison function for sorting
 */
PG_FUNCTION_INFO_V1(dna_cmp);
Datum
dna_cmp(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_INT32(dna_compare_internal(a, b));
}

/*
 * DNA contains operator (@>)
 * Returns true if the left DNA sequence contains the right DNA sequence
 */
PG_FUNCTION_INFO_V1(dna_contains);
Datum
dna_contains(PG_FUNCTION_ARGS)
{
    dna *haystack = PG_GETARG_DNA_P(0);
    dna *needle = PG_GETARG_DNA_P(1);
    char *haystack_str = dna_get_str(haystack);
    char *needle_str = dna_get_str(needle);
    bool result;
    
    result = (strstr(haystack_str, needle_str) != NULL);
    
    pfree(haystack_str);
    pfree(needle_str);
    
    PG_RETURN_BOOL(result);
}

/*
 * DNA contained by operator (<@)
 * Returns true if the left DNA sequence is contained in the right DNA sequence
 */
PG_FUNCTION_INFO_V1(dna_contained_by);
Datum
dna_contained_by(PG_FUNCTION_ARGS)
{
    dna *needle = PG_GETARG_DNA_P(0);
    dna *haystack = PG_GETARG_DNA_P(1);
    char *haystack_str = dna_get_str(haystack);
    char *needle_str = dna_get_str(needle);
    bool result;
    
    result = (strstr(haystack_str, needle_str) != NULL);
    
    pfree(haystack_str);
    pfree(needle_str);
    
    PG_RETURN_BOOL(result);
}

/*
 * DNA overlap operator (&&)
 * Returns true if two DNA sequences have any common subsequence
 */
PG_FUNCTION_INFO_V1(dna_overlap);
Datum
dna_overlap(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    char *seq_a = dna_get_str(a);
    char *seq_b = dna_get_str(b);
    int len_a = dna_get_length(a);
    int len_b = dna_get_length(b);
    bool result = false;
    int i, j, k;
    
    /* Check for common subsequences of length >= 3 */
    for (i = 0; i <= len_a - 3 && !result; i++)
    {
        for (j = 0; j <= len_b - 3 && !result; j++)
        {
            for (k = 3; k <= len_a - i && k <= len_b - j; k++)
            {
                if (memcmp(seq_a + i, seq_b + j, k) == 0)
                {
                    result = true;
                    break;
                }
            }
        }
    }
    
    pfree(seq_a);
    pfree(seq_b);
    
    PG_RETURN_BOOL(result);
}

/*
 * K-mer equality operator
 */
PG_FUNCTION_INFO_V1(kmer_eq);
Datum
kmer_eq(PG_FUNCTION_ARGS)
{
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    
    PG_RETURN_BOOL(kmer_compare_internal(a, b) == 0);
}

/*
 * K-mer inequality operator
 */
PG_FUNCTION_INFO_V1(kmer_ne);
Datum
kmer_ne(PG_FUNCTION_ARGS)
{
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    
    PG_RETURN_BOOL(kmer_compare_internal(a, b) != 0);
}

/*
 * K-mer comparison function for sorting
 */
PG_FUNCTION_INFO_V1(kmer_cmp);
Datum
kmer_cmp(PG_FUNCTION_ARGS)
{
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    
    PG_RETURN_INT32(kmer_compare_internal(a, b));
}

/*
 * DNA similarity operator (^@)
 * Returns similarity score between two DNA sequences
 */
PG_FUNCTION_INFO_V1(dna_similarity);
Datum
dna_similarity(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    char *seq_a = dna_get_str(a);
    char *seq_b = dna_get_str(b);
    int len_a = dna_get_length(a);
    int len_b = dna_get_length(b);
    int matches = 0;
    int max_len = (len_a > len_b) ? len_a : len_b;
    int min_len = (len_a < len_b) ? len_a : len_b;
    double similarity;
    int i;
    
    /* Count matching positions */
    for (i = 0; i < min_len; i++)
    {
        if (seq_a[i] == seq_b[i])
            matches++;
    }
    
    /* Calculate similarity as ratio of matches to maximum length */
    similarity = (max_len > 0) ? (double)matches / max_len : 0.0;
    
    pfree(seq_a);
    pfree(seq_b);
    
    PG_RETURN_FLOAT8(similarity);
}