#include "dna.h"
#include "utils/sortsupport.h"

/* Forward declarations for static functions */
static int dna_fastcmp(Datum x, Datum y, SortSupport ssup);
static int kmer_fastcmp(Datum x, Datum y, SortSupport ssup);

/*
 * B-tree support functions for DNA and K-mer types
 * Enables creation of B-tree indexes for ordering and range queries
 */

/*
 * B-tree comparison function for DNA type
 */
PG_FUNCTION_INFO_V1(dna_btree_cmp);
Datum
dna_btree_cmp(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_INT32(dna_compare_internal(a, b));
}

/*
 * B-tree less than function for DNA type
 */
PG_FUNCTION_INFO_V1(dna_btree_lt);
Datum
dna_btree_lt(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_BOOL(dna_compare_internal(a, b) < 0);
}

/*
 * B-tree less than or equal function for DNA type
 */
PG_FUNCTION_INFO_V1(dna_btree_le);
Datum
dna_btree_le(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_BOOL(dna_compare_internal(a, b) <= 0);
}

/*
 * B-tree greater than function for DNA type
 */
PG_FUNCTION_INFO_V1(dna_btree_gt);
Datum
dna_btree_gt(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_BOOL(dna_compare_internal(a, b) > 0);
}

/*
 * B-tree greater than or equal function for DNA type
 */
PG_FUNCTION_INFO_V1(dna_btree_ge);
Datum
dna_btree_ge(PG_FUNCTION_ARGS)
{
    dna *a = PG_GETARG_DNA_P(0);
    dna *b = PG_GETARG_DNA_P(1);
    
    PG_RETURN_BOOL(dna_compare_internal(a, b) >= 0);
}

/*
 * B-tree comparison function for K-mer type
 */
PG_FUNCTION_INFO_V1(kmer_btree_cmp);
Datum
kmer_btree_cmp(PG_FUNCTION_ARGS)
{
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    
    PG_RETURN_INT32(kmer_compare_internal(a, b));
}

/*
 * B-tree less than function for K-mer type
 */
PG_FUNCTION_INFO_V1(kmer_btree_lt);
Datum
kmer_btree_lt(PG_FUNCTION_ARGS)
{
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    
    PG_RETURN_BOOL(kmer_compare_internal(a, b) < 0);
}

/*
 * B-tree less than or equal function for K-mer type
 */
PG_FUNCTION_INFO_V1(kmer_btree_le);
Datum
kmer_btree_le(PG_FUNCTION_ARGS)
{
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    
    PG_RETURN_BOOL(kmer_compare_internal(a, b) <= 0);
}

/*
 * B-tree greater than function for K-mer type
 */
PG_FUNCTION_INFO_V1(kmer_btree_gt);
Datum
kmer_btree_gt(PG_FUNCTION_ARGS)
{
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    
    PG_RETURN_BOOL(kmer_compare_internal(a, b) > 0);
}

/*
 * B-tree greater than or equal function for K-mer type
 */
PG_FUNCTION_INFO_V1(kmer_btree_ge);
Datum
kmer_btree_ge(PG_FUNCTION_ARGS)
{
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    
    PG_RETURN_BOOL(kmer_compare_internal(a, b) >= 0);
}

/*
 * DNA sortsupport function for optimized sorting
 */
PG_FUNCTION_INFO_V1(dna_sortsupport);
Datum
dna_sortsupport(PG_FUNCTION_ARGS)
{
    SortSupport ssup = (SortSupport) PG_GETARG_POINTER(0);
    
    ssup->comparator = dna_fastcmp;
    ssup->ssup_extra = NULL;
    
    PG_RETURN_VOID();
}

/*
 * Fast comparison function for DNA sorting
 */
static int
dna_fastcmp(Datum x, Datum y, SortSupport ssup)
{
    dna *a = DatumGetDnaP(x);
    dna *b = DatumGetDnaP(y);
    
    return dna_compare_internal(a, b);
}

/*
 * K-mer sortsupport function for optimized sorting
 */
PG_FUNCTION_INFO_V1(kmer_sortsupport);
Datum
kmer_sortsupport(PG_FUNCTION_ARGS)
{
    SortSupport ssup = (SortSupport) PG_GETARG_POINTER(0);
    
    ssup->comparator = kmer_fastcmp;
    ssup->ssup_extra = NULL;
    
    PG_RETURN_VOID();
}

/*
 * Fast comparison function for K-mer sorting
 */
static int
kmer_fastcmp(Datum x, Datum y, SortSupport ssup)
{
    kmer *a = DatumGetKmerP(x);
    kmer *b = DatumGetKmerP(y);
    
    return kmer_compare_internal(a, b);
}

/*
 * DNA range query support functions
 */

/*
 * Check if DNA sequence is within a lexicographic range
 */
PG_FUNCTION_INFO_V1(dna_in_range);
Datum
dna_in_range(PG_FUNCTION_ARGS)
{
    dna *sequence = PG_GETARG_DNA_P(0);
    dna *lower_bound = PG_GETARG_DNA_P(1);
    dna *upper_bound = PG_GETARG_DNA_P(2);
    bool lower_inclusive = PG_GETARG_BOOL(3);
    bool upper_inclusive = PG_GETARG_BOOL(4);
    int lower_cmp, upper_cmp;
    
    lower_cmp = dna_compare_internal(sequence, lower_bound);
    upper_cmp = dna_compare_internal(sequence, upper_bound);
    
    if (lower_inclusive)
    {
        if (lower_cmp < 0)
            PG_RETURN_BOOL(false);
    }
    else
    {
        if (lower_cmp <= 0)
            PG_RETURN_BOOL(false);
    }
    
    if (upper_inclusive)
    {
        if (upper_cmp > 0)
            PG_RETURN_BOOL(false);
    }
    else
    {
        if (upper_cmp >= 0)
            PG_RETURN_BOOL(false);
    }
    
    PG_RETURN_BOOL(true);
}