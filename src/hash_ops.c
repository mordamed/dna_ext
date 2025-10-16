#include "dna.h"

/*
 * Hash support functions for DNA and K-mer types
 * Enables use of hash indexes and hash joins
 */

/*
 * Hash function for DNA type
 */
PG_FUNCTION_INFO_V1(dna_hash);
Datum
dna_hash(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int len = VARSIZE_ANY_EXHDR(d);
    
    /* Use PostgreSQL's hash_any function */
    PG_RETURN_UINT32(hash_any((unsigned char *) d->data, len));
}

/*
 * Extended hash function for DNA type (for hash partitioning)
 */
PG_FUNCTION_INFO_V1(dna_hash_extended);
Datum
dna_hash_extended(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    uint64 seed = PG_GETARG_INT64(1);
    int len = VARSIZE_ANY_EXHDR(d);
    
    /* Use PostgreSQL's hash_any_extended function */
    PG_RETURN_UINT64(hash_any_extended((unsigned char *) d->data, len, seed));
}

/*
 * Hash function for K-mer type
 */
PG_FUNCTION_INFO_V1(kmer_hash);
Datum
kmer_hash(PG_FUNCTION_ARGS)
{
    kmer *k = PG_GETARG_KMER_P(0);
    uint32 hash_value;
    
    /* Hash both the k value and the sequence data */
    hash_value = hash_any((unsigned char *) &k->k, sizeof(int32));
    hash_value ^= hash_any((unsigned char *) k->data, k->k);
    
    PG_RETURN_UINT32(hash_value);
}

/*
 * Extended hash function for K-mer type
 */
PG_FUNCTION_INFO_V1(kmer_hash_extended);
Datum
kmer_hash_extended(PG_FUNCTION_ARGS)
{
    kmer *k = PG_GETARG_KMER_P(0);
    uint64 seed = PG_GETARG_INT64(1);
    uint64 hash_value;
    
    /* Hash both the k value and the sequence data with seed */
    hash_value = hash_any_extended((unsigned char *) &k->k, sizeof(int32), seed);
    hash_value ^= hash_any_extended((unsigned char *) k->data, k->k, seed);
    
    PG_RETURN_UINT64(hash_value);
}

/*
 * Hash function for QKmer type
 */
PG_FUNCTION_INFO_V1(qkmer_hash);
Datum
qkmer_hash(PG_FUNCTION_ARGS)
{
    qkmer *qk = PG_GETARG_QKMER_P(0);
    uint32 hash_value;
    
    /* Hash the k value and both sequence and quality data */
    hash_value = hash_any((unsigned char *) &qk->k, sizeof(int32));
    hash_value ^= hash_any((unsigned char *) qk->sequence, qk->k * 2);
    
    PG_RETURN_UINT32(hash_value);
}

/*
 * Extended hash function for QKmer type
 */
PG_FUNCTION_INFO_V1(qkmer_hash_extended);
Datum
qkmer_hash_extended(PG_FUNCTION_ARGS)
{
    qkmer *qk = PG_GETARG_QKMER_P(0);
    uint64 seed = PG_GETARG_INT64(1);
    uint64 hash_value;
    
    /* Hash the k value and both sequence and quality data with seed */
    hash_value = hash_any_extended((unsigned char *) &qk->k, sizeof(int32), seed);
    hash_value ^= hash_any_extended((unsigned char *) qk->sequence, qk->k * 2, seed);
    
    PG_RETURN_UINT64(hash_value);
}

/*
 * Rolling hash function for efficient k-mer sliding window
 * Useful for generating k-mers from large sequences
 */
static uint32
rolling_hash_add_nucleotide(uint32 hash, char nucleotide, int k)
{
    uint32 base_value = nucleotide_to_int(nucleotide);
    
    if (base_value == -1)
        return hash; /* Invalid nucleotide, return unchanged hash */
    
    /* Simple rolling hash: shift left and add new base */
    return (hash << 2) | base_value;
}

static uint32
rolling_hash_remove_nucleotide(uint32 hash, char nucleotide, int k)
{
    uint32 base_value = nucleotide_to_int(nucleotide);
    uint32 mask = (1 << (k * 2)) - 1; /* Mask for k nucleotides */
    
    if (base_value == -1)
        return hash; /* Invalid nucleotide, return unchanged hash */
    
    /* Remove the leftmost nucleotide */
    return hash & mask;
}

/*
 * Generate hash values for all k-mers in a DNA sequence
 * Returns an array of hash values
 */
PG_FUNCTION_INFO_V1(dna_kmer_hashes);
Datum
dna_kmer_hashes(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int32 k = PG_GETARG_INT32(1);
    int seq_len = dna_get_length(d);
    char *seq = dna_get_str(d);
    ArrayType *result;
    Datum *elems;
    int num_kmers;
    uint32 rolling_hash = 0;
    int i;
    
    if (k <= 0 || k > seq_len)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k must be between 1 and sequence length")));
    
    if (k > 16) /* Limit to prevent overflow in rolling hash */
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k must be <= 16 for rolling hash")));
    
    num_kmers = seq_len - k + 1;
    elems = (Datum *) palloc(num_kmers * sizeof(Datum));
    
    /* Initialize rolling hash with first k-mer */
    for (i = 0; i < k; i++)
    {
        rolling_hash = rolling_hash_add_nucleotide(rolling_hash, seq[i], k);
    }
    elems[0] = UInt32GetDatum(rolling_hash);
    
    /* Generate remaining k-mer hashes using rolling hash */
    for (i = 1; i < num_kmers; i++)
    {
        /* Remove the leftmost nucleotide and add the new rightmost nucleotide */
        rolling_hash = rolling_hash_remove_nucleotide(rolling_hash, seq[i - 1], k);
        rolling_hash = rolling_hash_add_nucleotide(rolling_hash, seq[i + k - 1], k);
        
        elems[i] = UInt32GetDatum(rolling_hash);
    }
    
    result = construct_array(elems, num_kmers, INT4OID, 4, true, 'i');
    
    pfree(seq);
    pfree(elems);
    
    PG_RETURN_ARRAYTYPE_P(result);
}