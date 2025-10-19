#ifndef DNA_H
#define DNA_H

#include "postgres.h"
#include "fmgr.h"
#include "utils/builtins.h"
#include "libpq/pqformat.h"
#include "access/hash.h"
#include "access/spgist.h"
#include "utils/array.h"
#include "utils/lsyscache.h"
#include "catalog/pg_type.h"
#include "iupac.h"

/* DNA sequence type */
typedef struct
{
    int32 vl_len_;      /* Variable length header */
    char data[FLEXIBLE_ARRAY_MEMBER]; /* DNA sequence data */
} dna;

/* K-mer type */
typedef struct
{
    int32 vl_len_;      /* Variable length header */
    int32 k;            /* K-mer length */
    char data[FLEXIBLE_ARRAY_MEMBER]; /* K-mer sequence data */
} kmer;

/* Quality K-mer type (k-mer with quality scores) */
typedef struct
{
    int32 vl_len_;      /* Variable length header */
    int32 k;            /* K-mer length */
    char sequence[FLEXIBLE_ARRAY_MEMBER]; /* Sequence followed by quality scores */
} qkmer;

/* Macros for accessing DNA data */
#define DatumGetDnaP(X)         ((dna *) PG_DETOAST_DATUM(X))
#define DatumGetKmerP(X)        ((kmer *) PG_DETOAST_DATUM(X))
#define DatumGetQKmerP(X)       ((qkmer *) PG_DETOAST_DATUM(X))

#define PG_GETARG_DNA_P(n)      DatumGetDnaP(PG_GETARG_DATUM(n))
#define PG_GETARG_KMER_P(n)     DatumGetKmerP(PG_GETARG_DATUM(n))
#define PG_GETARG_QKMER_P(n)    DatumGetQKmerP(PG_GETARG_DATUM(n))

#define PG_RETURN_DNA_P(x)      PG_RETURN_POINTER(x)
#define PG_RETURN_KMER_P(x)     PG_RETURN_POINTER(x)
#define PG_RETURN_QKMER_P(x)    PG_RETURN_POINTER(x)

/* Function declarations */

/* Type input/output functions */
Datum dna_in(PG_FUNCTION_ARGS);
Datum dna_out(PG_FUNCTION_ARGS);
Datum dna_recv(PG_FUNCTION_ARGS);
Datum dna_send(PG_FUNCTION_ARGS);

Datum kmer_in(PG_FUNCTION_ARGS);
Datum kmer_out(PG_FUNCTION_ARGS);
Datum kmer_recv(PG_FUNCTION_ARGS);
Datum kmer_send(PG_FUNCTION_ARGS);

Datum qkmer_in(PG_FUNCTION_ARGS);
Datum qkmer_out(PG_FUNCTION_ARGS);
Datum qkmer_recv(PG_FUNCTION_ARGS);
Datum qkmer_send(PG_FUNCTION_ARGS);

/* Utility functions */
Datum dna_length(PG_FUNCTION_ARGS);
Datum generate_kmers(PG_FUNCTION_ARGS);
Datum dna_complement(PG_FUNCTION_ARGS);
Datum dna_reverse(PG_FUNCTION_ARGS);
Datum dna_reverse_complement(PG_FUNCTION_ARGS);
Datum dna_count(PG_FUNCTION_ARGS);
Datum dna_count_approx(PG_FUNCTION_ARGS);
Datum dna_to_string(PG_FUNCTION_ARGS);
Datum string_to_dna(PG_FUNCTION_ARGS);

/* Operators */
Datum dna_eq(PG_FUNCTION_ARGS);
Datum dna_ne(PG_FUNCTION_ARGS);
Datum dna_lt(PG_FUNCTION_ARGS);
Datum dna_le(PG_FUNCTION_ARGS);
Datum dna_gt(PG_FUNCTION_ARGS);
Datum dna_ge(PG_FUNCTION_ARGS);
Datum dna_cmp(PG_FUNCTION_ARGS);
Datum dna_contains(PG_FUNCTION_ARGS);
Datum dna_contained_by(PG_FUNCTION_ARGS);

Datum kmer_eq(PG_FUNCTION_ARGS);
Datum kmer_ne(PG_FUNCTION_ARGS);
Datum kmer_cmp(PG_FUNCTION_ARGS);

/* Hash support */
Datum dna_hash(PG_FUNCTION_ARGS);
Datum kmer_hash(PG_FUNCTION_ARGS);

/* B-tree support */
Datum dna_btree_cmp(PG_FUNCTION_ARGS);
Datum kmer_btree_cmp(PG_FUNCTION_ARGS);
Datum dna_btree_lt(PG_FUNCTION_ARGS);
Datum dna_btree_le(PG_FUNCTION_ARGS);
Datum dna_btree_gt(PG_FUNCTION_ARGS);
Datum dna_btree_ge(PG_FUNCTION_ARGS);
Datum kmer_btree_lt(PG_FUNCTION_ARGS);
Datum kmer_btree_le(PG_FUNCTION_ARGS);
Datum kmer_btree_gt(PG_FUNCTION_ARGS);
Datum kmer_btree_ge(PG_FUNCTION_ARGS);
Datum dna_sortsupport(PG_FUNCTION_ARGS);
Datum kmer_sortsupport(PG_FUNCTION_ARGS);
Datum dna_in_range(PG_FUNCTION_ARGS);

/* Function prototypes */
Datum dna_gc_content(PG_FUNCTION_ARGS);
Datum dna_count_nucleotide(PG_FUNCTION_ARGS);
Datum dna_find_subsequence(PG_FUNCTION_ARGS);
Datum dna_is_palindrome(PG_FUNCTION_ARGS);
Datum dna_translate(PG_FUNCTION_ARGS);
Datum dna_sliding_gc(PG_FUNCTION_ARGS);
Datum dna_overlap(PG_FUNCTION_ARGS);
Datum dna_similarity(PG_FUNCTION_ARGS);
Datum dna_hash_extended(PG_FUNCTION_ARGS);
Datum kmer_hash_extended(PG_FUNCTION_ARGS);
Datum qkmer_hash(PG_FUNCTION_ARGS);
Datum qkmer_hash_extended(PG_FUNCTION_ARGS);
Datum dna_kmer_hashes(PG_FUNCTION_ARGS);
Datum qkmer_avg_quality(PG_FUNCTION_ARGS);
Datum qkmer_min_quality(PG_FUNCTION_ARGS);
Datum qkmer_filter_quality(PG_FUNCTION_ARGS);

/* SP-GiST support */
Datum spgist_kmer_config(PG_FUNCTION_ARGS);
Datum spgist_kmer_choose(PG_FUNCTION_ARGS);
Datum spgist_kmer_picksplit(PG_FUNCTION_ARGS);
Datum spgist_kmer_inner_consistent(PG_FUNCTION_ARGS);
Datum spgist_kmer_leaf_consistent(PG_FUNCTION_ARGS);

/* Internal utility functions */
int dna_compare_internal(const dna *a, const dna *b);
int kmer_compare_internal(const kmer *a, const kmer *b);
char *dna_get_str(const dna *d);
char *kmer_get_str(const kmer *k);
int dna_get_length(const dna *d);
int kmer_get_k(const kmer *k);

#endif /* DNA_H */
