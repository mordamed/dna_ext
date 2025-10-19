#include "dna.h"
#include "iupac.h"
#include <string.h>
#include <ctype.h>
#include "catalog/namespace.h"

/*
 * DNA utility functions
 */

/*
 * Get the length of a DNA sequence
 */
PG_FUNCTION_INFO_V1(dna_length);
Datum
dna_length(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int32 len = VARSIZE_ANY_EXHDR(d);
    
    PG_RETURN_INT32(len);
}

/*
 * Get DNA sequence as C string
 */
char *
dna_get_str(const dna *d)
{
    int len = VARSIZE_ANY_EXHDR(d);
    char *result = palloc(len + 1);
    
    memcpy(result, d->data, len);
    result[len] = '\0';
    
    return result;
}

/*
 * Get the actual length of DNA sequence
 */
int
dna_get_length(const dna *d)
{
    return VARSIZE_ANY_EXHDR(d);
}

/*
 * Validate DNA sequence
 */
static bool
validate_dna_sequence(const char *seq, int len)
{
    int i;
    
    for (i = 0; i < len; i++)
    {
        char c = toupper(seq[i]);
        if (!is_valid_nucleotide(c))
            return false;
    }
    
    return true;
}

/*
 * Generate complement of a DNA sequence
 */
PG_FUNCTION_INFO_V1(dna_complement);
Datum
dna_complement(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int len = dna_get_length(d);
    dna *result;
    char *seq = dna_get_str(d);
    int i;
    
    result = (dna *) palloc(VARHDRSZ + len);
    SET_VARSIZE(result, VARHDRSZ + len);
    
    for (i = 0; i < len; i++)
    {
        result->data[i] = complement_nucleotide(seq[i]);
    }
    
    pfree(seq);
    
    PG_RETURN_DNA_P(result);
}

/*
 * Generate reverse of a DNA sequence
 */
PG_FUNCTION_INFO_V1(dna_reverse);
Datum
dna_reverse(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int len = dna_get_length(d);
    dna *result;
    char *seq = dna_get_str(d);
    int i;
    
    result = (dna *) palloc(VARHDRSZ + len);
    SET_VARSIZE(result, VARHDRSZ + len);
    
    for (i = 0; i < len; i++)
    {
        result->data[i] = seq[len - 1 - i];
    }
    
    pfree(seq);
    
    PG_RETURN_DNA_P(result);
}

/*
 * Generate reverse complement of a DNA sequence
 */
PG_FUNCTION_INFO_V1(dna_reverse_complement);
Datum
dna_reverse_complement(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int len = dna_get_length(d);
    dna *result;
    char *seq = dna_get_str(d);
    int i;
    
    result = (dna *) palloc(VARHDRSZ + len);
    SET_VARSIZE(result, VARHDRSZ + len);
    
    for (i = 0; i < len; i++)
    {
        result->data[i] = complement_nucleotide(seq[len - 1 - i]);
    }
    
    pfree(seq);
    
    PG_RETURN_DNA_P(result);
}

/*
 * Generate k-mers from a DNA sequence
 */
PG_FUNCTION_INFO_V1(generate_kmers);
Datum
generate_kmers(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int32 k = PG_GETARG_INT32(1);
    int seq_len = dna_get_length(d);
    char *seq = dna_get_str(d);
    ArrayType *result;
    Datum *elems;
    int num_kmers;
    int i;
    Oid kmer_type_oid;
    
    if (k <= 0 || k > seq_len)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k must be between 1 and sequence length")));
    
    /* Get the OID of the kmer type */
    kmer_type_oid = TypenameGetTypid("kmer");
    
    num_kmers = seq_len - k + 1;
    elems = (Datum *) palloc(num_kmers * sizeof(Datum));
    
    for (i = 0; i < num_kmers; i++)
    {
        kmer *km = (kmer *) palloc(VARHDRSZ + sizeof(int32) + k);
        SET_VARSIZE(km, VARHDRSZ + sizeof(int32) + k);
        km->k = k;
        memcpy(km->data, seq + i, k);
        
        elems[i] = PointerGetDatum(km);
    }
    
    result = construct_array(elems, num_kmers, kmer_type_oid,
                           -1, false, 'i');
    
    pfree(seq);
    pfree(elems);
    
    PG_RETURN_ARRAYTYPE_P(result);
}

/*
 * IUPAC nucleotide utility functions
 */

/* Complement mapping table */
const char complement_map[128] = {
    ['A'] = 'T', ['T'] = 'A', ['C'] = 'G', ['G'] = 'C',
    ['a'] = 't', ['t'] = 'a', ['c'] = 'g', ['g'] = 'c',
    ['R'] = 'Y', ['Y'] = 'R', ['S'] = 'S', ['W'] = 'W',
    ['K'] = 'M', ['M'] = 'K', ['B'] = 'V', ['V'] = 'B',
    ['D'] = 'H', ['H'] = 'D', ['N'] = 'N', ['-'] = '-'
};

/*
 * Check if character is a valid nucleotide
 */
bool
is_valid_nucleotide(char c)
{
    c = toupper(c);
    return (c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
            c == 'R' || c == 'Y' || c == 'S' || c == 'W' ||
            c == 'K' || c == 'M' || c == 'B' || c == 'D' ||
            c == 'H' || c == 'V' || c == 'N' || c == '-');
}

/*
 * Check if character is an ambiguous nucleotide
 */
bool
is_ambiguous_nucleotide(char c)
{
    c = toupper(c);
    return (c == 'R' || c == 'Y' || c == 'S' || c == 'W' ||
            c == 'K' || c == 'M' || c == 'B' || c == 'D' ||
            c == 'H' || c == 'V' || c == 'N');
}

/*
 * Get complement of a nucleotide
 */
char
complement_nucleotide(char c)
{
    if (c >= 0 && c < 128 && complement_map[(int)c] != 0)
        return complement_map[(int)c];
    
    return c; /* Return original if no complement defined */
}

/*
 * Convert nucleotide to integer representation
 */
int
nucleotide_to_int(char c)
{
    switch (toupper(c))
    {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1; /* Invalid or ambiguous */
    }
}

/*
 * Convert integer to nucleotide
 */
char
int_to_nucleotide(int i)
{
    switch (i)
    {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

/*
 * Count occurrences of a nucleotide in a DNA sequence
 */
PG_FUNCTION_INFO_V1(dna_count);
Datum
dna_count(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    text *nucl_text = PG_GETARG_TEXT_PP(1);
    char nucl = *VARDATA_ANY(nucl_text);
    int len = dna_get_length(d);
    char *seq = dna_get_str(d);
    int count = 0;
    int i;
    
    for (i = 0; i < len; i++)
    {
        if (toupper(seq[i]) == toupper(nucl))
            count++;
    }
    
    pfree(seq);
    
    PG_RETURN_INT32(count);
}

/*
 * Approximate count of nucleotides in a DNA sequence (GC content)
 */
PG_FUNCTION_INFO_V1(dna_count_approx);
Datum
dna_count_approx(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int len = dna_get_length(d);
    char *seq = dna_get_str(d);
    int count = 0;
    int i;
    
    for (i = 0; i < len; i++)
    {
        char c = toupper(seq[i]);
        if (c == 'G' || c == 'C')
            count++;
    }
    
    pfree(seq);
    
    PG_RETURN_INT32(count);
}

/*
 * Calculate GC content of a DNA sequence
 */
PG_FUNCTION_INFO_V1(dna_gc_content);
Datum
dna_gc_content(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int len = dna_get_length(d);
    char *seq = dna_get_str(d);
    int gc_count = 0;
    int i;
    
    for (i = 0; i < len; i++)
    {
        char c = toupper(seq[i]);
        if (c == 'G' || c == 'C')
            gc_count++;
    }
    
    pfree(seq);
    
    if (len == 0)
        PG_RETURN_FLOAT8(0.0);
    
    PG_RETURN_FLOAT8((double)gc_count / len);
}

/*
 * Convert DNA sequence to string (for output)
 */
PG_FUNCTION_INFO_V1(dna_to_string);
Datum
dna_to_string(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    char *str = dna_get_str(d);
    
    PG_RETURN_CSTRING(str);
}

/*
 * Convert string to DNA sequence (from input)
 */
PG_FUNCTION_INFO_V1(string_to_dna);
Datum
string_to_dna(PG_FUNCTION_ARGS)
{
    text *t = PG_GETARG_TEXT_P(0);
    char *str = VARDATA(t);
    int len = VARSIZE_ANY_EXHDR(t);
    dna *result;
    int i;
    
    result = (dna *) palloc(VARHDRSZ + len);
    SET_VARSIZE(result, VARHDRSZ + len);
    
    for (i = 0; i < len; i++)
    {
        result->data[i] = toupper(str[i]);
    }
    
    PG_RETURN_DNA_P(result);
}