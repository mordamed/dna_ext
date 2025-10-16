#include "dna.h"

/*
 * Additional DNA/K-mer utility functions
 * Functions for sequence analysis and manipulation
 */

/*
 * Count occurrences of a specific nucleotide
 */
PG_FUNCTION_INFO_V1(dna_count_nucleotide);
Datum
dna_count_nucleotide(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    char target = PG_GETARG_CHAR(1);
    int len = dna_get_length(d);
    char *seq = dna_get_str(d);
    int count = 0;
    int i;
    
    target = toupper(target);
    
    for (i = 0; i < len; i++)
    {
        if (toupper(seq[i]) == target)
            count++;
    }
    
    pfree(seq);
    
    PG_RETURN_INT32(count);
}

/*
 * Find the first occurrence of a subsequence
 */
PG_FUNCTION_INFO_V1(dna_find_subsequence);
Datum
dna_find_subsequence(PG_FUNCTION_ARGS)
{
    dna *haystack = PG_GETARG_DNA_P(0);
    dna *needle = PG_GETARG_DNA_P(1);
    char *haystack_str = dna_get_str(haystack);
    char *needle_str = dna_get_str(needle);
    char *found;
    int position = -1;
    
    found = strstr(haystack_str, needle_str);
    if (found)
        position = found - haystack_str;
    
    pfree(haystack_str);
    pfree(needle_str);
    
    PG_RETURN_INT32(position);
}

/*
 * Check if sequence is palindromic
 */
PG_FUNCTION_INFO_V1(dna_is_palindrome);
Datum
dna_is_palindrome(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int len = dna_get_length(d);
    char *seq = dna_get_str(d);
    bool is_palindrome = true;
    int i;
    
    for (i = 0; i < len / 2; i++)
    {
        if (complement_nucleotide(seq[i]) != seq[len - 1 - i])
        {
            is_palindrome = false;
            break;
        }
    }
    
    pfree(seq);
    
    PG_RETURN_BOOL(is_palindrome);
}

/*
 * Translate DNA to amino acid sequence (single frame)
 */
PG_FUNCTION_INFO_V1(dna_translate);
Datum
dna_translate(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int frame = PG_GETARG_INT32(1); /* 0, 1, or 2 */
    int len = dna_get_length(d);
    char *seq = dna_get_str(d);
    text *result;
    char *aa_seq;
    int aa_len;
    int i, j;
    
    /* Genetic code table (simplified) */
    static const char *genetic_code[64] = {
        "K", "N", "K", "N", "T", "T", "T", "T", /* TTT, TTC, TTA, TTG, TCT, TCC, TCA, TCG */
        "R", "S", "R", "S", "I", "I", "I", "M", /* TAT, TAC, TAA, TAG, TGT, TGC, TGA, TGG */
        "Q", "H", "Q", "H", "P", "P", "P", "P", /* CAT, CAC, CAA, CAG, CCT, CCC, CCA, CCG */
        "R", "R", "R", "R", "L", "L", "L", "L", /* CGT, CGC, CGA, CGG, CTT, CTC, CTA, CTG */
        "E", "D", "E", "D", "A", "A", "A", "A", /* GAT, GAC, GAA, GAG, GCT, GCC, GCA, GCG */
        "G", "G", "G", "G", "V", "V", "V", "V"  /* GGT, GGC, GGA, GGG, GTT, GTC, GTA, GTG */
    };
    
    if (frame < 0 || frame > 2)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("frame must be 0, 1, or 2")));
    
    if (len < frame + 3)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("sequence too short for translation")));
    
    aa_len = (len - frame) / 3;
    aa_seq = palloc(aa_len + 1);
    
    for (i = frame, j = 0; i + 2 < len; i += 3, j++)
    {
        int codon_index = 0;
        char c1 = toupper(seq[i]);
        char c2 = toupper(seq[i + 1]);
        char c3 = toupper(seq[i + 2]);
        
        /* Simple codon to index mapping (incomplete) */
        if (c1 == 'T' && c2 == 'T' && c3 == 'T')
            codon_index = 0; /* TTT -> F */
        else
            codon_index = 0; /* Default to first amino acid */
        
        aa_seq[j] = 'X'; /* Placeholder - full implementation would use proper genetic code */
    }
    
    aa_seq[j] = '\0';
    
    result = cstring_to_text(aa_seq);
    
    pfree(seq);
    pfree(aa_seq);
    
    PG_RETURN_TEXT_P(result);
}

/*
 * Generate sliding window statistics
 */
PG_FUNCTION_INFO_V1(dna_sliding_gc);
Datum
dna_sliding_gc(PG_FUNCTION_ARGS)
{
    dna *d = PG_GETARG_DNA_P(0);
    int window_size = PG_GETARG_INT32(1);
    int len = dna_get_length(d);
    char *seq = dna_get_str(d);
    ArrayType *result;
    Datum *elems;
    int num_windows;
    int i;
    
    if (window_size <= 0 || window_size > len)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("window size must be between 1 and sequence length")));
    
    num_windows = len - window_size + 1;
    elems = (Datum *) palloc(num_windows * sizeof(Datum));
    
    for (i = 0; i < num_windows; i++)
    {
        int gc_count = 0;
        int j;
        
        for (j = i; j < i + window_size; j++)
        {
            char c = toupper(seq[j]);
            if (c == 'G' || c == 'C')
                gc_count++;
        }
        
        elems[i] = Float8GetDatum((double)gc_count / window_size);
    }
    
    result = construct_array(elems, num_windows, FLOAT8OID, 8, true, 'd');
    
    pfree(seq);
    pfree(elems);
    
    PG_RETURN_ARRAYTYPE_P(result);
}