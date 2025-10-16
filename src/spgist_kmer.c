#include "dna.h"

/*
 * SP-GiST (Space Partitioned Generalized Search Tree) implementation for K-mers
 * This implements a trie-based index structure for efficient k-mer searches
 */

/* SP-GiST node for k-mer trie */
typedef struct KmerTrieNode
{
    int level;              /* Current level in the trie (0 to k-1) */
    char nucleotide;        /* Nucleotide at this level */
    bool is_leaf;          /* True if this is a leaf node */
} KmerTrieNode;

/*
 * SP-GiST config function for k-mer indexing
 */
PG_FUNCTION_INFO_V1(spgist_kmer_config);
Datum
spgist_kmer_config(PG_FUNCTION_ARGS)
{
    spgConfigOut *cfg = (spgConfigOut *) PG_GETARG_POINTER(1);
    
    cfg->prefixType = VOIDOID;      /* No prefix data */
    cfg->labelType = CHAROID;       /* Labels are nucleotide characters */
    cfg->leafType = InvalidOid;     /* Use the k-mer type itself */
    cfg->canReturnData = true;      /* Can return the indexed data */
    cfg->longValuesOK = true;       /* Handle variable-length k-mers */
    
    PG_RETURN_VOID();
}

/*
 * SP-GiST choose function for k-mer indexing
 * Decides which branch to follow when inserting/searching
 */
PG_FUNCTION_INFO_V1(spgist_kmer_choose);
Datum
spgist_kmer_choose(PG_FUNCTION_ARGS)
{
    spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);
    kmer *k = DatumGetKmerP(in->datum);
    KmerTrieNode *node = (KmerTrieNode *) DatumGetPointer(in->prefixDatum);
    int level = node ? node->level : 0;
    
    /* Check if we've reached the end of the k-mer */
    if (level >= k->k)
    {
        out->resultType = spgMatchNode;
        PG_RETURN_VOID();
    }
    
    /* Get the nucleotide at the current level */
    char nucleotide = k->data[level];
    
    /* Search for matching child node */
    for (int i = 0; i < in->nNodes; i++)
    {
        char label = DatumGetChar(in->nodeLabels[i]);
        if (label == nucleotide)
        {
            out->resultType = spgMatchNode;
            out->result.matchNode.nodeN = i;
            PG_RETURN_VOID();
        }
    }
    
    /* No matching child found, need to add new node */
    out->resultType = spgAddNode;
    out->result.addNode.nodeLabel = CharGetDatum(nucleotide);
    
    /* Create new node data */
    KmerTrieNode *new_node = (KmerTrieNode *) palloc(sizeof(KmerTrieNode));
    new_node->level = level + 1;
    new_node->nucleotide = nucleotide;
    new_node->is_leaf = (level + 1 >= k->k);
    
    out->result.addNode.nodeN = 0; /* Will be assigned by the system */
    
    PG_RETURN_VOID();
}

/*
 * SP-GiST picksplit function for k-mer indexing
 * Decides how to split an overfull node
 */
PG_FUNCTION_INFO_V1(spgist_kmer_picksplit);
Datum
spgist_kmer_picksplit(PG_FUNCTION_ARGS)
{
    spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
    spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);
    int level = in->level;
    int i;
    
    /* Count occurrences of each nucleotide at the current level */
    int counts[4] = {0}; /* A, C, G, T */
    char nucleotides[4] = {'A', 'C', 'G', 'T'};
    
    for (i = 0; i < in->nTuples; i++)
    {
        kmer *k = DatumGetKmerP(in->datums[i]);
        if (level < k->k)
        {
            char nucleotide = k->data[level];
            switch (nucleotide)
            {
                case 'A': counts[0]++; break;
                case 'C': counts[1]++; break;
                case 'G': counts[2]++; break;
                case 'T': counts[3]++; break;
            }
        }
    }
    
    /* Count how many different nucleotides we have */
    int num_nodes = 0;
    for (i = 0; i < 4; i++)
    {
        if (counts[i] > 0)
            num_nodes++;
    }
    
    /* Set up the split */
    out->nNodes = num_nodes;
    out->nodeLabels = (Datum *) palloc(sizeof(Datum) * num_nodes);
    out->mapTuplesToNodes = (int *) palloc(sizeof(int) * in->nTuples);
    out->leafTupleDatums = (Datum *) palloc(sizeof(Datum) * in->nTuples);
    
    /* Create node labels */
    int node_index = 0;
    for (i = 0; i < 4; i++)
    {
        if (counts[i] > 0)
        {
            out->nodeLabels[node_index] = CharGetDatum(nucleotides[i]);
            node_index++;
        }
    }
    
    /* Map tuples to nodes */
    for (i = 0; i < in->nTuples; i++)
    {
        kmer *k = DatumGetKmerP(in->datums[i]);
        char nucleotide = k->data[level];
        
        /* Find which node this tuple belongs to */
        node_index = 0;
        for (int j = 0; j < 4; j++)
        {
            if (counts[j] > 0)
            {
                if (nucleotides[j] == nucleotide)
                {
                    out->mapTuplesToNodes[i] = node_index;
                    break;
                }
                node_index++;
            }
        }
        
        /* Store the tuple data */
        out->leafTupleDatums[i] = in->datums[i];
    }
    
    /* Set prefix for child nodes */
    KmerTrieNode *new_node = (KmerTrieNode *) palloc(sizeof(KmerTrieNode));
    new_node->level = level + 1;
    new_node->nucleotide = '\0'; /* Not specific to any nucleotide */
    new_node->is_leaf = false;
    
    out->prefixDatum = PointerGetDatum(new_node);
    
    PG_RETURN_VOID();
}

/*
 * SP-GiST inner consistent function for k-mer indexing
 * Determines which child nodes to visit during a search
 */
PG_FUNCTION_INFO_V1(spgist_kmer_inner_consistent);
Datum
spgist_kmer_inner_consistent(PG_FUNCTION_ARGS)
{
    spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
    spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);
    int level = in->level;
    
    /* For now, implement a simple exact match search */
    /* In a full implementation, you would handle various search strategies */
    
    out->nNodes = 0;
    out->nodeNumbers = (int *) palloc(sizeof(int) * in->nNodes);
    out->levelAdds = (int *) palloc(sizeof(int) * in->nNodes);
    
    /* Check each query condition */
    for (int i = 0; i < in->nkeys; i++)
    {
        ScanKey key = &in->scankeys[i];
        
        if (key->sk_strategy == BTEqualStrategyNumber)
        {
            kmer *query_kmer = DatumGetKmerP(key->sk_argument);
            
            if (level < query_kmer->k)
            {
                char target_nucleotide = query_kmer->data[level];
                
                /* Find matching child nodes */
                for (int j = 0; j < in->nNodes; j++)
                {
                    char node_nucleotide = DatumGetChar(in->nodeLabels[j]);
                    if (node_nucleotide == target_nucleotide)
                    {
                        out->nodeNumbers[out->nNodes] = j;
                        out->levelAdds[out->nNodes] = 1;
                        out->nNodes++;
                    }
                }
            }
        }
    }
    
    PG_RETURN_VOID();
}

/*
 * SP-GiST leaf consistent function for k-mer indexing
 * Tests whether a leaf tuple matches the search conditions
 */
PG_FUNCTION_INFO_V1(spgist_kmer_leaf_consistent);
Datum
spgist_kmer_leaf_consistent(PG_FUNCTION_ARGS)
{
    spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
    spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);
    kmer *leaf_kmer = DatumGetKmerP(in->leafDatum);
    bool match = false;
    
    /* Check each query condition */
    for (int i = 0; i < in->nkeys; i++)
    {
        ScanKey key = &in->scankeys[i];
        
        switch (key->sk_strategy)
        {
            case BTEqualStrategyNumber:
            {
                kmer *query_kmer = DatumGetKmerP(key->sk_argument);
                match = (kmer_compare_internal(leaf_kmer, query_kmer) == 0);
                break;
            }
            case BTLessStrategyNumber:
            {
                kmer *query_kmer = DatumGetKmerP(key->sk_argument);
                match = (kmer_compare_internal(leaf_kmer, query_kmer) < 0);
                break;
            }
            case BTLessEqualStrategyNumber:
            {
                kmer *query_kmer = DatumGetKmerP(key->sk_argument);
                match = (kmer_compare_internal(leaf_kmer, query_kmer) <= 0);
                break;
            }
            case BTGreaterStrategyNumber:
            {
                kmer *query_kmer = DatumGetKmerP(key->sk_argument);
                match = (kmer_compare_internal(leaf_kmer, query_kmer) > 0);
                break;
            }
            case BTGreaterEqualStrategyNumber:
            {
                kmer *query_kmer = DatumGetKmerP(key->sk_argument);
                match = (kmer_compare_internal(leaf_kmer, query_kmer) >= 0);
                break;
            }
            default:
                match = false;
                break;
        }
        
        if (!match)
            break; /* AND logic: all conditions must be satisfied */
    }
    
    out->recheck = false; /* No need to recheck */
    out->leafValue = in->leafDatum;
    
    PG_RETURN_BOOL(match);
}