# DNA Extension for PostgreSQL

A PostgreSQL extension for handling DNA sequence data with specialized data types and operations.

## Repository Structure

```
dna_ext/
├── dna_ext.control              # Extension control file
├── makefile                     # Build configuration
├── README.md                    # This file
├── sql/
│   └── dna_ext--1.0.sql        # SQL definitions for types, functions, and operators
└── src/
    ├── dna.h                   # Main header file with type definitions
    ├── iupac.h                 # IUPAC nucleotide codes and utilities
    ├── dna_utils.c             # Core DNA utility functions
    ├── type_dna.c              # DNA type input/output functions
    ├── type_kmer.c             # K-mer type input/output functions
    ├── type_qkmer.c            # Quality k-mer type functions
    ├── funcs.c                 # Extended DNA analysis functions
    ├── ops.c                   # Comparison and containment operators
    ├── hash_ops.c              # Hash support for indexing
    ├── btree_ops.c             # B-tree support for indexing
    └── spgist_kmer.c           # SP-GiST trie implementation for k-mers
```

## Data Types

### DNA
- Stores variable-length DNA sequences
- Validates IUPAC nucleotide codes
- Supports complement and reverse complement operations

### K-mer
- Stores DNA subsequences of length k
- Optimized for pattern matching and indexing
- Supports trie-based indexing via SP-GiST

### QKmer (Quality K-mer)
- Stores k-mers with associated quality scores
- Useful for sequencing data analysis
- Supports quality filtering operations

## Features

### Core Functions
- `dna_length()` - Get sequence length
- `dna_complement()` - Generate complement sequence
- `dna_reverse_complement()` - Generate reverse complement
- `generate_kmers()` - Extract all k-mers from a sequence
- `dna_gc_content()` - Calculate GC content percentage

### Analysis Functions
- `dna_count_nucleotide()` - Count specific nucleotides
- `dna_find_subsequence()` - Find subsequence positions
- `dna_is_palindrome()` - Check for palindromic sequences
- `dna_translate()` - Translate DNA to amino acids
- `dna_sliding_gc()` - Sliding window GC analysis

### Quality Functions
- `qkmer_avg_quality()` - Average quality score
- `qkmer_min_quality()` - Minimum quality score
- `qkmer_filter_quality()` - Quality-based filtering

### Operators
- `=`, `<>`, `<`, `<=`, `>`, `>=` - Standard comparisons
- `@>` - Contains (sequence contains subsequence)
- `<@` - Contained by (subsequence contained in sequence)
- `&&` - Overlap (sequences share common subsequences)
- `^@` - Similarity (similarity score between sequences)

### Indexing Support
- **B-tree**: Standard ordering and range queries
- **Hash**: Equality comparisons and hash joins
- **SP-GiST**: Trie-based indexing for efficient k-mer searches

## Building

```bash
make
sudo make install
```

## Usage

```sql
-- Create extension
CREATE EXTENSION dna_ext;

-- Create table with DNA data
CREATE TABLE sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA
);

-- Insert DNA sequences
INSERT INTO sequences (name, sequence) VALUES 
    ('seq1', 'ATCGATCGATCG'),
    ('seq2', 'GCTAGCTAGCTA');

-- Query examples
SELECT name, dna_length(sequence) FROM sequences;
SELECT name, dna_gc_content(sequence) FROM sequences;
SELECT generate_kmers(sequence, 3) FROM sequences WHERE name = 'seq1';

-- Create indexes
CREATE INDEX idx_sequences_dna ON sequences USING btree (sequence);
CREATE INDEX idx_sequences_hash ON sequences USING hash (sequence);
```

## IUPAC Support

The extension supports all IUPAC nucleotide codes:
- Standard: A, C, G, T
- Ambiguous: R, Y, S, W, K, M, B, D, H, V, N
- Gap: -

## License

This extension is provided as-is for educational and research purposes.


## lancer demo

docker ps
docker run --name dna_demo_db -p 5432:5432 -d dna_ext_demo #lancer le conteneur s'il est dans la liste
docker exec -it dna_demo_db psql -U postgres -d dna_demo #se connecter 


src/
├── module.c          → Point d'entrée de l'extension
├── type_dna.c        → Type DNA (séquences ADN)
├── type_kmer.c       → Type K-mer (sous-séquences)
├── type_qkmer.c      → Type Q-Kmer (k-mers avec qualité)
├── dna_utils.c       → Fonctions utilitaires pour ADN
├── funcs.c           → Fonctions d'analyse avancées
├── ops.c             → Opérateurs de comparaison
├── btree_ops.c       → Support d'index B-tree
├── hash_ops.c        → Support d'index Hash
└── spgist_kmer.c     → Support d'index SP-GiST (trie)
| --- HEADERS ---
├── dna.h             → Fichier d'en-tête (définitions communes)
├── iupac.h           → Codes IUPAC pour nucléotides