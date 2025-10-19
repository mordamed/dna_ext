-- DNA Extension Demo
-- This script demonstrates the basic functionality of the DNA extension

-- Create the extension
CREATE EXTENSION IF NOT EXISTS dna_ext;

-- Create a table for DNA sequences
CREATE TABLE dna_sequences (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255),
    sequence dna,
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Insert some sample DNA sequences
INSERT INTO dna_sequences (name, sequence, description) VALUES 
    ('Human Beta-globin fragment', 'ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACT', 'Part of human beta-globin gene'),
    ('E.coli promoter', 'TTGACAATTAATCATCGGCTCGTATAATGTGTGG', 'Bacterial promoter sequence'),
    ('Ribosome binding site', 'AGGAGGTT', 'Shine-Dalgarno sequence'),
    ('Start codon', 'ATG', 'Translation initiation codon'),
    ('Stop codon', 'TAA', 'Translation termination codon');

-- Show all sequences
SELECT 'All DNA sequences:' as demo;
SELECT id, name, sequence, description FROM dna_sequences ORDER BY id;

-- Demonstrate basic DNA operations
SELECT 'DNA Operations Demo:' as demo;

-- Length of sequences
SELECT name, sequence, length(sequence::text) as length 
FROM dna_sequences 
ORDER BY length DESC;

-- Find sequences containing ATG (start codon)
SELECT 'Sequences containing ATG:' as demo;
SELECT name, sequence 
FROM dna_sequences 
WHERE sequence::text LIKE '%ATG%';

-- Simple pattern matching
SELECT 'Pattern matching demo:' as demo;
SELECT name, sequence 
FROM dna_sequences 
WHERE sequence::text ~ '[AT]{3,}';  -- 3 or more A's or T's in a row

-- Test DNA utility functions
SELECT 'DNA Utility Functions:' as demo;

-- Test dna_length
SELECT 'Length of sequences:' as demo;
SELECT name, sequence, dna_length(sequence) as length 
FROM dna_sequences 
ORDER BY dna_length(sequence) DESC;

-- Test dna_complement
SELECT 'Complement of start codon ATG:' as demo;
SELECT dna_complement('ATG'::dna) as complement;

-- Test dna_reverse
SELECT 'Reverse of ATG:' as demo;
SELECT dna_reverse('ATG'::dna) as reversed;

-- Test dna_reverse_complement
SELECT 'Reverse complement of sequences:' as demo;
SELECT name, sequence, dna_reverse_complement(sequence) as reverse_complement 
FROM dna_sequences 
WHERE id <= 3;

-- Test dna_gc_content
SELECT 'GC content of sequences:' as demo;
SELECT name, sequence, ROUND(dna_gc_content(sequence)::numeric, 4) as gc_content 
FROM dna_sequences 
ORDER BY dna_gc_content(sequence) DESC;

-- Test dna_count_nucleotide
SELECT 'Count of A nucleotides:' as demo;
SELECT name, sequence, dna_count_nucleotide(sequence, 'A') as count_a 
FROM dna_sequences 
ORDER BY dna_count_nucleotide(sequence, 'A') DESC;

-- Test dna_is_palindrome
SELECT 'Check palindromic sequences:' as demo;
SELECT name, sequence, dna_is_palindrome(sequence) as is_palindrome 
FROM dna_sequences;

-- Test dna_find_subsequence
SELECT 'Find ATG in sequences:' as demo;
SELECT name, sequence, dna_find_subsequence(sequence, 'ATG'::dna) as atg_position 
FROM dna_sequences 
WHERE dna_find_subsequence(sequence, 'ATG'::dna) >= 0;

-- Test generate_kmers
SELECT 'Generate 3-mers from start codon:' as demo;
SELECT generate_kmers('ATGCAT'::dna, 3) as kmers;

-- Test dna_sliding_gc
SELECT 'Sliding GC content (window=3):' as demo;
SELECT name, dna_sliding_gc(sequence, 3) as sliding_gc 
FROM dna_sequences 
WHERE dna_length(sequence) >= 3 
LIMIT 3;

-- Test DNA operators
SELECT 'DNA Operators:' as demo;

-- Test contains operator @>
SELECT 'Sequences containing ATG:' as demo;
SELECT name, sequence 
FROM dna_sequences 
WHERE sequence @> 'ATG'::dna;

-- Test equality
SELECT 'Test equality:' as demo;
SELECT ('ATG'::dna = 'ATG'::dna) as equal, 
       ('ATG'::dna = 'TAA'::dna) as not_equal;

-- Test comparison operators
SELECT 'Comparison operators:' as demo;
SELECT name, sequence 
FROM dna_sequences 
ORDER BY sequence 
LIMIT 3;

-- Test similarity operator
SELECT 'Similarity between sequences:' as demo;
SELECT a.name as seq1, b.name as seq2, 
       ROUND((a.sequence ^@ b.sequence)::numeric, 4) as similarity 
FROM dna_sequences a, dna_sequences b 
WHERE a.id = 4 AND b.id = 5;

-- Test type conversion functions
SELECT 'Type conversion:' as demo;
SELECT dna_to_string('ATGC'::dna) as to_string,
       string_to_dna('GCTA') as from_string;

SELECT 'Extension successfully installed and working!' as status;