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

SELECT 'Extension successfully installed and working!' as status;