#ifndef IUPAC_HEADER_H
#define IUPAC_HEADER_H

/*
 * IUPAC nucleotide codes for DNA sequences
 * International Union of Pure and Applied Chemistry standard
 */

/* Standard nucleotides */
#define IUPAC_A 'A'  /* Adenine */
#define IUPAC_C 'C'  /* Cytosine */
#define IUPAC_G 'G'  /* Guanine */
#define IUPAC_T 'T'  /* Thymine */

/* Ambiguous nucleotides */
#define IUPAC_R 'R'  /* A or G (puRine) */
#define IUPAC_Y 'Y'  /* C or T (pYrimidine) */
#define IUPAC_S 'S'  /* G or C (Strong) */
#define IUPAC_W 'W'  /* A or T (Weak) */
#define IUPAC_K 'K'  /* G or T (Keto) */
#define IUPAC_M 'M'  /* A or C (aMino) */
#define IUPAC_B 'B'  /* C or G or T (not A) */
#define IUPAC_D 'D'  /* A or G or T (not C) */
#define IUPAC_H_CODE 'H'  /* A or C or T (not G) */
#define IUPAC_V 'V'  /* A or C or G (not T) */
#define IUPAC_N 'N'  /* Any nucleotide */
#define IUPAC_GAP '-' /* Gap character */

/* Function declarations */
bool is_valid_nucleotide(char c);
bool is_ambiguous_nucleotide(char c);
char complement_nucleotide(char c);
int nucleotide_to_int(char c);
char int_to_nucleotide(int i);

/* Nucleotide validation macros */
#define IS_STANDARD_DNA(c) ((c) == 'A' || (c) == 'C' || (c) == 'G' || (c) == 'T')
#define IS_VALID_IUPAC(c) (is_valid_nucleotide(c))

/* Complement mapping */
extern const char complement_map[128];

#endif /* IUPAC_H */