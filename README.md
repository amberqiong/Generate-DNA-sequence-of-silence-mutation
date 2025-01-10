# Generate-silence-mutated-DNA-sequence

This process focuses on generating mutated DNA sequences for a given protein-coding DNA sequence while preserving the encoded amino acid sequence. The goal is to introduce all possible silent mutations—changes in the DNA sequence that do not alter the resulting protein—and ensure the mutated sequences have comparable %GC content to the original sequence.

Steps involved:
* Input: a protein coding DNA sequence starting from ATG
* Translation: DNA sequence is translated into amino acid sequence using codon
* Mutant DNA sequence generation: codon used in the original DNA sequence is excluded from replacements. A random synonymous codon is selected for each amino acid from available options. Repeated to generate multiple unique silent mutant sequences.
* GC content calculation: GC% is calculated to retain those with GC% within +-10% range of the original DNA sequence.

Two example sequences are included, as inx1 and inx3 from fly.
