original_inx1="ATGGCGGTCTTTGGCATGGTCTCGGCTGTGTCCGGCTTTATCAAGATACGCTATCTGCTGGACAAGGCGGTCATCGACAACATGGTTTTCCGCTGCCACTACAGGATCACGACGGCCATTCTGTTCACTTGTTGCATCATCGTCACCGCAAACAATCTGATCGGCGATCCCATAAGTTGCATCAACGACGGCGCCATTCCGATGCACGTAATCAACACCTTCTGCTGGATCACCTACACGTACACGATACCTGGCCAGCAGCATCGGCAAATCGGAACGGACGTGGCTGGTCCGGGATTGGGCAATGAGTACGGCCAGGAGAAGCGCTATCACAGCTACTATCAATGGGTGCCATTCGTGCTATTCTTCCAGGGGCTCATGTTCTACGTGCCCCACTGGGTTTGGAAGAACATGGAAGACGGCAAGATCCGCATGATCACCGATGGACTGCGCGGCATGGTCAGTGTGCCGGATGACTATCGGCGTGATCGCCAGGACCGGATCCTCAAGTATTTCGTGAACAGTTTGAACACCCACAACGGCTACTCGTTCGCATACTTTTTCTGCGAACTGCTTAACTTTATCAACGTGATTGTGAATATCTTTATGGTGGATAAGTTTCTGGGCGGTGCTTTCATGTCCTACGGTACGGATGTGCTCAAGTTCTCCAACATGGATCAGGACAAGCGCTTCGATCCGATGATCGAGATCTTTCCCCGGCTCACCAAGTGCACGTTCCACAAGTTCGGTCCCAGTGGATCGGTCCAGAAACACGACACTCTCTGTGTGCTGGCACTGAATATCCTCAACGAGAAGATCTACATCTTCCTGTGGTTCTGGTTTATCATCTTGGCCACCATCTCCGGCGTGGCTGTGCTCTATTCACTGGTGGTTATCATGATGCCCACCACTCGGGAGACCATCATCAAGCGTTCTTATCGCTCGGCACAACGGAAGGAAATCGCCGGACTGGTTAGGCGTTTGGAGATTGGTGACTTCCTGATCTTGCACTTCCTGAGCCAGAACCTCAGCACAAGATCGTACAGTGACATGCTGCAGCAGCTCTGCGGTCTGCTAGGAGCATCCCGAACACCATCCGCACCATCGACCCTTGAAATGAACCGCATCAGCCATCCGATTTATCCGCCCGTGGAGACCTTCGGCGGTGGCAAGGAGACGGAGACATGA"
original_inx3="ATGTATAAGTTGCTGGGTAGCCTGAAGAGCTACCTCAAGTGGCAGGACATCCAGACGGACAACGCCGTCTTCCGGCTGCACAACTCCTTTACCACGGTGCTCCTGCTAACCTGCAGCCTGATCATCACCGCCACCCAGTACGTGGGCCAGCCGATTAGCTGCATCGTCAATGGCGTACCGCCGCACGTGGTCAACACGTTCTGCTGGATCCACAGCACTTTCACCATGCCGGACGCTTTTCGCAGACAGGTTGGCCGAGAGGTGGCTCATCCCGGTGTGGCCAATGATTTTGGCGACGAG
GATGCCAAGA AGTACTACAC CTACTACCAG TGGGTGTGCT TCGTGCTTTT CTTCCAGGCC
ATGGCCTGTT ATACGCCCAA ATTCCTGTGG AATAAATTCG AGGGCGGACT GATGCGCATG
ATTGTGATGG GTCTGAATAT CACGATCTGC ACCCGCGAGG AGAAGGAGGC CAAACGCGAT
GCCCTGCTGG ACTATCTGAT CAAGCACGTG AAGCGCCACA AGCTGTACGC CATTCGGTAC
TGGGCCTGCG AATTTCTCTG CTGCATCAAC ATTATCGTGC AGATGTATCT GATGAATCGC
TTTTTCGATG GCGAGTTCCT CTCGTACGGT ACGAATATCA TGAAGCTTTC GGATGTGCCG
CAGGAGCAAA GGGTGGATCC CATGGTCTAT GTGTTCCCCC GGGTGACCAA GTGCACCTTC
CACAAGTATG GTCCCTCTGG TTCGCTGCAG AAGCACGACT CACTCTGCAT CCTGCCGCTG
AACATTGTGA ACGAGAAGAC GTACGTGTTC ATCTGGTTCT GGTTCTGGAT CCTGCTCGTC
CTGCTCATCG GACTGATAGT GTTCCGTGGC TGCATTATCT TTATGCCGAA ATTCCGACCC
CGCCTCCTGA ACGCCAGCAA TCGCATGATT CCGATGGAGA TCTGTCGCTC GCTGTCCCGC
AAACTGGACA TCGGTGACTG GTGGCTAATC TATATGCTGG GTCGCAATCT TGATCCGGTC
ATCTACAAGG ACGTGATGAG CGAGTTTGCC AAGCAGGTGG AGCCCTCCAA GCACGACCGT
GCCAAGTAG
"
original_inx3 <- gsub(" ", "", original_inx3)
original_inx3 <- gsub("\n", "", original_inx3)
nchar(original_inx3)

codon_table=list("TTT"="F","TTC"="F",
                 "TTA"="L","TTG"="L","CTT"="L","CTC"="L","CTA"="L","CTG"="L",
                 "ATT"="I","ATC"="I","ATA"="I",
                 "ATG"="M",
                 "GTT"="V","GTC"="V","GTA"="V","GTG"="V",
                 "TCT"="S","TCC"="S","TCA"="S","TCG"="S","AGT"="S","AGC"="S",
                 "CCT"="P","CCC"="P","CCA"="P","CCG"="P",
                 "ACT"="T","ACC"="T","ACA"="T","ACG"="T",
                 "GCT"="A","GCC"="A","GCA"="A","GCG"="A",
                 "TAT"="Y","TAC"="Y",
                 "TAA"="*","TAG"="*","TGA"="*",
                 "CAT"="H","CAC"="H",
                 "CAA"="Q","CAG"="Q",
                 "AAT"="N","AAC"="N",
                 "AAA"="K","AAG"="K",
                 "GAT"="D","GAC"="D",
                 "GAA"="E","GAG"="E",
                 "TGT"="C","TGC"="C",
                 "TGG"="W",
                 "CGT"="R","CGC"="R","CGA"="R","CGG"="R","AGA"="R","AGG"="R",
                 "GGT"="G","GGC"="G","GGA"="G","GGG"="G")
DNA_to_AA=function(DNAsequence){
  n=nchar(DNAsequence)
  AA=character(n %/% 3)
  
  for (i in seq_along(AA)){
    start=3*(i-1)+1
    end=start+2
    codon=substr(DNAsequence,start,end)
    print(i)
    print(codon)
    print(codon_table[[codon]])
    AA[i]=codon_table[[codon]]
    
  }
  #return(paste(AA,collapse = ""))
  return(AA)
}

inx1_AA=test
calculate_GC_content <- function(sequence) {
  num_GC <- sum(strsplit(sequence, "")[[1]] %in% c("G", "C"))
  gc_content <- num_GC / nchar(sequence)
  return(gc_content)
}

AA_to_DNA=function(AA_sequence,original_DNA,num_sequence=1000){
  mutant_sequences <- replicate(num_sequence,{
    mutant_sequence=""
    for (i in 1:length(AA_sequence)){
      #print(i)
      codons = names(codon_table)[codon_table == AA_sequence[i]]
      if (length(codons)==1) {
        selected_codon=codons
        #print(paste("the only avaiilable codon for",AA_sequence[i], "is", toString(selected_codon)))
      }else{
        available_codons=codons[!codons %in% substr(original_DNA,3*(i-1)+1,3*i)]
        #print(paste("Available codons for AA", i, ":", toString(available_codons)))
        selected_codon=sample(available_codons,1)
        #print(paste("Selected codon for AA", i, ":", selected_codon))
      }
      mutant_sequence=paste(mutant_sequence,selected_codon,sep = "")
    }
    return(mutant_sequence)
  })
  return(mutant_sequences)
}

mutant_inx1=AA_to_DNA(AA_sequence = test,original_DNA = original_inx1,num_sequence = 50000 )
mutant_GC_contents <- sapply(mutant_inx1,calculate_GC_content)

inx3_AA=DNA_to_AA(original_inx3)
mutant_inx3=AA_to_DNA(AA_sequence = inx3_AA,original_DNA = original_inx3,num_sequence = 50000 )
mutant_inx3_GC_contents <- sapply(mutant_inx3,calculate_GC_content)
original_inx3_GC_content<-calculate_GC_content(original_inx3)

GC_content_min=original_inx1_GC-0.1
GC_content_max=original_inx1_GC+0.1

valid_indices <- which(mutant_inx3_GC_contents >=0.412)

mutant_inx1[valid_indices]

writeLines(mutant_inx3[valid_indices], "inx3_mutant_sequences.txt")
