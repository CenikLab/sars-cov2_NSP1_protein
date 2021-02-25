
# A NOTE ON NSP1 AND NSP2 SEQUENCES

For NSP1 and NSP2 genes, we used two sequences.

## VIRAL SEQUENCE

These are the original sequences existing in the virus genome.
We obtained these sequences from the virus via NSBI website : https://www.ncbi.nlm.nih.gov/nuccore/NC_045512

For NSP1, we took the nucleotides from position  265 to 805
For NSP2, we took the nucleotides from position  805 to 2719
We added the start codon "ATG" to the beginning of the NSP2 sequence.

For the above coordinates, we used 0 based indexing (the first nucleotide is at position 0).
The last nucleotides are excluded whereas the first nucleotide is included.

The reference files having viral sequences are indicated by "virus" in the file names.


## VECTOR SEQUENCES

The NSP1 and NSP2 RNA sequences, we made the cells express, have different nucleotide composition.
If the reference files don't have "virus" in their names as a substring, they contain these sequences.

As a sanity check, we verified that both RNA (technicxally DNA) sequences produce the same protein sequence.
More precisely, using https://web.expasy.org/translate/, we verified that both NSP1 sequences translate to

**MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRSGETLGVLVPHVGEIPVAYRKVLLRKNGNKGAGGHSYGADLKSFDLGDELGTDPYEDFQENWNTKHSSGVTRELMRELNGG**

Both NSP2 sequences translate to
**MAYTRYVDNNFCGPDGYPLECIKDLLARAGKASCTLSEQLDFIDTKRGVYCCREHEHEIAWYTERSEKSYELQTPFEIKLAKKFDTFNGECPNFVFPLNSIIKTIQPRVEKKKLDGFMGRIRSVYPVASPNECNQMCLSTLMKCDHCGETSWQTGDFVKATCEFCGTENLTKEGATTCGYLPQNAVVKIYCPACHNSEVGPEHSLAEYHNESGLKTILRKGGRTIAFGGCVFSYVGCHNKCAYWVPRASANIGCNHTGVVGEGSEGLNDNLLEILQKEKVNINIVGDFKLNEEIAIILASFSASTSAFVETVKGLDYKAFKQIVESCGNFKVTKGKAKKGAWNIGEQKSILSPLYAFASEAARVVRSIFSRTLETAQNSVRVLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYITGGVVQLTSQWLTNIFGTVYEKLKPVLDWLEEKFKEGVEFLRDGWEIVKFISTCACEIVGGQIVTCAKEIKESVQTFFKLVNKFLALCADSIIIGGAKLKALNLGETFVTHSKGLYRKCVKSREETGLLMPLKAPKEIIFLEGETLPTEVLTEEVVLKTGDLQPLEQPTSEAVEAPLVGTPVCINGLMLLEIKDTEKYCALAPNMMVTNNTFTLKGG**


