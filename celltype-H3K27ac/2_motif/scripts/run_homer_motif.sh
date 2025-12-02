# HOMER Motif Analysis
# Note: The same command works for both intersect and non-intersect analysis.
#       Just change the input BED file accordingly.

# ============================================================
# Installation (using conda)
# ============================================================
# conda create -n homer -c bioconda homer
# conda activate homer
# which findMotifsGenome.pl

# ============================================================
# Usage
# ============================================================
# findMotifsGenome.pl <input.bed> <genome> <output_dir> [options]
#   -size 200      : Region size (±100bp around peak center)
#   -len 8,10,12   : Motif lengths to search
#   -S 10          : Number of motifs to find

# ============================================================
# Run motif analysis for each cell type
# ============================================================

# Neg (Microglia)
findMotifsGenome.pl \
  ./merged_final/Neg_clean.bed \
  rn7 ./motif/Neg -size 200 -len 8,10,12 -S 10

# NeuN (Neurons)
findMotifsGenome.pl \
  ./merged_final/NeuN_clean.bed \
  rn7 ./motif/NeuN -size 200 -len 8,10,12 -S 10

# Nurr (Dopaminergic)
findMotifsGenome.pl \
  ./merged_final/Nurr_clean.bed \
  rn7 ./motif/Nurr -size 200 -len 8,10,12 -S 10

# Olig (Oligodendrocyte)
findMotifsGenome.pl \
  ./merged_final/Olig_clean.bed \
  rn7 ./motif/Olig -size 200 -len 8,10,12 -S 10
