# Inferring areas of introgressed Neanderthal ancestry using a Hidden Markov Model
___
## Questions
1. Not sure if path in shell script will work in OSCAR, was copied over from CS project
2. Could it be more effective if I stochastically modified the chances in the naive HMM a little first?

## Easy Testing command
`python3 extract_obs.py ../cs282_sim_data/rep_id_1_intro_pos.csv.gz ../cs282_sim_data/rep_id_1_polarized_geno_mat.csv.gz
`
## Input
###ancestries
a list of 4 haplotype sequences (ancestries)
- each haplotype sequence will be of the same length
- each element of each ancestry will display a 0 or 1, indicating an ancestral or derived variable site, respectively
- similar indices between ancestries indicate the same locus
- all variable sites are biallelic
- the order of the ancestries will be initialized as follows:
    - AFR1, one african ancestry
    - AFR2, another african ancestry
    - TEST, the ancestry that we are investigating for introgression
    - NEAN, the neanderthal ancestry that may have introgressed
- example:
- [
- (AFR1)[1, 1, 0, 0, 0, 0],
- (AFR2)[1, 0, 0, 0, 0, 1],
- (TEST)[0, 1, 1, 1, 1, 0],
- (NEAN)[0, 0, 0, 1, 1, 1]
- ]

###loci
a list of the genome positions of each SNP
- each index corresponds to a SNP position in an ancestry
- values are measured in kb
- example: [10(kb), 15(kb), 40(kb), 43(kb), 45(kb)]
___
## Steps
1. Label all sites in ancestries
   * PERHAPS NOT NECESSARY - JUST KEEP A COUNTER WHEN ITERATING THROUGH THE BINS
   * "C" for "Consistent with introgression"
   * "N" for "Not consistent with introgression"


2. Bin each phased halpoid genome based on proximity
   * bin_size = 0.00005 Morgans (non-overlapping windows)
   * Chosen to be a tenth of the typical genetic span of archaic segments
   * 0.0005 = 1/2000 Morgans, the length expected from admixture 2,000 generations old
   1. Classifying bins
      1. a bin is considered informative about ancestry if
         1. the number of consistent sites (nC) > 0 AND the number of non-consistent sites (nN) = 0 ==> bin is C
            1. (this supports the test haplotype being in a clade with archaic haplotypes)
         2. OR nC=0, nN>0 ==> bin is N
            1. (contradicting such a history)
      2. a bin is considered uninformative about ancestry if
         1. nC=0, nN=0 (data cleaning will ensure this never happens)
         2. nC > 0, nN > 0 => N
         3. Uninformative sites are classified as Non-consistent
   2. Creating Observed Sequence
      1. Each labeled bin is appended to list O, which will be the observed sequence in the Hidden Markov Model
      2. example: O = [N, C, N]
   3. Storing bin information
      1. To understand which loci get placed in which bin, whenever a bin is added to the observed sequence, we store its infomration
      2. Create the dictionary Bin_dict
         1. key t stores information about O[t]
         2. example: Bin_dict[0] ==> [10, 15, 40, 45]
         3. (Optional) - this data structure can also store other information
            1. Variant site labels within bins [N, N, N, C]
            2. Ancestral variant patterns [1100, 1010, 0010, 0011]
            
3. Create the Hidden Markov Model!
   1. Assumptions from Prufer 2014 supplemental materials
      1. Ancestry switch rate is 0.0005 Morgans
      2. Prior probability for archaic ancestry at any locus is 0.01
         1. Corresponds to Neanderthal admixture proportion of 1%
         2. Conservatively at the low end of the range of what has been previously reported
      3. Probability of archaic ancestry conditional on all SNPs in the window being of state "C" is .99
         1. For the sake of reducing our parameter space it will also be the probability of modern human ancestry conditional on all SNPs in a window being of state "N"
         2. Allows for some probability of a test haplotype not being introgressed even if locally it is a clade with an archaic haplotype due to incomplete lineage sorting
         3. CONFIRM - This setting means that a locus is only called as confidently introgressed if at least two windows in a row are "C", thus discriminating against haplotypes that owe their origin to incomplete lineage sorting and that we do not wish to call as Neandertal introgressed.
   2. 