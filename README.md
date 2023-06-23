# softconsensus

Data sets and example code accompaning the manuscript

*Phylogenetic Information as Soft Constraints in RNA Secondary Structure Prediction*

Packages `numpy`, `biopython`, and `RNAlib` from ViennaRNA are required for `softconsensus.py`. The file provides a class object called `Alignment` which computes internally the pseudo-energies and proper soft constraints for given alignment and sequence to fold.

For the usage example, we invite users to look at the notebook `RNAsoftconsensus.ipynb`. Additional package `tabulate` is needed to disply the MCC table.
