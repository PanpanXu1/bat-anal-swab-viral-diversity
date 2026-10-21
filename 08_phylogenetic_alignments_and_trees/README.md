# Phylogenetic Alignments and Trees

This archive contains the supplied sequence alignments and phylogenetic trees associated with Figures 4–5 and Supplementary Figures 1, 6–22 and 24. It contains **60 paired datasets**: 60 FASTA files and 60 tree files across 21 figure directories, totaling **91,642,323 bytes** (approximately 91.6 MB).

## Contents

| Figure directory | Dataset description | Alignment/tree pairs |
|---|---|---:|
| [Figure4/](Figure4/) | Viral-group phylogenies for Figure 4a | 19 |
| [Figure5/](Figure5/) | Coronavirus whole-genome and Merbecovirus/Sarbecovirus RBD datasets for Figure 5b–d | 3 |
| [Sfigure1/](Sfigure1/) | Bat-species cytb dataset for Supplementary Figure 1a | 1 |
| [Sfigure6/](Sfigure6/) | Cluster-specific datasets | 11 |
| [Sfigure7/](Sfigure7/) | Alphatetraviridae, Hepeviridae, Lispiviridae, Mymonaviridae, Paramyxoviridae, Phasmaviridae and Phenuiviridae | 7 |
| [Sfigure8/](Sfigure8/) | Rhabdoviridae | 1 |
| [Sfigure9/](Sfigure9/) | Amarillovirales and Rotavirus A | 2 |
| [Sfigure10/](Sfigure10/) | Merbecovirus | 1 |
| [Sfigure11/](Sfigure11/) | Sarbecovirus | 1 |
| [Sfigure12/](Sfigure12/) | Coronavirus whole genome | 1 |
| [Sfigure13/](Sfigure13/) | Coronavirus ORF1ab | 1 |
| [Sfigure14/](Sfigure14/) | Coronavirus ORF1ab | 1 |
| [Sfigure15/](Sfigure15/) | Coronavirus S | 1 |
| [Sfigure16/](Sfigure16/) | Coronavirus S | 1 |
| [Sfigure17/](Sfigure17/) | Coronavirus E | 1 |
| [Sfigure18/](Sfigure18/) | Coronavirus E | 1 |
| [Sfigure19/](Sfigure19/) | Coronavirus M | 1 |
| [Sfigure20/](Sfigure20/) | Coronavirus M | 1 |
| [Sfigure21/](Sfigure21/) | Coronavirus N | 1 |
| [Sfigure22/](Sfigure22/) | Coronavirus N | 1 |
| [Sfigure24/](Sfigure24/) | ORF1a, ORF1b and Region3 datasets | 3 |

The [file manifest](file_manifest.tsv) lists every alignment/tree pair, its relative paths and file sizes in bytes. Paths are relative to this directory.

## File Organization and Formats

Files are organized by figure directory and named by dataset. Files with the same stem form a supplied alignment/tree pair:

- `.fasta`: sequence alignment in FASTA format.
- `.treefile`: phylogenetic tree in Newick format.

For example, [Figure5b_Coronavirus_Whole_genome.fasta](Figure5/Figure5b/Figure5b_Coronavirus_Whole_genome.fasta) is paired with [Figure5b_Coronavirus_Whole_genome.treefile](Figure5/Figure5b/Figure5b_Coronavirus_Whole_genome.treefile). Figure 4a, Figure 5b–d and Supplementary Figure 1a retain their panel-level subdirectories. Sequence identifiers, sequence content and tree content are preserved as supplied.

## Shared Figure Datasets

The following figure directories contain byte-identical alignment/tree pairs. Each copy is retained under its figure directory to preserve the figure mapping:

- Figure 5b and Supplementary Figure 12: coronavirus whole genome.
- Supplementary Figures 13 and 14: coronavirus ORF1ab.
- Supplementary Figures 15 and 16: coronavirus S.
- Supplementary Figures 17 and 18: coronavirus E.
- Supplementary Figures 19 and 20: coronavirus M.
- Supplementary Figures 21 and 22: coronavirus N.

## Access and Use

Download the repository or the desired files directly from GitHub. FASTA and Newick files can be opened with software supporting those formats. No repository Python or R environment is required to download or inspect these files. This archive is outside the batch analysis runners.

## Limitations

Alignment-generation commands, tree-inference commands, substitution models, software versions and branch-support settings are not supplied in this archive. These files document the supplied alignments and trees; the archive alone does not provide a complete tree-inference workflow.

## License and Citation

Reuse follows the repository [data license](../DATA_LICENSE.md) and [notice](../NOTICE.md). Reference sequences retain their original attribution and terms. Cite this repository using [CITATION.cff](../CITATION.cff) and the associated manuscript when available.
