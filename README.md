
# a2m2conserv
This script calculates conservation scores for each position in a multiple sequence alignment (MSA) in the **a2m format**. It uses **Shannon entropy** to evaluate residue frequency distributions and visualizes the results through: a heatmap of per-position conservation scores and a heatmap of Amino acid propensities at each alignment position. Positions that are mostly gaps (i.e., columns with fewer than 6 amino acids) are penalized and assigned a conservation score of 0.

# Usage
To run the script, make sure to replace the file paths for msa_path, output_csv_path, and output_plot_path with your actual paths. 
Then, simply execute the script:
> python3 a2m2conserv.py

# Requirements
Before running this script, you need to install the following libraries:
- Biopython: For handling sequence alignments.
- NumPy: For numerical computations.
- Pandas: For data manipulation and saving results.
- Seaborn: For creating heatmaps and visualizations.
- Matplotlib: For plotting and saving the heatmap images.

You can install these dependencies using pip:
> pip install biopython numpy pandas seaborn matplotlib

# Input
The script requires an MSA file in a2m format to analyze. The path to this file should be specified in the msa_path variable.

Example:
msa_path = "/path/to/your/alignment.a2m"
msa_format = "fasta"

# Output
The script generates two output files:

1. CSV File: A CSV file containing conservation scores and residue frequencies for each position in the alignment. The file is saved at the path specified by the output_csv_path variable.
Example:
output_csv_path = "/path/to/output/conservation_scores.csv"

2. Heatmap Image: A heatmap showing the conservation scores and amino acid frequencies across positions in the alignment. The plot is saved at the path specified by the output_plot_path variable.

output_plot_path = "/path/to/output/conservation_heatmap.png"
Conservation Strip Heatmap: A strip heatmap focused solely on conservation scores for each position. This is saved at a path derived from the original heatmap path.


Calculate Conservation Score
The conservation score is computed as:
conservation = 1 - (entropy / max_entropy)
This score is then penalized for positions with high gap fractions.


# License
This script is open source and available under the MIT license. See the LICENSE file for more details.
