# === IMPORTS ===
# Required libraries for sequence handling, data manipulation, and plotting
from Bio import AlignIO
from collections import Counter
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# === CONFIGURATION ===
# Paths for input MSA file and output files (CSV and plot image)
msa_path = "/media/elhabashy/Elements/collaboration/Andrei_project/ARCFU/evcomplex/msa/Y1502_ARCFU_1-68_b0.3/align/Y1502_ARCFU_1-68_b0.3.a2m"
msa_format = "fasta"

output_csv_path = "/media/elhabashy/Elements/collaboration/Andrei_project/ARCFU/plots/conservation_with_frequencies_fixed.csv"
output_plot_path = "/media/elhabashy/Elements/collaboration/Andrei_project/ARCFU/plots/conservation_with_frequencies2.png"

# === LOAD MSA ALIGNMENT ===
# Read the alignment file in FASTA format
alignment = AlignIO.read(msa_path, msa_format)
n_seqs = len(alignment)  # Total number of sequences
alignment_length = alignment.get_alignment_length()  # Number of positions (columns)

# === DEFINE SYMBOL SETS ===
# Standard amino acids and gap symbols (both - and .)
aa_list = list("ACDEFGHIKLMNPQRSTVWY")
gap_symbols = ['-', '.']
all_symbols = aa_list + gap_symbols

# Container for results
records = []

# === ANALYZE EACH COLUMN OF THE ALIGNMENT ===
for i in range(alignment_length):
    col = alignment[:, i]  # Get column i (residues at position i across all sequences)
    counts = Counter(col)  # Count occurrences of each residue/gap

    # Calculate the number of non-gap residues
    total_non_gap = sum(counts.get(aa, 0) for aa in aa_list)
    gap_count = sum(counts.get(gap, 0) for gap in gap_symbols)
    gap_fraction = gap_count / n_seqs  # Fraction of sequences with a gap at this position

    # === CALCULATE FREQUENCIES ===
    freqs = {}
    for aa in aa_list:
        # Normalize each residue count to the total number of non-gap residues
        freqs[aa] = counts.get(aa, 0) / total_non_gap if total_non_gap > 0 else 0.0
    freqs['gap'] = gap_fraction  # Gap frequency normalized to all sequences

    # === CALCULATE CONSERVATION SCORE ===
    if total_non_gap < 6:
        # Not enough data to assess conservation
        conservation = 0.0
    else:
        # Use Shannon entropy to measure variability at this position
        non_gap_freqs = {aa: counts.get(aa, 0) / total_non_gap for aa in aa_list}
        entropy = -sum(p * np.log2(p) for p in non_gap_freqs.values() if p > 0)
        max_entropy = np.log2(len(aa_list))  # Max possible entropy with 20 amino acids
        conservation = 1 - (entropy / max_entropy)  # Higher value = more conserved
        conservation *= (1 - gap_fraction)  # Penalize positions with many gaps

    # === RECORD RESULTS ===
    row = {
        'i': i + 1,  # 1-based indexing for positions
        'A_i': col[0],  # Reference residue from first sequence
        'conservation': round(conservation, 3)  # Round for cleaner output
    }
    # Add rounded frequencies for each amino acid
    row.update({aa: round(freqs[aa], 3) for aa in aa_list})
    row['-'] = round(freqs['gap'], 3)  # Include gap frequency

    records.append(row)

# === CREATE A DATAFRAME AND SAVE AS CSV ===
df = pd.DataFrame(records)
df.to_csv(output_csv_path, index=False)
print(f"CSV saved to: {output_csv_path}")
print(df.head())  # Print first few rows for verification

# === PREPARE DATA FOR HEATMAP ===
# Select only residue columns (excluding conservation and A_i)
residue_columns = ['-', 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

# Transpose so residues are on the y-axis and positions on the x-axis
heatmap_data = df.set_index('i')[residue_columns].T

# === PLOTTING HEATMAP ===
plt.figure(figsize=(20, 6))  # Set the size of the figure
sns.set_theme(style='white')  # Minimalist background
sns.set_context("talk", font_scale=1.2)  # Font scaling for presentation or publication

# Draw the heatmap using seaborn
ax = sns.heatmap(
    heatmap_data,
    cmap='Blues',  # Color palette
    cbar_kws={'label': 'Score'},  # Label for colorbar
    linewidths=0.5,  # Thin lines between cells
    linecolor='gray',
    xticklabels=True,  # Show position numbers
    yticklabels=residue_columns  # Show residue labels on y-axis
)

# === LABELS AND TITLES ===
plt.xlabel("Position i", fontsize=14)
plt.ylabel("Residue", fontsize=14)
plt.title("Conservation Score", fontsize=16)

# Rotate tick labels for better readability
plt.xticks(rotation=45, ha='center', fontsize=12)
plt.yticks(rotation=0, fontsize=12)


# Final layout tweaks
plt.tight_layout()
plt.savefig(output_plot_path, dpi=300, bbox_inches='tight')  # High-resolution output
plt.close()

print(f"Plot saved to: {output_plot_path}")

# === PLOT STRIP HEATMAP FOR CONSERVATION ONLY ===
# Prepare data: 2D array with shape (1, N positions)
conservation_data = df[['conservation']].T  # shape: (1, alignment_length)

plt.figure(figsize=(20, 1.5))  # Shorter height for strip
sns.set_theme(style='white')
sns.set_context("talk", font_scale=1.2)

ax = sns.heatmap(
    conservation_data,
    cmap='Blues',  # Use a distinct color map for contrast
    cbar_kws={'label': 'Conservation Score'},
    xticklabels=True,
    yticklabels=['Conservation'],
    linewidths=0.5,
    linecolor='gray'
)

plt.xlabel("Position i", fontsize=12)
#plt.xticks(rotation=45, ha='center', fontsize=10)
plt.yticks(rotation=0, fontsize=12)

positions = df['i'].values
ax.set_xticks(np.arange(len(positions)) + 0.5)  # positions centered
ax.set_xticklabels(positions, rotation=45, ha='center', fontsize=10)

plt.tight_layout()
strip_plot_path = output_plot_path.replace(".png", "_strip.png")
plt.savefig(strip_plot_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Conservation strip heatmap saved to: {strip_plot_path}")
