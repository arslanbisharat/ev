"""
Threshold Analysis Script
Analyzes different combinations of FDR and Fold Change thresholds
to determine optimal cutoffs for significant protein identification.
"""

import csv
import os

def load_results(filepath):
    """Load results from CSV file"""
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Convert to floats
    data = []
    for row in rows:
        data.append({
            'Gene': row['Gene'],
            'Accession': row['Accession'],
            'Description': row['Description'],
            'Log2FC': float(row['Log2FC']),
            'P_value': float(row['P_value']),
            'FDR': float(row['FDR']),
            'Cohens_d': float(row['Cohens_d'])
        })
    return data


def analyze_thresholds(data, fdr_thresholds, fc_thresholds):
    """
    Analyze how many proteins pass different threshold combinations

    Parameters:
    - data: list of protein results
    - fdr_thresholds: list of FDR cutoffs to test
    - fc_thresholds: list of |Log2FC| cutoffs to test

    Returns: list of results for each combination
    """
    results = []

    for fdr_cutoff in fdr_thresholds:
        for fc_cutoff in fc_thresholds:
            # Count proteins passing both criteria
            sig_proteins = [
                p for p in data
                if p['FDR'] < fdr_cutoff and abs(p['Log2FC']) > fc_cutoff
            ]

            # Also count directional changes
            upregulated = [p for p in sig_proteins if p['Log2FC'] > 0]
            downregulated = [p for p in sig_proteins if p['Log2FC'] < 0]

            results.append({
                'FDR_cutoff': fdr_cutoff,
                'FC_cutoff': fc_cutoff,
                'Total_significant': len(sig_proteins),
                'Upregulated': len(upregulated),
                'Downregulated': len(downregulated),
                'proteins': sig_proteins
            })

    return results


def print_threshold_matrix(results, title):
    """Print a matrix showing protein counts for different thresholds"""
    print(f"\n{'='*80}")
    print(f"{title}")
    print(f"{'='*80}\n")

    # Get unique FDR and FC thresholds
    fdr_vals = sorted(list(set(r['FDR_cutoff'] for r in results)))
    fc_vals = sorted(list(set(r['FC_cutoff'] for r in results)))

    # Print header
    header = "FDR \\ |Log2FC|"
    print(f"{header:<15}", end='')
    for fc in fc_vals:
        print(f"{fc:>10.2f}", end='')
    print()
    print("-" * (15 + 10*len(fc_vals)))

    # Print rows
    for fdr in fdr_vals:
        print(f"{fdr:<15.3f}", end='')
        for fc in fc_vals:
            # Find matching result
            match = [r for r in results if r['FDR_cutoff']==fdr and r['FC_cutoff']==fc]
            if match:
                count = match[0]['Total_significant']
                print(f"{count:>10}", end='')
        print()

    print()


def print_detailed_results(results, fdr_target, fc_target):
    """Print detailed info for a specific threshold combination"""
    match = [r for r in results
             if r['FDR_cutoff']==fdr_target and r['FC_cutoff']==fc_target]

    if not match:
        print(f"No results found for FDR={fdr_target}, FC={fc_target}")
        return

    result = match[0]

    print(f"\n{'='*80}")
    print(f"DETAILED RESULTS: FDR < {fdr_target}, |Log2FC| > {fc_target}")
    print(f"{'='*80}")
    print(f"\nTotal significant proteins: {result['Total_significant']}")
    print(f"  Upregulated (Sarc > Control): {result['Upregulated']}")
    print(f"  Downregulated (Sarc < Control): {result['Downregulated']}")

    # Show top 10 by effect size
    proteins = sorted(result['proteins'], key=lambda x: abs(x['Log2FC']), reverse=True)

    print(f"\nTop 10 by |Log2FC|:")
    print(f"{'Rank':<6} {'Gene':<12} {'Log2FC':<10} {'FDR':<12} {'Direction':<12}")
    print("-" * 60)

    for i, p in enumerate(proteins[:10], 1):
        direction = "Up" if p['Log2FC'] > 0 else "Down"
        print(f"{i:<6} {p['Gene']:<12} {p['Log2FC']:<10.3f} {p['FDR']:<12.2e} {direction:<12}")

    return result


def save_filtered_results(result, output_path):
    """Save filtered proteins to CSV"""
    with open(output_path, 'w', newline='') as f:
        if not result['proteins']:
            return

        fieldnames = ['Gene', 'Accession', 'Description', 'Log2FC', 'P_value', 'FDR', 'Cohens_d']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for p in result['proteins']:
            writer.writerow({
                'Gene': p['Gene'],
                'Accession': p['Accession'],
                'Description': p['Description'],
                'Log2FC': p['Log2FC'],
                'P_value': p['P_value'],
                'FDR': p['FDR'],
                'Cohens_d': p['Cohens_d']
            })


def main():
    # Paths
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    results_dir = os.path.join(project_dir, 'results')

    set1_path = os.path.join(results_dir, 'results_set1_old.csv')
    set2_path = os.path.join(results_dir, 'results_set2_fresh.csv')

    # Load data
    print("\nLoading results...")
    set1_data = load_results(set1_path)
    set2_data = load_results(set2_path)

    # Define thresholds to test
    fdr_thresholds = [0.001, 0.01, 0.05, 0.1, 0.2]
    fc_thresholds = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]

    # Analyze both sets
    set1_results = analyze_thresholds(set1_data, fdr_thresholds, fc_thresholds)
    set2_results = analyze_thresholds(set2_data, fdr_thresholds, fc_thresholds)

    # Print matrices
    print_threshold_matrix(set1_results, "SET 1 (OLD SAMPLES) - Significant Protein Counts")
    print_threshold_matrix(set2_results, "SET 2 (FRESH SAMPLES) - Significant Protein Counts")

    # Recommended thresholds to explore in detail
    recommended = [
        (0.05, 0.0),   # Standard FDR only
        (0.05, 0.5),   # Standard with moderate FC
        (0.05, 1.0),   # Standard with 2-fold change
        (0.1, 0.5),    # Relaxed FDR with moderate FC
        (0.1, 1.0),    # Relaxed FDR with 2-fold change
    ]

    print("\n" + "="*80)
    print("RECOMMENDED THRESHOLD COMBINATIONS")
    print("="*80)

    for fdr, fc in recommended:
        print(f"\n--- FDR < {fdr}, |Log2FC| > {fc} ---")
        set1_match = [r for r in set1_results if r['FDR_cutoff']==fdr and r['FC_cutoff']==fc][0]
        set2_match = [r for r in set2_results if r['FDR_cutoff']==fdr and r['FC_cutoff']==fc][0]

        print(f"Set 1: {set1_match['Total_significant']} proteins "
              f"({set1_match['Upregulated']} up, {set1_match['Downregulated']} down)")
        print(f"Set 2: {set2_match['Total_significant']} proteins "
              f"({set2_match['Upregulated']} up, {set2_match['Downregulated']} down)")

        # Calculate overlap
        set1_genes = set(p['Gene'] for p in set1_match['proteins'])
        set2_genes = set(p['Gene'] for p in set2_match['proteins'])
        overlap = len(set1_genes & set2_genes)

        print(f"Overlap: {overlap} proteins")

    # Show detailed results for a specific combination
    print_detailed_results(set1_results, 0.05, 0.5)
    print_detailed_results(set2_results, 0.05, 0.5)

    # Save filtered results with recommended threshold
    print("\n" + "="*80)
    print("SAVING FILTERED RESULTS")
    print("="*80)

    # Use FDR < 0.05, |Log2FC| > 0.5 as default
    set1_filtered = [r for r in set1_results if r['FDR_cutoff']==0.05 and r['FC_cutoff']==0.5][0]
    set2_filtered = [r for r in set2_results if r['FDR_cutoff']==0.05 and r['FC_cutoff']==0.5][0]

    save_filtered_results(set1_filtered, os.path.join(results_dir, 'results_set1_filtered_FDR0.05_FC0.5.csv'))
    save_filtered_results(set2_filtered, os.path.join(results_dir, 'results_set2_filtered_FDR0.05_FC0.5.csv'))

    print(f"\nSaved filtered results to:")
    print(f"  - results_set1_filtered_FDR0.05_FC0.5.csv ({set1_filtered['Total_significant']} proteins)")
    print(f"  - results_set2_filtered_FDR0.05_FC0.5.csv ({set2_filtered['Total_significant']} proteins)")

    # Calculate and save overlap
    set1_genes_filtered = set(p['Gene'] for p in set1_filtered['proteins'])
    set2_genes_filtered = set(p['Gene'] for p in set2_filtered['proteins'])
    overlap_filtered = set1_genes_filtered & set2_genes_filtered

    overlap_proteins = [p for p in set1_filtered['proteins'] if p['Gene'] in overlap_filtered]

    if overlap_proteins:
        save_filtered_results({'proteins': overlap_proteins},
                            os.path.join(results_dir, 'storage_stable_proteins_filtered_FDR0.05_FC0.5.csv'))
        print(f"  - storage_stable_proteins_filtered_FDR0.05_FC0.5.csv ({len(overlap_proteins)} proteins)")

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print("\nRECOMMENDATION:")
    print("  For proteomics, common thresholds are:")
    print("    - FDR < 0.05 AND |Log2FC| > 0.5 (1.4-fold change)")
    print("    - FDR < 0.05 AND |Log2FC| > 1.0 (2-fold change)")
    print("\n  Choose based on your research question:")
    print("    - Discovery phase: Use FDR < 0.1, |Log2FC| > 0.5")
    print("    - Validation phase: Use FDR < 0.05, |Log2FC| > 1.0")
    print("="*80)


if __name__ == "__main__":
    main()
