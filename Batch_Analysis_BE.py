#!/usr/bin/env python
# coding: utf-8

# In[1]:


# ============================================================================
# CELL 1: IMPORT LIBRARIES
# ============================================================================
# Run this cell first to load all necessary libraries

import pandas as pd
import re
from pathlib import Path
import os

# Make display() safe outside Jupyter
try:
    from IPython.display import display  # works if running in a notebook
except Exception:                        # plain python: define a fallback
    def display(x):
        try:
            print(x.to_string())
        except Exception:
            print(x)

print("✓ Libraries imported successfully!")
print("Ready to analyze base editing data.")


# In[13]:


# ============================================================================
# CELL 2: CONFIGURATION - EDIT THIS SECTION FOR YOUR EXPERIMENT
# ============================================================================
# This is the ONLY cell you need to modify for your analysis

# 1. SET YOUR DATA FOLDER PATH
#    Replace this with the path to your batch folder containing CRISPResso results
main_folder = r"/pub/madelk1/LargeBatch15/RawData/Renamed_Files/dedup_fastq/CRISPRessoBatch15Array"

# 2. SET YOUR CUTOFF THRESHOLD
#    Alleles below this percentage will be excluded from analysis
CUTOFF_PCT = 0.20

# 3. SET YOUR TARGET POSITIONS AND STRAND FOR EACH SAMPLE
#    Format: {'sample_identifier': (strand, target_position)}
#    - sample_identifier: unique text in the sample name (e.g., 'Site1', 'OFF1', etc.)
#    - strand: 'forward' for A→G editing, 'reverse' for T→C editing
#    - target_position: A-number from PAM-distal end (e.g., 5 means A5), or None for no target
#
#    Examples:
#    TARGET_CONFIG = {
#        'Site3': ('reverse', 5),    # Site1 samples: reverse strand, target at A5
#        'Site2': ('forward', 7),    # Site2 samples: forward strand, target at A7
#        'OFF1': ('reverse', None),  # OFF1 samples: reverse strand, no specific target
#        '5OT1': ('forward', None),  # 5OT1 samples: forward strand, no specific target
#    }
#
#    For purely off-target analysis with known strand:
#    TARGET_CONFIG = {
#        'sample': ('reverse', None),  # Will match all samples, reverse strand, no target
#    }
#
TARGET_CONFIG = {
    '_Site_3_OT_3OT1_':    ('forward', None),
    '_Site_5_OT_5OT1_':    ('reverse', None),
    '_Site_5_OT_5OT2_':    ('forward', None),
    '_Site_5_OT_5OT3_':    ('forward', None),
    '_Site_6_OT_6OT4_':    ('reverse', None),
    '_Site_6_OT_6OT6_':    ('reverse', None),
    '_Site_6_OT_6OT12_':    ('reverse', None),
    '_Site_EMX1_OT_EMX1OT1_':    ('reverse', None),
    '_Site_EMX1_OT_EMX1OT2_':    ('forward', None),
    '_Site_EMX1_OT_EMX1OT3_':    ('forward', None),
    '_Site_7.10_HEK_Site_3_OT_OT1_':    ('forward', None),
    '_Site_7.10_HEK_Site_3_OT_OT2_':    ('forward', None),
    '_Site_7.10_HEK_Site_3_OT_OT3_':    ('forward', None),
    '_Site_7.10_HEK_Site_4_OT_OT1_':    ('forward', None),
    '_Site_7.10_HEK_Site_4_OT_OT3_':    ('forward', None),
    '_Site_7.10_HEK_Site_4_OT_OT4_':    ('forward', None),
    'R_Loop_Site_1_':    ('reverse', None),
    'R_Loop_Site_2_':    ('reverse', None),
    'R_Loop_Site_3_':    ('forward', None),
    'R_Loop_Site_4_':    ('forward', None),
    'R_Loop_Site_5_':    ('forward', None),
    'ABE7.10_Site_17_':    ('forward', None),
}

# 4. EDIT TYPE - ABE (Adenine Base Editor)
EDIT_TYPE = "ABE"

# 5. OUTPUT FILE NAMES (optional - you can change these)
SUMMARY_OUTPUT = "Batch15editing_analysis_summary.csv"
POSITION_OUTPUT = "Batch15editing_positions_for_prism.csv"

print("="*80)
print("CONFIGURATION - ABE ANALYSIS")
print("="*80)
print(f"Data folder: {main_folder}")
print(f"Cutoff threshold: ≥{CUTOFF_PCT}%")
print(f"\nSample configurations:")
if TARGET_CONFIG:
    for identifier, (strand, position) in TARGET_CONFIG.items():
        target_str = f"A{position}" if position else "no target"
        print(f"  • '{identifier}' → {strand} strand, {target_str}")
else:
    print("  • No configuration specified!")
print(f"\nOutput files: {SUMMARY_OUTPUT}, {POSITION_OUTPUT}")
print("="*80)


# In[14]:


# ============================================================================
# CELL 3: HELPER FUNCTIONS
# ============================================================================
# These functions do the analysis - you don't need to modify this cell

def reverse_complement(seq):
    """Convert DNA sequence to reverse complement."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def find_sgrna_in_reference(reference_seq, sgrna_filename):
    """
    Find sgRNA position in reference sequence.
    Tries both forward and reverse complement.
    Handles both uppercase and lowercase.
    Returns: (actual_sgrna, start_position, end_position, is_reverse_complement) or None
    """
    # Convert to uppercase for matching
    reference_seq_upper = reference_seq.upper()
    sgrna_filename_upper = sgrna_filename.upper()
    
    # Try reverse complement first (common for reverse strand editing)
    actual_sgrna = reverse_complement(sgrna_filename_upper)
    if actual_sgrna in reference_seq_upper:
        start = reference_seq_upper.index(actual_sgrna)
        return actual_sgrna, start, start + len(actual_sgrna), True
    
    # Try forward strand
    if sgrna_filename_upper in reference_seq_upper:
        start = reference_seq_upper.index(sgrna_filename_upper)
        return sgrna_filename_upper, start, start + len(sgrna_filename_upper), False
    
    return None

def get_pam_distal_position(pos_from_5prime, sgrna_length, is_reverse_strand):
    """
    Convert 5'-based position to PAM-distal position (A numbering).
    
    For forward strand: PAM at 3' (right), count from left (5' end)
        Position 1 from left = A1, position 2 = A2, etc.
    
    For reverse strand: PAM at 3' (right), count from right (reverse 5' end)
        Position 1 from right = A1, position 2 from right = A2, etc.
    """
    if is_reverse_strand:
        # Count from right (reverse 5' end)
        return sgrna_length - pos_from_5prime + 1
    else:
        # Count from left (forward 5' end)  
        return pos_from_5prime

def find_config_for_sample(sample_name, config_dict):
    """
    Find strand and target position for a sample based on name patterns.
    Returns (strand, target_position) tuple or (None, None) if not found.
    """
    if not config_dict:
        return None, None
        
    for pattern, (strand, position) in config_dict.items():
        if pattern.lower() in sample_name.lower():
            return strand, position
    return None, None

def count_edits_in_window(aligned_seq, reference_seq, start, end, edit_from, edit_to):
    """
    Count specific base edits in a window.
    Returns: dictionary with position (1-indexed from 5' end) as key, True if edited
    """
    edits = {}
    for i in range(start, end):
        ref_base = reference_seq[i].upper()
        aligned_base = aligned_seq[i].upper()
        if ref_base == edit_from and aligned_base == edit_to:
            pos_1indexed = i - start + 1
            edits[pos_1indexed] = True
    return edits

def count_all_differences(aligned_seq, reference_seq, start, end):
    """
    Count all differences (substitutions or indels) in a window.
    Returns: number of differences
    """
    differences = 0
    for i in range(start, end):
        ref_base = reference_seq[i].upper()
        aligned_base = aligned_seq[i].upper()
        if ref_base != aligned_base or ref_base == '-' or aligned_base == '-':
            differences += 1
    return differences

print("✓ Helper functions loaded successfully!")


# In[15]:


# ============================================================================
# CELL 4: FIND AND VALIDATE INPUT FILES
# ============================================================================
# This cell finds all CRISPResso sample folders in your batch directory

os.chdir(main_folder)
sample_folders = [f for f in os.listdir(main_folder) if os.path.isdir(f)]

print("="*80)
print("SAMPLE DETECTION")
print("="*80)
print(f"Looking in: {main_folder}\n")
print(f"Found {len(sample_folders)} sample folder(s):\n")

for i, folder in enumerate(sample_folders, 1):
    print(f"  {i}. {folder}")

if len(sample_folders) == 0:
    print("\n⚠ WARNING: No sample folders found!")
    print("Make sure you've extracted the CRISPResso ZIP files.")
else:
    print(f"\n✓ Ready to analyze {len(sample_folders)} sample(s)")


# In[16]:


# ============================================================================
# CELL 5: RUN ANALYSIS ON ALL SAMPLES
# ============================================================================
# This cell processes all samples and calculates editing metrics

print("="*80)
print("BASE EDITING ANALYSIS - ABE")
print("="*80)
print(f"Cutoff: ≥{CUTOFF_PCT}%")
print("="*80)

# Storage for results
all_results = []
all_position_results = []
all_unintended_edits = []

# Process each sample folder
for folder in sample_folders:
    print(f"\n{'='*80}")
    print(f"PROCESSING: {folder}")
    print(f"{'='*80}")
    
    # Parse sample name - remove "CRISPResso_on_" prefix
    if folder.startswith("CRISPResso_on_"):
        sample_name = folder.replace("CRISPResso_on_", "")
    else:
        sample_name = folder
    
    print(f"Sample name: {sample_name}")
    
    # Find configuration for this sample
    detected_strand, target_pos_A = find_config_for_sample(sample_name, TARGET_CONFIG)
    
    if not detected_strand:
        print(f"  ⚠ No configuration found for this sample - SKIPPING")
        print(f"     Add a matching pattern to TARGET_CONFIG in Cell 2")
        continue
    
    # Find allele frequency table
    sample_path = os.path.join(main_folder, folder)
    allele_files = [f for f in os.listdir(sample_path) 
                   if f.startswith('Alleles_frequency_table_around_sgRNA_') 
                   and f.endswith('.txt')]
    
    if not allele_files:
        print(f"  ⚠ No allele frequency file found - SKIPPING")
        continue
    
    # Read data
    file_path = os.path.join(sample_path, allele_files[0])
    df = pd.read_csv(file_path, sep='\t')
    
    # Extract sgRNA sequence from filename (handle both upper and lowercase)
    sgrna_match = re.search(r'sgRNA_([ATCGatcg]+)', allele_files[0], re.IGNORECASE)
    if not sgrna_match:
        print(f"  ⚠ Could not extract sgRNA from filename - SKIPPING")
        continue
    
    sgrna_filename = sgrna_match.group(1)
    reference_seq = df.iloc[0]['Reference_Sequence']
    
    # Find sgRNA in reference
    sgrna_info = find_sgrna_in_reference(reference_seq, sgrna_filename)
    if not sgrna_info:
        print(f"  ⚠ Could not find sgRNA in reference sequence - SKIPPING")
        continue
    
    actual_sgrna, sgrna_start, sgrna_end, found_as_rc = sgrna_info
    sgrna_length = len(actual_sgrna)
    
    # Use configured strand (not auto-detected)
    if detected_strand == 'reverse':
        edit_from = 'T'
        edit_to = 'C'
        is_reverse_strand = True
        strand_label = "Reverse (T→C, template A→G)"
    elif detected_strand == 'forward':
        edit_from = 'A'
        edit_to = 'G'
        is_reverse_strand = False
        strand_label = "Forward (A→G)"
    else:
        print(f"  ⚠ Invalid strand configuration: {detected_strand} - SKIPPING")
        continue
    
    print(f"sgRNA: {actual_sgrna}")
    print(f"Length: {sgrna_length} nt")
    print(f"Configured strand: {strand_label}")
    print(f"Looking for: {edit_from}→{edit_to} edits")
    
    # Convert target position from A-numbering to 5' position
    if target_pos_A:
        if is_reverse_strand:
            target_pos_5prime = sgrna_length - target_pos_A + 1
        else:
            target_pos_5prime = target_pos_A
        
        print(f"Target: A{target_pos_A} (position {target_pos_5prime} from 5' end)")
        print(f"Target base: {actual_sgrna[target_pos_5prime-1]}→{edit_to}")
    else:
        target_pos_5prime = None
        print(f"Target: None (off-target analysis)")
    
    # Find all editable positions
    editable_positions_5prime = []
    for i in range(sgrna_length):
        if actual_sgrna[i].upper() == edit_from:
            editable_positions_5prime.append(i + 1)
    
    # PAM-distal positions (A numbering)
    editable_positions_A = [get_pam_distal_position(pos, sgrna_length, is_reverse_strand) 
                            for pos in editable_positions_5prime]
    print(f"Editable positions: {['A' + str(pos) for pos in sorted(editable_positions_A)]}")
    
    # Filter by cutoff
    total_reads_all = df['#Reads'].sum()
    df_filtered = df[df['%Reads'] >= CUTOFF_PCT].copy()
    total_reads_filtered = df_filtered['#Reads'].sum()
    
    print(f"Total reads (all alleles): {total_reads_all:,}")
    print(f"Total reads (≥{CUTOFF_PCT}%): {total_reads_filtered:,}")
    print(f"Alleles analyzed (≥{CUTOFF_PCT}%): {len(df_filtered)}/{len(df)}")
    
    # Initialize counters
    overall_editing_pct = 0.0
    overall_editing_reads = 0
    on_target_pct = 0.0
    on_target_reads = 0
    pure_target_pct = 0.0
    pure_target_reads = 0
    unintended_edits_pct = 0.0
    unintended_edits_reads = 0
    position_edit_counts = {pos: 0 for pos in editable_positions_5prime}
    position_edit_pcts = {pos: 0.0 for pos in editable_positions_5prime}
    
    # Track unintended edits by position and type
    unintended_by_position = {}
    
    # Analyze each allele
    for idx, row in df_filtered.iterrows():
        aligned_seq = row['Aligned_Sequence']
        n_reads = row['#Reads']
        pct_reads = row['%Reads']
        
        # Find all intended edits
        edits = count_edits_in_window(aligned_seq, reference_seq, 
                                      sgrna_start, sgrna_end, 
                                      edit_from, edit_to)
        
        # Check for unintended edits
        has_unintended = False
        for i in range(sgrna_start, sgrna_end):
            ref_base = reference_seq[i].upper()
            aligned_base = aligned_seq[i].upper()
            
            if ref_base == aligned_base:
                continue
            
            if ref_base == '-' or aligned_base == '-':
                continue
            
            if ref_base == edit_from and aligned_base == edit_to:
                continue
            
            has_unintended = True
            
            pos_5prime = i - sgrna_start + 1
            pos_A = get_pam_distal_position(pos_5prime, sgrna_length, is_reverse_strand)
            edit_type = f"{ref_base}to{aligned_base}"
            
            if pos_A not in unintended_by_position:
                unintended_by_position[pos_A] = {}
            if edit_type not in unintended_by_position[pos_A]:
                unintended_by_position[pos_A][edit_type] = {'reads': 0, 'pct': 0.0}
            
            unintended_by_position[pos_A][edit_type]['reads'] += n_reads
            unintended_by_position[pos_A][edit_type]['pct'] += pct_reads
        
        if has_unintended:
            unintended_edits_reads += n_reads
            unintended_edits_pct += pct_reads
        
        # Count position-specific edits
        for pos in edits.keys():
            if pos in position_edit_counts:
                position_edit_counts[pos] += n_reads
                position_edit_pcts[pos] += pct_reads
        
        # Overall editing
        if edits:
            overall_editing_reads += n_reads
            overall_editing_pct += pct_reads
        
        # On-target editing
        if target_pos_5prime and target_pos_5prime in edits:
            on_target_reads += n_reads
            on_target_pct += pct_reads
            
            # Purity
            total_diffs = count_all_differences(aligned_seq, reference_seq, 
                                               sgrna_start, sgrna_end)
            if total_diffs == 1:
                pure_target_reads += n_reads
                pure_target_pct += pct_reads
    
    print(f"\n--- Results ---")
    print(f"Overall editing (any {edit_from}→{edit_to}): {overall_editing_reads:,} reads ({overall_editing_pct:.2f}%)")
    if target_pos_5prime:
        print(f"On-target (A{target_pos_A}): {on_target_reads:,} reads ({on_target_pct:.2f}%)")
        print(f"Purity (only A{target_pos_A}): {pure_target_reads:,} reads ({pure_target_pct:.2f}%)")
    else:
        print(f"On-target: N/A (no target specified)")
        print(f"Purity: N/A (no target specified)")
    print(f"Unintended edits: {unintended_edits_reads:,} reads ({unintended_edits_pct:.2f}%)")
    
    # Display unintended edits
    if unintended_by_position:
        print(f"\n--- Unintended Edits Detected ---")
        print(f"Position | Edit Type | Reads    | %")
        print("---------|-----------|----------|-------")
        for pos_A in sorted(unintended_by_position.keys()):
            for edit_type, data in unintended_by_position[pos_A].items():
                print(f"   A{pos_A:2d}   | {edit_type:9s} | {data['reads']:8,} | {data['pct']:5.2f}%")
    
    # Display position profile
    print(f"\n--- Position Profile (PAM-distal numbering) ---")
    print(f"Position | Base | Type      | Editing %")
    print("---------|------|-----------|----------")
    
    position_data = []
    for pos_5prime in editable_positions_5prime:
        pos_A = get_pam_distal_position(pos_5prime, sgrna_length, is_reverse_strand)
        if target_pos_5prime:
            pos_type = "TARGET" if pos_5prime == target_pos_5prime else "BYSTANDER"
        else:
            pos_type = "N/A"
        edit_pct = position_edit_pcts[pos_5prime]
        position_data.append((pos_A, pos_5prime, pos_type, edit_pct))
    
    position_data.sort(key=lambda x: x[0])
    
    for pos_A, pos_5prime, pos_type, edit_pct in position_data:
        print(f"   A{pos_A:2d}   |  {actual_sgrna[pos_5prime-1]}   | {pos_type:9s} | {edit_pct:7.2f}%")
    
    # Store summary results
    all_results.append({
        'Sample': sample_name,
        'sgRNA': actual_sgrna,
        'sgRNA_Length': sgrna_length,
        'Strand': detected_strand,
        'Edit_Type': f"{edit_from}→{edit_to}",
        'Total_Reads_Analyzed': total_reads_filtered,
        'Overall_Editing_Reads': overall_editing_reads,
        'Overall_Editing_%': round(overall_editing_pct, 2),
        'On_Target_Reads': on_target_reads if target_pos_5prime else 'N/A',
        'On_Target_%': round(on_target_pct, 2) if target_pos_5prime else 'N/A',
        'Pure_Target_Reads': pure_target_reads if target_pos_5prime else 'N/A',
        'Purity_%': round(pure_target_pct, 2) if target_pos_5prime else 'N/A',
        'Target_Position': f"A{target_pos_A}" if target_pos_A else 'N/A'
    })
    
    # Store position-specific results
    for pos_5prime in editable_positions_5prime:
        pos_A = get_pam_distal_position(pos_5prime, sgrna_length, is_reverse_strand)
        if target_pos_5prime:
            pos_type = "TARGET" if pos_5prime == target_pos_5prime else "BYSTANDER"
        else:
            pos_type = "N/A"
        edit_pct = position_edit_pcts[pos_5prime]
        
        all_position_results.append({
            'Sample': sample_name,
            'Position': f"A{pos_A}",
            'Base': actual_sgrna[pos_5prime-1],
            'Strand': detected_strand,
            'Type': pos_type,
            'Editing_%': round(edit_pct, 2)
        })
    
    # Store unintended edits details
    for pos_A in unintended_by_position:
        for edit_type, data in unintended_by_position[pos_A].items():
            all_unintended_edits.append({
                'Sample': sample_name,
                'Position': f"A{pos_A}",
                'Unintended_Edit': edit_type,
                'Reads': data['reads'],
                'Percentage': round(data['pct'], 2)
            })

print(f"\n{'='*80}")
print(f"✓ ANALYSIS COMPLETE - Processed {len(all_results)} samples")
print(f"{'='*80}")


# In[17]:


# ============================================================================
# CELL 6: DISPLAY AND SAVE RESULTS
# ============================================================================
# View results and export to CSV files

# Create summary dataframe
print("="*80)
print("SUMMARY TABLE - OVERALL METRICS")
print("="*80)
print("This table shows overall editing, on-target, and purity for each sample\n")

results_df = pd.DataFrame(all_results)

# Fix the Edit_Type column to be CSV-friendly (replace arrow with "to")
if 'Edit_Type' in results_df.columns:
    results_df['Edit_Type'] = results_df['Edit_Type'].str.replace('→', 'to')

display(results_df)

# Separate position data by strand
print(f"\n{'='*80}")
print("POSITION EDITING PROFILES - FORMATTED FOR GRAPHPAD PRISM")
print("="*80)
print("Positions are PAM-distal (A-numbering)")
print("Separated by strand for easier analysis\n")

position_df = pd.DataFrame(all_position_results)

# Split into reverse and forward strand
reverse_positions = position_df[position_df['Strand'] == 'reverse'].copy()
forward_positions = position_df[position_df['Strand'] == 'forward'].copy()

# Create pivot tables for each strand
reverse_pivot = None
forward_pivot = None

if not reverse_positions.empty:
    print("--- REVERSE STRAND (T to C) ---")
    reverse_pivot = reverse_positions.pivot(index='Position', columns='Sample', values='Editing_%')
    reverse_pivot = reverse_pivot.round(2)
    # Sort by A-number (extract number from "A5" format)
    reverse_pivot.index = reverse_pivot.index.map(lambda x: int(x[1:]) if isinstance(x, str) and x.startswith('A') else x)
    reverse_pivot = reverse_pivot.sort_index()
    reverse_pivot.index = reverse_pivot.index.map(lambda x: f"A{x}")
    display(reverse_pivot)
    print()

if not forward_positions.empty:
    print("--- FORWARD STRAND (A to G) ---")
    forward_pivot = forward_positions.pivot(index='Position', columns='Sample', values='Editing_%')
    forward_pivot = forward_pivot.round(2)
    # Sort by A-number (extract number from "A5" format)
    forward_pivot.index = forward_pivot.index.map(lambda x: int(x[1:]) if isinstance(x, str) and x.startswith('A') else x)
    forward_pivot = forward_pivot.sort_index()
    forward_pivot.index = forward_pivot.index.map(lambda x: f"A{x}")
    display(forward_pivot)
    print()

# Display unintended edits table
print(f"\n{'='*80}")
print("UNINTENDED EDITS TABLE")
print("="*80)
print("Details of non-intended base substitutions\n")

if all_unintended_edits:
    unintended_df = pd.DataFrame(all_unintended_edits)
    unintended_df = unintended_df.sort_values(['Sample', 'Position'])
    display(unintended_df)
else:
    print("No unintended edits detected!")
    unintended_df = pd.DataFrame()

# Save files separately

# FILE 1: Editing analysis summary (overall metrics + unintended edits)
with open(SUMMARY_OUTPUT, 'w') as f:
    # Section 1: Overall metrics
    f.write("OVERALL METRICS\n")
    results_df.to_csv(f, index=False)
    f.write("\n\n")
    
    # Section 2: Unintended edits
    f.write("UNINTENDED EDITS DETAILS\n")
    if not unintended_df.empty:
        unintended_df.to_csv(f, index=False)
    else:
        f.write("No unintended edits detected\n")

# FILE 2: Position profiles for Prism (separate tables by strand)
with open(POSITION_OUTPUT, 'w') as f:
    if reverse_pivot is not None and not reverse_positions.empty:
        f.write("REVERSE STRAND (T to C)\n")
        reverse_pivot.to_csv(f)
        f.write("\n\n")
    
    if forward_pivot is not None and not forward_positions.empty:
        f.write("FORWARD STRAND (A to G)\n")
        forward_pivot.to_csv(f)

print(f"\n{'='*80}")
print("FILES SAVED")
print(f"{'='*80}")
print(f"✓ {SUMMARY_OUTPUT}")
print(f"  → Section 1: Overall metrics (total editing, on-target, purity)")
print(f"  → Section 2: Unintended edits details")
if not unintended_df.empty:
    print(f"      • Found {len(unintended_df)} unintended edit entries")
else:
    print(f"      • No unintended edits detected")

print(f"\n✓ {POSITION_OUTPUT}")
print(f"  → Position-by-position editing (PAM-distal A-numbering, separated by strand)")
if reverse_pivot is not None and not reverse_positions.empty:
    print(f"    • Reverse strand samples: {len(reverse_pivot.columns)}")
if forward_pivot is not None and not forward_positions.empty:
    print(f"    • Forward strand samples: {len(forward_pivot.columns)}")

print(f"\nFiles saved in: {os.getcwd()}")
print(f"{'='*80}")


# In[ ]:




