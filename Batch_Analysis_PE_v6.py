#!/usr/bin/env python
# coding: utf-8

# ============================================================================
# CELL 1: IMPORT LIBRARIES
# ============================================================================

import pandas as pd
import re
import os

try:
    from IPython.display import display
except Exception:
    def display(x):
        try:    print(x.to_string())
        except: print(x)

print("✓ Libraries imported successfully!")


# ============================================================================
# CELL 2: CONFIGURATION — EDIT THIS SECTION
# ============================================================================

# 1. PATH TO YOUR CRISPResso OUTPUT FOLDER
main_folder = r"/pub/madelk1/LargeBatch15/RawData/Renamed_Files/CRISPRessoBatch15PaperAnalysis_PEsamples"

# 2. OUTPUT FILE NAMES
SUMMARY_OUTPUT   = "Batch15updated_analysis_summary.csv"
PRISM_OUTPUT     = "Batch15updated_outcomes_for_prism.csv"
import datetime
_ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
FIDELITY_OUTPUT  = f"RTT_fidelity_analysis_{_ts}.csv"

# 3. NOT USED CURRENTLY (was leftover from previous script to analyze fidelity) RTT SEQUENCE (DNA 5'->3' as in pegRNA; RC will be searched in amplicon)
RTT_SEQ = "GGAAAAGCGATCAAGGT"

print("=" * 80)
print("CONFIGURATION — PRIME EDITING ANALYSIS (CRISPResso2 native PE output)")
print("=" * 80)
print(f"Data folder: {main_folder}")
print(f"Output:      {SUMMARY_OUTPUT}, {PRISM_OUTPUT}")
print("=" * 80)


# ============================================================================
# CELL 3: HELPER — PARSE CRISPResso_quantification_of_editing_frequency.txt
# ============================================================================

def parse_quantification_file(filepath):
    """
    Parse CRISPResso2's CRISPResso_quantification_of_editing_frequency.txt.

    In PE mode this is a tab-separated table with one row per amplicon class:
        Amplicon | Unmodified% | Modified% | Reads_in_input |
        Reads_aligned_all_amplicons | Reads_aligned | Unmodified |
        Modified | Discarded | Insertions | Deletions | Substitutions | ...

    Row labels: Reference, Prime-edited, Scaffold-incorporated

    Returns a dict of extracted values, or None if the file can't be parsed.
    """
    try:
        df = pd.read_csv(filepath, sep='\t', index_col=0)
        # Normalize index just in case of whitespace differences
        df.index = df.index.str.strip()
    except Exception:
        return None

    def get(row_label, col, as_int=True):
        try:
            val = df.loc[row_label, col]
            return int(val) if as_int else float(val)
        except Exception:
            return None

    result = {}

    # Total input reads and aligned reads (same for all rows — take from Reference)
    result['reads_total']       = get('Reference', 'Reads_in_input')
    result['reads_aligned_all'] = get('Reference', 'Reads_aligned_all_amplicons')

    # Per-amplicon aligned read counts
    result['ref_reads'] = get('Reference',              'Reads_aligned')
    result['pe_reads']  = get('Prime-edited',           'Reads_aligned')
    result['sc_reads']  = get('Scaffold-incorporated',  'Reads_aligned')

    # Per-amplicon unmodified (no indels/subs within that amplicon's window)
    result['ref_unmod'] = get('Reference',              'Unmodified')
    result['pe_unmod']  = get('Prime-edited',           'Unmodified')

    # Per-amplicon modified (= has indels or substitutions beyond the intended edit)
    result['ref_modified'] = get('Reference',     'Modified')
    result['pe_modified']  = get('Prime-edited',  'Modified')

    # Unmodified% and Modified% directly from file (as % of that amplicon's aligned reads)
    result['ref_unmod_pct']  = get('Reference',     'Unmodified%', as_int=False)
    result['ref_mod_pct']    = get('Reference',     'Modified%',   as_int=False)
    result['pe_unmod_pct']   = get('Prime-edited',  'Unmodified%', as_int=False)
    result['pe_mod_pct']     = get('Prime-edited',  'Modified%',   as_int=False)

    # Indel breakdown on reference amplicon (= NHEJ)
    result['ref_indel_reads'] = get('Reference', 'Modified')   # modified on ref = indels

    # Indel breakdown on PE amplicon (= PE + extra indel)
    result['pe_indel_reads']  = get('Prime-edited', 'Modified')

    return result


def pct_of_aligned(reads, reads_aligned_all):
    """Calculate percentage relative to all aligned reads."""
    if reads is None or reads_aligned_all is None or reads_aligned_all == 0:
        return None
    return round(100.0 * reads / reads_aligned_all, 2)


print("✓ Helper functions loaded.")


# ============================================================================
# CELL 4: FIND SAMPLE FOLDERS
# ============================================================================

os.chdir(main_folder)
sample_folders = sorted([
    f for f in os.listdir(main_folder)
    if os.path.isdir(os.path.join(main_folder, f))
])

print("=" * 80)
print("SAMPLE DETECTION")
print("=" * 80)
print(f"Looking in: {main_folder}\n")
print(f"Found {len(sample_folders)} folder(s):\n")
for i, f in enumerate(sample_folders, 1):
    print(f"  {i:>3}. {f}")
print(f"\n✓ Ready to process {len(sample_folders)} folder(s)")


# ============================================================================
# CELL 5: RUN ANALYSIS
# ============================================================================

print("\n" + "=" * 80)
print("PRIME EDITING ANALYSIS — reading CRISPResso2 native PE outputs")
print("=" * 80)

all_results = []

for folder in sample_folders:
    print(f"\n{'─'*80}")
    print(f"FOLDER: {folder}")

    sample_name = (folder.replace("CRISPResso_on_", "")
                   if folder.startswith("CRISPResso_on_") else folder)
    print(f"Sample: {sample_name}")

    folder_path = os.path.join(main_folder, folder)

    # ── 1. Parse quantification file ─────────────────────────────────────────
    quant_file = os.path.join(folder_path,
                              "CRISPResso_quantification_of_editing_frequency.txt")
    if not os.path.exists(quant_file):
        print(f"  ⚠ CRISPResso_quantification_of_editing_frequency.txt not found — SKIPPING")
        continue

    q = parse_quantification_file(quant_file)
    if q is None:
        print(f"  ⚠ Could not parse quantification file — SKIPPING")
        continue

    reads_total   = q['reads_total']   or 0
    reads_aligned = q['reads_aligned_all'] or 0

    pe_reads  = q['pe_reads']  or 0
    sc_reads  = q['sc_reads']  or 0
    ref_reads = q['ref_reads'] or 0

    # Recompute % of aligned reads (more consistent across samples)
    pe_pct  = pct_of_aligned(pe_reads,  reads_aligned)
    sc_pct  = pct_of_aligned(sc_reads,  reads_aligned)
    ref_pct = pct_of_aligned(ref_reads, reads_aligned)

    # Indels:
    #   - ref_indel: indels that landed on the reference amplicon  → NHEJ/unintended
    #   - pe_indel:  indels within prime-edited reads              → PE + indel byproduct
    ref_indel_reads = q['ref_indel_reads'] or 0
    pe_indel_reads  = q['pe_indel_reads']  or 0

    ref_indel_pct = pct_of_aligned(ref_indel_reads, reads_aligned)
    pe_indel_pct  = pct_of_aligned(pe_indel_reads,  reads_aligned)

    # Total indels (sum of both classes)
    total_indel_reads = ref_indel_reads + pe_indel_reads
    total_indel_pct   = pct_of_aligned(total_indel_reads, reads_aligned)

    # Antoniou et al. 2025 definitions:
    #   Prime edits all = precise PE + PE reads co-occurring with indels
    prime_edits_all_reads = pe_reads + pe_indel_reads
    prime_edits_all_pct   = pct_of_aligned(prime_edits_all_reads, reads_aligned)

    #   Scaffold of PE = scaffold reads / (PE reads + scaffold reads)
    #   i.e. scaffold contamination rate among all editing products
    sc_of_pe_denom = pe_reads + sc_reads
    sc_of_pe_pct   = round(100.0 * sc_reads / sc_of_pe_denom, 2) if sc_of_pe_denom > 0 else None

    # ── 2. Print summary ──────────────────────────────────────────────────────
    print(f"\n  Total reads:        {reads_total:>10,}")
    print(f"  Aligned reads:      {reads_aligned:>10,}")
    print(f"\n  {'Category':<35}  {'Reads':>9}  {'% aligned':>10}")
    print(f"  {'─'*35}  {'─'*9}  {'─'*10}")
    print(f"  {'Prime-edited (perfect PE)':<35}  {pe_reads:>9,}  {pe_pct or 0:>9.2f}%")
    print(f"  {'Scaffold-incorporated':<35}  {sc_reads:>9,}  {sc_pct or 0:>9.2f}%")
    print(f"  {'Reference (unmodified)':<35}  {ref_reads:>9,}  {ref_pct or 0:>9.2f}%")
    print(f"  {'  └─ indels on reference':<35}  {ref_indel_reads:>9,}  {ref_indel_pct or 0:>9.2f}%")
    print(f"  {'  └─ indels on PE allele':<35}  {pe_indel_reads:>9,}  {pe_indel_pct or 0:>9.2f}%")
    print(f"  {'  └─ total indels':<35}  {total_indel_reads:>9,}  {total_indel_pct or 0:>9.2f}%")

    # ── 3. Store result ───────────────────────────────────────────────────────
    all_results.append({
        'Sample':                    sample_name,
        'Reads_total':               reads_total,
        'Reads_aligned':             reads_aligned,

        # Core PE outcomes (% of aligned reads)
        'PE_%':                      pe_pct,
        'Scaffold_incorporated_%':   sc_pct,
        'Reference_unmodified_%':    ref_pct,

        # Indels (% of aligned reads)
        'Indels_on_reference_%':     ref_indel_pct,
        'Indels_on_PE_allele_%':     pe_indel_pct,
        'Total_indels_%':            total_indel_pct,

        # Antoniou et al. 2025 definitions
        'Prime_edits_all_%':         prime_edits_all_pct,
        'Scaffold_of_PE_%':          sc_of_pe_pct,

        # Raw read counts
        'PE_reads':                  pe_reads,
        'Scaffold_reads':            sc_reads,
        'Reference_reads':           ref_reads,
        'Indel_ref_reads':           ref_indel_reads,
        'Indel_PE_reads':            pe_indel_reads,
    })

print(f"\n{'='*80}")
print(f"✓ ANALYSIS COMPLETE — {len(all_results)} sample(s) processed")
print(f"{'='*80}")


# ============================================================================
# CELL 6: DISPLAY AND SAVE RESULTS
# ============================================================================

results_df = pd.DataFrame(all_results)

# ── Full summary table ────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("FULL SUMMARY TABLE")
print("=" * 80)
display(results_df)

# ── Prism-ready tables (% of aligned reads) ───────────────────────────────────
pct_cols = ['PE_%', 'Scaffold_incorporated_%', 'Reference_unmodified_%',
            'Total_indels_%', 'Indels_on_reference_%', 'Indels_on_PE_allele_%',
            'Prime_edits_all_%', 'Scaffold_of_PE_%']

print(f"\n{'='*80}")
print("PRISM-READY TABLE  (rows = outcomes, columns = samples, values = % aligned reads)")
print("=" * 80)

prism_df = None
if not results_df.empty:
    prism_df = (results_df[['Sample'] + pct_cols]
                .set_index('Sample')
                .T
                .round(2))
    display(prism_df)

# ── Save files ────────────────────────────────────────────────────────────────
with open(SUMMARY_OUTPUT, 'w', newline='') as f:
    f.write("PRIME EDITING OUTCOME SUMMARY\n")
    results_df.to_csv(f, index=False)

with open(PRISM_OUTPUT, 'w', newline='') as f:
    f.write("OUTCOME PERCENTAGES  (rows = outcomes, columns = samples)\n")
    if prism_df is not None:
        prism_df.to_csv(f)

print(f"\n{'='*80}")
print("FILES SAVED")
print(f"{'='*80}")
print(f"✓ {SUMMARY_OUTPUT}  — full metrics table (all samples × all columns)")
print(f"✓ {PRISM_OUTPUT}     — Prism-ready pivot (outcomes × samples, % aligned reads)")
print(f"\nSaved in: {os.getcwd()}")
print("=" * 80)


# ============================================================================
# CELL 7: RTT FIDELITY ANALYSIS — skipped for this batch (per-sample RTT seqs)
# ============================================================================


print("\n✓ Analysis complete.")
