import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

def analyze_sequence_spikes(filepath, threshold, fr, sample):
    """Scans for spikes in second highest peaks"""
    print(f"Scanning for mutation spikes, threshold: {threshold:.0%}...\n")
    record = SeqIO.read(filepath, "abi")
    trace_data = record.annotations['abif_raw']
    
    fwo = trace_data['FWO_1'].decode('utf-8')
    channels = {base: trace_data[f'DATA{9+i}'] for i, base in enumerate(fwo)}
    peak_locations = trace_data['PLOC1']
    
    positions = []
    noise_ratios = []
    flagged_mutations = []
    
    # Trim messy ends
    if sample == "C7":
        start_trim = 15
        end_trim = 2000
        # if fr == "F":
        #     start_trim = 30
        #     end_trim = 670
        # else:
        #     start_trim = 50
        #     end_trim = 653
    else:
        start_trim = 55
        end_trim = 714
    end_index = min(end_trim, len(peak_locations) - 1)
    valid_peak_locations = peak_locations[start_trim:end_index-15]
    
    for bp_index, loc in enumerate(valid_peak_locations, start=start_trim+1):
        signals = {
            'A': channels['A'][loc],
            'C': channels['C'][loc],
            'G': channels['G'][loc],
            'T': channels['T'][loc]
        }
        
        sorted_bases = sorted(signals.items(), key=lambda item: item[1], reverse=True)
        primary_base, primary_signal = sorted_bases[0]
        secondary_base, secondary_signal = sorted_bases[1]
        
        if primary_signal > 0:
            ratio = secondary_signal / primary_signal
            positions.append(bp_index)
            noise_ratios.append(ratio)
            
            # Flag if the secondary peak exceeds ratio threshold
            if ratio >= threshold:
                flagged_mutations.append({
                    'Base_Pair_Num': bp_index,
                    'Primary_Base': primary_base,
                    'Secondary_Base': secondary_base,
                    'Ratio': ratio
                })

    return positions, noise_ratios, flagged_mutations

def plot_spike_analysis(base_path, target_file, duration, threshold, fr, sample):
    full_path = os.path.join(base_path, target_file)

    if os.path.exists(full_path):
        positions, ratios, mutations = analyze_sequence_spikes(full_path, threshold, fr, sample)
        
        print(f"--- Results for {duration}s UV Exposure ({target_file}) ---")
        if not mutations:
            print("No spikes found exceeding threshold")
        else:
            print(f"Found {len(mutations)} potential mutation sites!")
            for m in mutations:
                print(f" - BP #{m['Base_Pair_Num']}: 1) {m['Primary_Base']}, 2) {m['Secondary_Base']}, Ratio) {m['Ratio']:.2%})")
        
        # Plot noise ratio across the entire sequence
        plt.figure(figsize=(12, 5))
        plt.plot(positions, ratios, color='gray', linewidth=0.8, zorder=1)

        # Highlight flagged mutations in red
        spike_x = [m['Base_Pair_Num'] for m in mutations]
        spike_y = [m['Ratio'] for m in mutations]
        plt.scatter(spike_x, spike_y, color='red', zorder=5, label='Flagged Spikes')
        plt.axhline(y=threshold, color='r', linestyle='--', alpha=0.5, label=f'{threshold:.0%} Threshold')
        plt.title(f"Noise Ratio Across Sequence Length (Sample {sample}{fr} {duration}s UV Exposure)")
        plt.xlabel("Base Pair Position")
        plt.ylabel("Secondary/Primary Peak Ratio")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"figs/{sample}{fr}-control-{threshold:.0%}.svg")
        plt.show()
    else:
        print(f"Could not find {full_path}")

def plot_mutation_trends(file_mapping, base_path, threshold, fr, sample):
    """
    Iterates through all sequence files, counts the number of mutations 
    detected at the given threshold, and plots the trend.
    """
    results = []
    print(f"Scanning all files for mutations (Threshold set to {threshold:.0%})...\n")

    for filename, exposure_time in file_mapping.items():
        full_path = os.path.join(base_path, filename)
        if not os.path.exists(full_path):
            print(f"Warning: Could not find {filename}. Skipping...")
            continue
            
        _, _, mutations = analyze_sequence_spikes(full_path, threshold, fr, sample)
        mutation_count = len(mutations)
        results.append({
            "File": filename,
            "UV_Time_Seconds": exposure_time,
            "Mutation_Count": mutation_count
        })
        print(f"Processed {filename} (UV: {exposure_time}s) | Mutations Found: {mutation_count}")

    df = pd.DataFrame(results)
    df = df.sort_values(by="UV_Time_Seconds")

    plt.figure(figsize=(10, 6))
    plt.plot(df['UV_Time_Seconds'], df['Mutation_Count'], marker='o', linestyle='-', color='purple', markersize=8)
    plt.fill_between(df['UV_Time_Seconds'], df['Mutation_Count'], color='purple', alpha=0.1)
    plt.title(f"Mutation Count vs. UV Light Exposure Time\n(Sample {sample}{fr}, Detection Threshold: {threshold:.0%})", fontsize=14)
    plt.xlabel("UV Exposure Time (Seconds)", fontsize=12)
    plt.ylabel("Number of Detected Mutations (Secondary Peak Spikes)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.ylim(bottom=0) 
    plt.tight_layout()
    plt.savefig(f"figs/trend-{sample}{fr}-{threshold:.0%}-untrimmed.svg")
    plt.show()


if __name__ == "__main__":
    FR = "F"
    SAMPLE = "C8"
    file_first_mapping = {
        f"Kangas7925_G11-C_R1-16S-rRNA-seq{FR}.ab1": 0,
        f"Kangas7925_G11-1_R1-16S-rRNA-seq{FR}.ab1": 1,
        # f"Kangas7925_G11-2_R1-16S-rRNA-seq{FR}.ab1": 2,
        f"Kangas7925_G11-4_R1-16S-rRNA-seq{FR}.ab1": 4,
        f"Kangas7925_G11-8_R2-16S-rRNA-seq{FR}.ab1": 8,
        f"Kangas7925_G11-16_R1-16S-rRNA-seq{FR}.ab1": 16
    }
    file_second_mapping = {
        f"Kangas5125_C7_G11_C_R1-16S-rRNA-seq{FR}.ab1": 0,
        f"Kangas5125_C7_G11_3_R2-16S-rRNA-seq{FR}.ab1": 3,
        f"Kangas5125_C7_G11_7_R2-16S-rRNA-seq{FR}.ab1": 7,
        f"Kangas5125_C7_G11_15_R1-16S-rRNA-seq{FR}.ab1": 15,
        f"Kangas5125_C7_G11_30_R1-16S-rRNA-seq{FR}.ab1": 30,
        f"Kangas5125_C7_G11_60_R1-16S-rRNA-seq{FR}.ab1": 60
    }

    # Change these two variables if we want to analyze First trials or Second trials
    if SAMPLE == "C8":
        file_mapping = file_first_mapping
        base_path = "FirstSeq/ab1_files/"
    else:
        file_mapping = file_second_mapping
        base_path = "SecondSeq/ab1_files/"

    DETECTION_THRESHOLD = 0.25

    # Uncomment to analyze one file
    # filename = f"Kangas7925_G11-C_R1-16S-rRNA-seq{FR}.ab1"
    # plot_spike_analysis(base_path, filename, 0, threshold=DETECTION_THRESHOLD, fr=FR, sample=SAMPLE)

    # Analyze all files for one mapping
    # for filename, exposure_time in file_mapping.items():
    #     plot_spike_analysis(base_path, filename, exposure_time, threshold=DETECTION_THRESHOLD, fr=FR, sample=SAMPLE)

    plot_mutation_trends(file_mapping, base_path, threshold=DETECTION_THRESHOLD, fr=FR, sample=SAMPLE)