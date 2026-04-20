import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

def analyze_variance(filepath):
    record = SeqIO.read(filepath, "abi")
    trace_data = record.annotations['abif_raw']
    
    #FWO_1 contains a string like b'GATC' showing order of DATA9-DATA12
    fwo = trace_data['FWO_1'].decode('utf-8')
    channels = {}
    for i, base in enumerate(fwo):
        channels[base] = trace_data[f'DATA{9+i}']
        
    peak_locations = trace_data['PLOC1']
    
    noise_ratios = []
    
    for loc in peak_locations:
        #Fetch raw data :)
        signals = [
            channels['A'][loc],
            channels['C'][loc],
            channels['G'][loc],
            channels['T'][loc]
        ]
        
        sorted_signals = sorted(signals, reverse=True)
        primary_signal = sorted_signals[0]
        secondary_signal = sorted_signals[1]
        
        if primary_signal > 0:
            ratio = secondary_signal / primary_signal
            noise_ratios.append(ratio)
            
    #Sanger is very messy at the ends, so mess with this "trim" value a good bit
    trim = 10
    if len(noise_ratios) > 2 * trim:
        valid_ratios = noise_ratios[trim:-trim]
    else:
        valid_ratios = noise_ratios
        
    return np.mean(valid_ratios)

#Change FR to "F" or "R" depending on how data looks...
FR = "R"
file_first_mapping = {
    f"Kangas7925_G11-C_R1-16S-rRNA-seq{FR}.ab1": 0,
    f"Kangas7925_G11-1_R1-16S-rRNA-seq{FR}.ab1": 1,
    f"Kangas7925_G11-2_R1-16S-rRNA-seq{FR}.ab1": 2,
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

#Change these two variables if we want to analyze First trials or Second trials
file_mapping = file_second_mapping
base_path = "SecondSeq/ab1_files/"
results = []

print("Processing files and calculating mutation varience...\n")

for filename, exposure_time in file_mapping.items():
    full_path = os.path.join(base_path, filename)
    
    if not os.path.exists(full_path):
        print(f"Warning: Could not find {filename}")
        continue
        
    variance_score = analyze_variance(full_path)
    results.append({
        "File": filename,
        "UV_Time_Seconds": exposure_time,
        "Var_Score": variance_score
    })
    
    print(f"Processed {filename} (UV: {exposure_time}s) | Score: {variance_score:.4f}")

df = pd.DataFrame(results)

df = df.sort_values(by="UV_Time_Seconds")

plt.figure(figsize=(10, 6))
plt.plot(df['UV_Time_Seconds'], df['Var_Score'], marker='o', linestyle='-', color='b')
plt.title("DNA Mutation Variance vs. UV Light Exposure Time (Sample C8)")
plt.xlabel("UV Exposure Time (Seconds)")
plt.ylabel("Noise Ratio Average across all channels")
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()