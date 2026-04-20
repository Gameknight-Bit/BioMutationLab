import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO


def analyze_peak_height(filepath):
    # Use BioPython to read the data from the ab1 files
    record = SeqIO.read(filepath, "abi")
    trace_data = record.annotations['abif_raw']

    # FWO_1 contains a string like b'GATC' showing order of DATA9-DATA12
    fwo = trace_data['FWO_1'].decode('utf-8')
    channels = {}
    for i, base in enumerate(fwo):
        channels[base] = trace_data[f'DATA{9+i}']

    peak_locations = trace_data['PLOC1']

    peak_heights = []

    for loc in peak_locations:
        signals = [
            channels['A'][loc],
            channels['C'][loc],
            channels['G'][loc],
            channels['T'][loc]
        ]
        primary = max(signals)
        peak_heights.append(primary)

    # Sanger is very messy at the ends, trim same as Jace did when analyzing variance 
    trim = 10
    if len(peak_heights) > 2 * trim:
        valid_heights = peak_heights[trim:-trim]
    else:
        valid_heights = peak_heights

    return np.mean(valid_heights)


# Change FR to "F" or "R" depending on how data looks
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

# Change these two variables to switch between First and Second trials
file_mapping = file_second_mapping
base_path = "SecondSeq/ab1_files/"

results = []

print("Processing files and calculating average peak heights...\n")

for filename, exposure_time in file_mapping.items():
    full_path = os.path.join(base_path, filename)

    if not os.path.exists(full_path):
        print(f"Warning: Could not find {filename}")
        continue

    avg_height = analyze_peak_height(full_path)
    results.append({
        "File": filename,
        "UV_Time_Seconds": exposure_time,
        "Avg_Peak_Height": avg_height
    })

    print(f"Processed {filename} (UV: {exposure_time}s) | Avg Peak Height: {avg_height:.2f}")

df = pd.DataFrame(results)
df = df.sort_values(by="UV_Time_Seconds")

plt.figure(figsize=(10, 6))
plt.plot(df['UV_Time_Seconds'], df['Avg_Peak_Height'], marker='o', linestyle='-', color='b')
plt.title("Average Primary Peak Height vs. UV Light Exposure Time")
plt.xlabel("UV Exposure Time (Seconds)")
plt.ylabel("Average Primary Peak Height (RFU)")
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()
