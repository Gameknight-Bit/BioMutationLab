from Bio import SeqIO

FirstPath = "FirstSeq"

# Load the .ab1 file
filename = FirstPath + "/ab1_files/Kangas7925_G11-1_R1-16S-rRNA-seqF.ab1"
record = SeqIO.read(filename, "abi")

# Access raw trace data
# abif_raw contains the raw binary data
trace_data = record.annotations['abif_raw']

# DATA9-12 typically contain the raw signals for A, C, G, T (order can vary)
# Use FWO_ to determine which channel corresponds to which base
channel_A = trace_data['DATA9']
channel_C = trace_data['DATA10']
channel_G = trace_data['DATA11']
channel_T = trace_data['DATA12']

# Access peak locations (PLOC)
peak_locations = trace_data['PLOC1']
print("Peak locations:", peak_locations)

# So if we can copy this code and do it for all "valid" location peaks across all files, then we can
# compare the signal strength differences across locations (and the most different from our control gives
# use the most likely mutation location)
location = 4
strength = channel_A[peak_locations[location]]
print(f"Signal strength for A at {location}th peak:", strength)
strength = channel_C[peak_locations[location]]
print(f"Signal strength for C at {location}th peak:", strength)
strength = channel_G[peak_locations[location]]
print(f"Signal strength for G at {location}th peak:", strength)
strength = channel_T[peak_locations[location]]
print(f"Signal strength for T at {location}th peak:", strength)