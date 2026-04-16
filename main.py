from Bio import SeqIO



# Load the .ab1 file
filename = ""
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

# Example: Get signal strength for A at the 10th peak
# location = peak_locations[9]
# strength = channel_A[location]