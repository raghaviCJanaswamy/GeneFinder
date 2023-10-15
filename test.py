def getCoordinates(orf, DNA):
    """Get the beginning and end coordinates of an ORF in the DNA sequence."""
    start_index = DNA.find(orf)
    if start_index == -1:
        return None  # ORF not found in DNA
    end_index = start_index + len(orf)
    return (start_index, end_index)

# Example usage:
input_DNA = "ATGAAATGCTAAATGAAGGGTGA"
stopList = ["TAA", "TAG", "TGA"]
orf_list = findORFs(input_DNA)

for orf in orf_list:
    coordinates = getCoordinates(orf, input_DNA)
    if coordinates is not None:
        print(f"ORF: {orf}, Start: {coordinates[0]}, End: {coordinates[1]}")
