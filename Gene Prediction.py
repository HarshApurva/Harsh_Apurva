from Bio import Entrez, SeqIO

# Function to fetch the sequence from NCBI
def fetch_sequence(accession):
    """Fetch sequence from NCBI using accession number."""
    Entrez.email = "har.apu.bt22@dypatil.edu"
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching data: {e}")
        return None

# Function to analyze promoter region
def predict_promoter(sequence):
    """Predict promoter region based on common motifs."""
    promoter_region = sequence[:100]
    if "TATA" in promoter_region:
        print("Promoter motif (TATA box) found!")
    else:
        print("No common promoter motif found.")
    return promoter_region

# Function to analyze terminator region
def predict_terminator(sequence):
    """Predict terminator region based on common motifs."""
    terminator_region = sequence[-100:]  
    if "AATAAA" in terminator_region:
        print("Terminator motif (poly-A signal) found!")
    else:
        print("No common terminator motif found.")
    return terminator_region

# Function to predict gene region using start and stop codons
def predict_gene(sequence):
    """Predict the gene region using start (ATG) and stop codons (TAA, TAG, TGA)."""
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]

    # Find the first occurrence of the start codon
    start_index = sequence.find(start_codon)
    if start_index == -1:
        print("No start codon (ATG) found.")
        return None


    # Find the first occurrence of any stop codon after the start codon
    gene_sequence = sequence[start_index:]
    stop_index = -1
    for stop_codon in stop_codons:
        stop_index = gene_sequence.find(stop_codon)
        if stop_index != -1:
            break

    if stop_index == -1:
        print("No stop codon found after the start codon.")
        return None

    # Extract the gene region
    gene_region = sequence[start_index:start_index + stop_index + 3]  # Including the stop codon
    print(f"Predicted gene region: {gene_region}")
    return gene_region

# Main function
def main():
    # Predefined accession number for non-interactive environments
    accession = "NM_001200001"  # Replace with the desired accession number
    print(f"Using accession number: {accession}")
    print("Fetching sequence...")
    record = fetch_sequence(accession)

    if record:
        print(f"Sequence fetched: {record.description}")
        sequence = str(record.seq)
        
        # Analyze promoter region
        print("\nAnalyzing promoter region...")
        promoter = predict_promoter(sequence)
        print(f"Promoter region: {promoter}\n")

        # Analyze terminator region
        print("Analyzing terminator region...")
        terminator = predict_terminator(sequence)
        print(f"Terminator region: {terminator}\n")

        # Predict gene region using start and stop codons
        print("Predicting gene region using start and stop codons...")
        gene_region = predict_gene(sequence)
        if gene_region:
            print(f"Gene region: {gene_region}\n")

        print("Gene analysis complete!")
    else:
        print("Failed to fetch sequence. Please check the accession number.")

if __name__ == "__main__":
    main()
