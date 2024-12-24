import requests

def calculate_gc_content(dna_sequence):
    """
    Calculate the GC content of a given DNA sequence.

    :param dna_sequence: The DNA sequence as a string.
    :return: GC content as a percentage (float).
    """
    g_count = dna_sequence.count('G')
    c_count = dna_sequence.count('C')
    total_length = len(dna_sequence)
    gc_content = (g_count + c_count) / total_length
    return gc_content


def fetch_uniprot_data(size, query, label, make_fasta=None, make_fasta_name=None):
    """
    Fetch data from UniProt based on a query and optionally write results to a FASTA file.

    :param size: Number of hits to fetch.
    :param query: Search query combining organism and protein information.
    :param label: Label for the query to categorize results.
    :param make_fasta: Whether to create a FASTA file ('Y' or 'N').
    :param make_fasta_name: Name of the FASTA file if creating one.
    :return: A tuple of two lists: GC content and sequence lengths, and corresponding labels.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "size": size,
        "format": "json"
    }

    response = requests.get(url, params=params)

    if response.status_code != 200:
        print(f"Error: {response.status_code}")
        print(response.text)
        return [], []

    data = response.json()
    gc_content_seq_len = []  # List to accumulate GC content and sequence lengths
    labels = []

    if "results" in data:
        if make_fasta and make_fasta.upper() == "Y":
            with open(make_fasta_name, "w") as fasta_file:
                for entry in data["results"]:
                    accession = entry.get("primaryAccession", "N/A")
                    protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown")
                    organism = entry.get("organism", {}).get("scientificName", "Unknown")
                    sequence = entry.get("sequence", {}).get("value", "N/A")

                    fasta_file.write(f">{protein_name} | {organism} | Accession: {accession}\n")
                    for i in range(0, len(sequence), 80):
                        fasta_file.write(sequence[i:i+80] + "\n")
            print(f"FASTA file '{make_fasta_name}' has been created.")
        else:
            for entry in data["results"]:
                sequence = entry.get("sequence", {}).get("value", "N/A")
                length = entry.get("sequence", {}).get("length", "N/A")
                # Calculate GC content and store the pair (GC content, length)
                gc_content_seq_len.append([calculate_gc_content(sequence), length])
                labels.append(label)

    return gc_content_seq_len, labels


def fetch_data_for_multiple_queries(size, queries, make_fasta=None, make_fasta_name=None):
    """
    Fetch data for multiple queries, and return the results as a tuple of two lists: GC content/sequence lengths and labels.

    :param size: Number of hits to fetch for each query.
    :param queries: A list of queries (strings).
    :param make_fasta: Whether to create a FASTA file ('Y' or 'N').
    :param make_fasta_name: Name of the FASTA file if creating one.
    :return: A tuple of two lists: GC content and sequence lengths, and corresponding labels.
    """
    part_X = []
    part_Y = []

    for label, query in enumerate(queries):
        gc_content_seq_len, labels = fetch_uniprot_data(size, query, label, make_fasta, make_fasta_name)
        part_X.extend(gc_content_seq_len)
        part_Y.extend(labels)

    return part_X, part_Y



