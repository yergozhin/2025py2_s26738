#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever
Extended script to connect to NCBI, retrieve genetic sequence records for a given taxonomic ID,
filter by sequence length, and generate a CSV report and visualization.
"""

from Bio import Entrez, SeqIO
import time
import os
import pandas as pd
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key

        # Set up Entrez
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid):
        """Search for all records associated with a taxonomic ID."""
        print(f"Searching for records with taxID: {taxid}")
        try:
            # Get taxonomy information first
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            # Search for records
            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y", retmax=10000)
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name}")
                return None

            print(f"Found {count} records")

            # Store search results for later use
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count
            self.taxid = taxid
            self.organism = organism_name

            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=500):
        """Fetch a batch of records using the stored search results."""
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            # Limit to prevent overloading
            batch_size = min(max_records, 500)

            # Pause to avoid hitting rate limits
            time.sleep(0.4)

            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )

            # Parse GenBank records
            records = list(SeqIO.parse(handle, "genbank"))
            print(f"Fetched {len(records)} records.")
            return records

        except Exception as e:
            print(f"Error fetching records: {e}")
            return []

    def filter_records(self, records, min_len, max_len):
        """Filter records by sequence length."""
        return [
            (rec.id, len(rec.seq), rec.description)
            for rec in records if min_len <= len(rec.seq) <= max_len
        ]

    def generate_csv(self, filtered_data):
        """Generate a CSV report with accession, length, and description."""
        df = pd.DataFrame(filtered_data, columns=["Accession", "Length", "Description"])
        csv_path = os.path.join(os.getcwd(), f"taxid_{self.taxid}_report.csv")
        df.to_csv(csv_path, index=False)
        print(f"Saved CSV report to {csv_path}")
        return df

    def plot_length_chart(self, df):
        """Generate a line chart showing sequence lengths."""
        df_sorted = df.sort_values("Length", ascending=False)
        plt.figure(figsize=(10, 5))
        plt.plot(df_sorted["Accession"], df_sorted["Length"], marker='o')
        plt.xticks([], [])  # Hide x-axis labels to reduce clutter
        plt.ylabel("Sequence Length")
        plt.title(f"Sequence Lengths for TaxID {self.taxid} ({self.organism})")
        plt.tight_layout()
        plot_path = os.path.join(os.getcwd(), f"taxid_{self.taxid}_plot.png")
        plt.savefig(plot_path)
        print(f"Saved plot to {plot_path}")

def main():
    # Get user credentials
    while True:
        email = input("Enter your email address for NCBI: ").strip()
        if "@" in email and "." in email:
            break
        print("Invalid email format. Please enter a valid email.")

    api_key = input("Enter your NCBI API key (or leave blank if none): ").strip()

    # Create retriever object
    retriever = NCBIRetriever(email, api_key)

    # Get taxid and validate
    while True:
        taxid = input("Enter taxonomic ID (taxid) of the organism: ").strip()
        if taxid.isdigit():
            break
        print("TaxID must be a number. Try again.")

    # Get and validate sequence length inputs
    while True:
        try:
            min_len = int(input("Enter minimum sequence length: "))
            max_len = int(input("Enter maximum sequence length: "))
            if min_len <= 0 or max_len <= 0:
                print("Lengths must be positive integers.")
            elif min_len > max_len:
                print("Minimum length cannot be greater than maximum length.")
            else:
                break
        except ValueError:
            print("Please enter valid integer values for lengths.")

    # Search for records
    count = retriever.search_taxid(taxid)

    if not count:
        print("No records found. Exiting.")
        return

    # Fetch records
    records = retriever.fetch_records(start=0, max_records=500)

    # Filter by length
    filtered = retriever.filter_records(records, min_len, max_len)
    if not filtered:
        print("No records matched the length criteria.")
        return

    # Generate CSV and plot
    df = retriever.generate_csv(filtered)
    retriever.plot_length_chart(df)

if __name__ == "__main__":
    main()
