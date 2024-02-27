import sqlite3
import cyvcf2
import re

# Open the VEP output file
with open('/Users/tombennett/ADMIXTURE/vep_output.txt', 'r') as file:
    # Read all lines
    lines = file.readlines()[33:]

# Initialize dictionaries to store ClinVar_CLNSIG and filtered ClinVar_CLNDN values
clinvar_clnsig_values = {}
filtered_clinvar_clndn_values = {}

# Process each line in the VEP output
for line in lines:
    # Check if the line contains ClinVar_CLNSIG and ClinVar_CLNDN
    if 'ClinVar_CLNSIG' in line:
        # Extract SNP ID and ClinVar_CLNSIG value using regex
        snp_id_match = re.search(r'^(\S+)', line)
        clnsig_match = re.search(r'ClinVar_CLNSIG=([^;]+)', line)

        # Append the values to the dictionaries if found
        if snp_id_match and clnsig_match:
            snp_id = snp_id_match.group(1)
            # Replace underscores with spaces and remove additional characters
            clnsig_value = clnsig_match.group(1).replace('_', ' ').split(',')[0]
            # Replace "|" with commas
            clnsig_value = clnsig_value.split('|')
            clinvar_clnsig_values[snp_id] = clnsig_value

    if 'ClinVar_CLNDN' in line:
        # Extract SNP ID and ClinVar_CLNDN value using regex
        snp_id_match = re.search(r'^(\S+)', line)
        clndn_match = re.search(r'ClinVar_CLNDN=([^;]+)', line)

        # Append the values to the dictionary if found
        if snp_id_match and clndn_match:
            snp_id = snp_id_match.group(1)
            clndn_value = clndn_match.group(1)
            # Replace underscores with spaces and remove additional characters
            clndn_value = clndn_value.replace('_', ' ').split(',')[0]
            # Replace "|" with commas
            clndn_value = clndn_value.split('|')

            filtered_clinvar_clndn_values[snp_id] = clndn_value

# Connect to the SQLite database (or create if not exists)
conn = sqlite3.connect('Mandriva/mandriva_database.db')
cursor = conn.cursor()

# Create a table for clinical significance
cursor.execute('''
    CREATE TABLE IF NOT EXISTS Clinical_Significance (
        ID INTEGER PRIMARY KEY AUTOINCREMENT,
        SNP_ID TEXT,
        Significance TEXT,
        FOREIGN KEY (SNP_ID) REFERENCES SNPs(SNP_ID)
    )
''')

# Insert data into the clinical significance table
for snp_id, significance_list in clinvar_clnsig_values.items():
    for significance in significance_list:
        cursor.execute('''
            INSERT INTO Clinical_Significance (SNP_ID, Significance) VALUES (?, ?)
        ''', (snp_id, significance))

# Create a table for disease names
cursor.execute('''
    CREATE TABLE IF NOT EXISTS Disease_Names (
        ID INTEGER PRIMARY KEY AUTOINCREMENT,
        SNP_ID TEXT,
        Disease_Name TEXT,
        FOREIGN KEY (SNP_ID) REFERENCES SNPs(SNP_ID)
    )
''')

# Insert data into the disease names table
for snp_id, disease_names_list in filtered_clinvar_clndn_values.items():
    for disease_name in disease_names_list:
        cursor.execute('''
            INSERT INTO Disease_Names (SNP_ID, Disease_Name) VALUES (?, ?)
        ''', (snp_id, disease_name))

# Commit changes
conn.commit()

###########################################################################################################

# Open your VCF file and skip the header lines
with open('/Users/tombennett/ADMIXTURE/variant_effect_output.txt', 'r') as file:
    lines = file.readlines()[15:]  # Skip the first 15 lines

snp_genes = {}

# Process all lines
for line_number, line in enumerate(lines, start=1):
    uploaded_variation, symbol = line.strip().split('\t')

    # Check if the SNP is already in the dictionary
    if uploaded_variation in snp_genes:
        # Add the gene to the SNP entry if it's not already present
        if symbol not in snp_genes[uploaded_variation]:
            snp_genes[uploaded_variation].append(symbol)
    else:
        # Initialize the SNP entry with the gene
        snp_genes[uploaded_variation] = [symbol]


# Create a single table with auto-incremented ID for GenesData
cursor.execute('''
    CREATE TABLE IF NOT EXISTS GenesData (
        ID INTEGER PRIMARY KEY AUTOINCREMENT,
        SNP_ID TEXT,
        Gene TEXT,
        FOREIGN KEY (SNP_ID) REFERENCES SNPs(SNP_ID)
    )
''')

# Insert data into GenesData table
for snp_id, genes_list in snp_genes.items():
    for gene in genes_list:
        cursor.execute('''
            INSERT INTO GenesData (SNP_ID, Gene) VALUES (?, ?)
        ''', (snp_id, gene))

# Commit changes
conn.commit()

###########################################################################################################

# Create SNPs table (if not exists)
cursor.execute('''
    CREATE TABLE IF NOT EXISTS SNPs (
        SNP_ID TEXT,
        RS_ID TEXT,
        Chromosome TEXT,
        Position INTEGER,
        Ref_Allele TEXT,
        Alt_Allele TEXT
    )
''')

# Commit the changes
conn.commit()

# Load VCF data into SNPs table
vcf_reader = cyvcf2.VCF('/Users/tombennett/ADMIXTURE/filtered_mandriva.vcf.gz')

second_vcf_reader = cyvcf2.VCF('/Users/tombennett/ADMIXTURE/variants_annotated.vcf.gz')

rsids = []
for record in second_vcf_reader:
    rsids.append(record.ID)

count = 0
for record in vcf_reader:
    # Insert data into the SQL table
    cursor.execute('''
        INSERT INTO SNPs (SNP_ID, RS_ID, Chromosome, Position, Ref_Allele, Alt_Allele)
        VALUES (?, ?, ?, ?, ?, ?)
    ''', (
        str(record.ID), str(rsids[count]), str(record.CHROM), int(record.POS),
        str(record.REF), str(record.ALT)))
    count += 1

# Commit the changes
conn.commit()

# Close the connection when done
conn.close()
