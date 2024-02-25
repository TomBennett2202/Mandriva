import sqlite3
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
from flask import Flask, render_template, request, send_file
from plotly.offline import plot
matplotlib.use('agg')
import base64
from io import BytesIO
import numpy as np
import seaborn as sns

DATABASE = r'mandriva_database.db' #this is the file path to read the database

app = Flask(__name__)

# Pre-made functions
def get_snp_ids_for_genes(gene_names):
    connection = sqlite3.connect(DATABASE)
    cursor = connection.cursor()
    query = """
        SELECT DISTINCT GenesData.SNP_ID
        FROM GenesData
        JOIN SNPs ON GenesData.SNP_ID = SNPs.SNP_ID
        WHERE GenesData.Gene IN ({});
    """.format(','.join(['?'] * len(gene_names)))

    cursor.execute(query, gene_names)
    result = cursor.fetchall()
    connection.close()
    # Flatten the list of tuples into a single list of snp ids
    unique_snp_ids = [row[0] for row in result]
    # Return the list of SNP_IDs
    return unique_snp_ids


def get_snp_ids_for_rs_ids(rs_ids):
    connection = sqlite3.connect(
        DATABASE)
    cursor = connection.cursor()
    query = "SELECT SNP_ID FROM SNPs WHERE RS_ID IN ({})"
    query = query.format(','.join(['?'] * len(rs_ids)))
    cursor.execute(query, rs_ids)
    result = cursor.fetchall()
    connection.close()
    # Flatten the list of tuples into a single list of snp ids
    snp_ids = [row[0] for row in result]
    # Return the list of SNP IDs
    return snp_ids


def fetch_allele_frequency_from_database(snp_ids, populations):
    conn = sqlite3.connect(DATABASE)
    cursor = conn.cursor()

    query = f"SELECT SNP_ID, {', '.join(populations)} FROM Allele_Count WHERE SNP_ID IN ({', '.join(['?']*len(snp_ids))})"
    cursor.execute(query, snp_ids)

    results = cursor.fetchall()

    conn.close()

    allele_frequencies = {}

    for result in results:
        snp_id, *allele_counts = result

        allele_counts = list(map(lambda x: list(map(int, x.split(','))), allele_counts))
        total_counts = [sum(counts) for counts in zip(*allele_counts)]

        frequencies = [count / sum(total_counts) if sum(total_counts) > 0 else 0 for count in total_counts]

        allele_frequencies[snp_id] = {'freqs': frequencies}

    return allele_frequencies


def fetch_genotype_frequency_from_database(snp_ids, populations):
    conn = sqlite3.connect(DATABASE)
    cursor = conn.cursor()

    query = f"SELECT SNP_ID, {', '.join(populations)} FROM Genotype_Count WHERE SNP_ID IN ({', '.join(['?']*len(snp_ids))})"
    cursor.execute(query, snp_ids)

    results = cursor.fetchall()

    conn.close()

    genotype_frequencies = {}

    for result in results:
        snp_id, *genotype_counts = result

        genotype_counts = list(map(lambda x: list(map(int, x.split(','))), genotype_counts))
        total_counts = [sum(counts) for counts in zip(*genotype_counts)]

        frequencies = [count / sum(total_counts) if sum(total_counts) > 0 else 0 for count in total_counts]

        genotype_frequencies[snp_id] = {'freqs': frequencies}

    return genotype_frequencies


def get_rs_ids_for_snp_ids(snp_ids):
    conn = sqlite3.connect(DATABASE)
    cursor = conn.cursor()

    # Create a placeholder for the IN clause in the SQL query
    placeholders = ', '.join(['?' for _ in snp_ids])

    # Query to retrieve RS IDs for specified SNP IDs from the SNPs table
    query = f"SELECT SNP_ID, RS_ID FROM SNPs WHERE SNP_ID IN ({placeholders});"

    cursor.execute(query, snp_ids)
    results = cursor.fetchall()

    # Create a dictionary to store the mapping of SNP_ID to RS_ID
    rs_ids_for_snp_ids = {snp_id: rs_id for snp_id, rs_id in results}

    # Close the database connection
    conn.close()

    return rs_ids_for_snp_ids


def get_snps_in_region(chromosome, start_position, end_position):
    # Connect to the SQLite database
    conn = sqlite3.connect(DATABASE)
    cursor = conn.cursor()
    # Query for SNP IDs within the specified region
    query = """
    SELECT SNP_ID
    FROM SNPs
    WHERE Chromosome = ? AND Position BETWEEN ? AND ?
    """
    cursor.execute(query, (chromosome, start_position, end_position))

    snp_ids = [result[0] for result in cursor.fetchall()]

    conn.close()

    return snp_ids


def get_unique_significance_and_diseases_for_snps(snp_ids):
    conn = sqlite3.connect(DATABASE)
    cursor = conn.cursor()

    # Initialize dictionary to store unique significances and disease names for each SNP
    results_dict = {}

    for snp_id in snp_ids:
        # Fetch unique clinical significance for the current SNP
        cursor.execute("""
            SELECT DISTINCT Significance FROM Clinical_Significance
            WHERE SNP_ID = ?
        """, (snp_id,))
        clinical_significances = [row[0] for row in cursor.fetchall()]

        # Fetch unique disease names for the current SNP
        cursor.execute("""
            SELECT DISTINCT Disease_Name FROM Disease_Names
            WHERE SNP_ID = ?
        """, (snp_id,))
        disease_names = [row[0] for row in cursor.fetchall()]

        # Combine results into a tuple and add to the dictionary
        results_dict[snp_id] = (clinical_significances, disease_names)

    conn.close()

    return results_dict


# Function to calculate FST for a single SNP
def calculate_fst(hom_ref, het, hom_alt, ref_count, alt_count):
    # Calculate the total count of alleles for each SNP
    total_count = [sum(x) for x in zip(hom_ref, het, hom_alt)]

    # Calculate the frequency of the reference allele (p) for each SNP
    p = [(2 * hr + h) / (2 * tc) if tc > 0 else 0 for hr, h, tc in zip(hom_ref, het, total_count)]

    # Calculate the frequency of the alternate allele (q) for each SNP
    q = [1 - pi for pi in p]

    # Calculate the FST value for each SNP
    fst = [1 - (pi ** 2 + qi ** 2) if tc > 0 else 0 for pi, qi, tc in zip(p, q, total_count)]

    return fst

def calculate_fst_for_selection(populations, snps):
    # Connect to the SQLite database
    conn = sqlite3.connect(DATABASE)
    # Query genotype and allele counts for selected populations and SNPs
    genotype_counts_query = f"SELECT SNP_ID, {', '.join(populations)} FROM Genotype_Count WHERE SNP_ID IN ({', '.join(['?' for _ in snps])})"
    allele_counts_query = f"SELECT SNP_ID, {', '.join(populations)} FROM Allele_Count WHERE SNP_ID IN ({', '.join(['?' for _ in snps])})"

    # Retrieve genotype and allele count data as DataFrames
    genotype_counts_df = pd.read_sql_query(genotype_counts_query, conn, params=snps)
    allele_counts_df = pd.read_sql_query(allele_counts_query, conn, params=snps)

    # Create an empty DataFrame to store FST values
    fst_matrix = pd.DataFrame(index=genotype_counts_df['SNP_ID'], columns=populations)

    # Iterate through each population and calculate FST for each SNP
    for population in populations:
        # Extract genotypic counts for each SNP and population
        hom_ref, het, hom_alt = zip(*genotype_counts_df[population].str.split(',').apply(
            lambda x: (int(x[0]), int(x[1]), int(x[2])) if x[0] != '' else (0, 0, 0)))
        ref_count, alt_count = zip(*allele_counts_df[population].str.split(',').apply(
            lambda x: (int(x[0]), int(x[1])) if x[0] != '' else (0, 0)))

        # Calculate FST for each SNP using the previously defined function
        fst_values = calculate_fst(hom_ref, het, hom_alt, ref_count, alt_count)

        # Store FST values in the matrix
        fst_matrix[population] = fst_values

    conn.close()

    return fst_matrix


# Function to calculate pairwise FST between populations
def calculate_pairwise_fst(fst_matrix):
    # Get the list of populations from the columns of the FST matrix
    populations = fst_matrix.columns

    # Create an empty DataFrame to store pairwise FST values
    pairwise_fst_matrix = pd.DataFrame(index=populations, columns=populations)

    # Iterate through each pair of populations and calculate pairwise FST
    for pop1 in populations:
        for pop2 in populations:
            # Ensure that we're comparing different populations
            if pop1 != pop2:
                # Calculate pairwise FST as the mean absolute difference between population FST values
                pairwise_fst = np.mean(np.abs(fst_matrix[pop1] - fst_matrix[pop2]))

                # Store the pairwise FST value in the matrix
                pairwise_fst_matrix.loc[pop1, pop2] = pairwise_fst
                pairwise_fst_matrix.loc[pop2, pop1] = pairwise_fst

    # Return the matrix of pairwise FST values
    return pairwise_fst_matrix


@app.route('/')
def home():
    return render_template('home.html')

@app.route('/analysis')
def analysis():
    return render_template('analysis.html')

@app.route('/snp', methods=['GET', 'POST'])
def snp():
    if request.method != 'POST':
        return render_template('snp.html')

    else:
        try:
            # Make an empty list of snp ids to store for later
            snp_ids = []

            # Extract all user inputs
            selected_populations = request.form.getlist('pop[]')
            gene_name = request.form.get('gene_name')
            rs_id = request.form.get('snp_id')
            chromosome = request.form.get('chromosome')
            start_position = request.form.get('start')
            end_position = request.form.get('end')

            # If gene names are provided find unique snp ids for them
            if gene_name is not None and gene_name.strip() != '':
                gene_list = gene_name.split(',')
                # Add to snp ids list
                snp_ids.extend(get_snp_ids_for_genes(gene_list))


            # If rs ids are provided find snp ids for them
            if rs_id is not None and rs_id.strip() != '':
                rs_id_list = rs_id.split(',')
                # Add to snp ids list
                snp_ids.extend(get_snp_ids_for_rs_ids(rs_id_list))

            # Check if all three text boxes have values and get snps for the location
            if chromosome and start_position and end_position:
                # Add to snp ids list
                snp_ids.extend(get_snps_in_region(chromosome, start_position, end_position))

            # Get allele frequencies for snp ids
            allele_frequencies = fetch_allele_frequency_from_database(snp_ids, selected_populations)
            # Get genotype frequencies for snp ids
            genotype_frequencies = fetch_genotype_frequency_from_database(snp_ids, selected_populations)
            # Get rs ids to combine into the table
            rs_ids = get_rs_ids_for_snp_ids(snp_ids)

            clinical_data = get_unique_significance_and_diseases_for_snps(snp_ids)

            # Merge the data
            combined_data = {}
            for snp, allele_data in allele_frequencies.items():
                rs_id = rs_ids.get(snp, '')  # Get RS ID for the current SNP ID
                clinical_significances, disease_names = clinical_data.get(snp, ([], []))

                combined_data[snp] = {
                    'rs_id': rs_id,
                    'allele_freqs': allele_data['freqs'],
                    'clinical_significances': clinical_significances,
                    'disease_names': disease_names
                }

                if snp in genotype_frequencies:
                    combined_data[snp]['genotype_freqs'] = genotype_frequencies[snp]['freqs']

            # Calculate FST only if 2 or more populations are selected
            if len(selected_populations) >= 2:
                # Calculate FST for the selected populations and SNPs
                result_matrix = calculate_fst_for_selection(selected_populations, snp_ids)

                # Check if the result_matrix has data
                if not result_matrix.empty:
                    # Calculate pairwise FST matrix
                    pairwise_fst_matrix = calculate_pairwise_fst(result_matrix)

                    # Save pairwise FST matrix to a text file
                    output_file_path = "pairwise_fst_matrix.txt"
                    pairwise_fst_matrix.to_csv(output_file_path, sep='\t')

                    # Visual representation of the matrix using a heatmap
                    plt.figure(figsize=(10, 8))
                    sns.heatmap(pairwise_fst_matrix.astype(float), annot=True, cmap="coolwarm", linewidths=.5,
                                annot_kws={"fontsize": 6})
                    plt.title("Pairwise FST Matrix")
                    plt.tight_layout()

                    # Save the heatmap as an image
                    img_buf = BytesIO()
                    plt.savefig(img_buf, format='png')
                    img_buf.seek(0)
                    img_data = base64.b64encode(img_buf.getvalue()).decode('utf-8')
                    plt.close()

                    # Pass the image data to the template
                    return render_template('result_snp.html', combined_data=combined_data, heatmap_data=img_data)
                else:
                    # Handle the case where result_matrix is empty
                    return render_template('result_snp.html', combined_data=combined_data, heatmap_data=None)
            else:
                # Handle the case where less than 2 populations are selected
                return render_template('result_snp.html', combined_data=combined_data, heatmap_data=None)
        except (ValueError, sqlite3.OperationalError):
            # Handle the case where no data is available
            error_message = "No data found matching the provided input. Please try again with different values."
            return render_template('error_page.html', error_message=error_message)

@app.route('/cluster', methods=['GET', 'POST'])
def cluster():
    if request.method != 'POST':
        return render_template('cluster.html')

    else:
        selected_populations = request.form.getlist('pop[]')  # Ensure the form has a 'population' input
        placeholders = ','.join('?' for population in selected_populations)
        query = f'SELECT Population, PC1, PC2 FROM RESULTS WHERE Population IN ({placeholders})'
        with sqlite3.connect(DATABASE) as con:  # Use context manager
            df = pd.read_sql_query(query, con, params=selected_populations)

        color_map = {
            'ACB': '#1f77b4',  # muted blue
            'ASW': '#ff7f0e',  # safety orange
            'ESN': '#2ca02c',  # cooked asparagus green
            'GWD': '#d62728',  # brick red
            'LWK': '#9467bd',  # muted purple
            'MSL': '#8c564b',  # chestnut brown
            'YRI': '#e377c2',  # raspberry yogurt pink
            'CLM': '#7f7f7f',  # middle gray
            'MXL': '#bcbd22',  # curry yellow-green
            'PEL': '#17becf',  # blue-teal
            'PUR': '#1fa77b',  # jade green
            'CDX': '#b41f77',  # dark magenta
            'CHB': '#77b41f',  # olive green
            'CHS': '#1f77a5',  # steel blue
            'JPT': '#a51f77',  # dark pink
            'KHV': '#77a51f',  # avocado green
            'CEU': '#a5771f',  # bronze
            'FIN': '#7fa51f',  # fern green
            'GBR': '#1f7fa5',  # peacock blue
            'IBS': '#a51f7f',  # cerise
            'TSI': '#5f5f5f',  # dim gray
            'BEB': '#2f4f4f',  # dark slate gray
            'GIH': '#ff69b4',  # hot pink
            'ITU': '#b4ff69',  # lime
            'PJL': '#69b4ff',  # sky blue
            'STU': '#ff69b4',  # light pink
            'SIB': '#4f2f4f'  # purple taupe
        }

        fig = px.scatter(df, x='PC1', y='PC2', color='Population',
                         color_discrete_map=color_map)
        plot_div = plot(fig, output_type='div', include_plotlyjs=False)

        return render_template('plot.html', plot_div=plot_div)


@app.route('/admixture', methods=['GET', 'POST'])
def admixture():
    if request.method != 'POST':
        return render_template('admixture.html')

    else:
        selected_populations = request.form.getlist('pop[]')
        placeholders = ','.join('?' * len(selected_populations))
        query = f'''
            SELECT K1, K2, K3, K4, K5
            FROM RESULTS 
            WHERE POPULATION IN ({placeholders});
        '''

        # Connect to the SQLite database
        conn = sqlite3.connect(DATABASE)
        cursor = conn.cursor()

        # Execute the query and fetch the results into a pandas DataFrame
        cursor.execute(query, selected_populations)
        df = pd.DataFrame(cursor.fetchall(), columns=['K1', 'K2', 'K3', 'K4', 'K5'])

        # Close the cursor and connection
        cursor.close()
        conn.close()

        # Prepare the data for plotting
        k_values = df.values

        # Plotting
        plt.figure(figsize=(10, 6))
        colors = ['r', 'g', 'b', 'c', 'm']  # Colors for K1 to K5

        # Create a stacked bar chart
        bar_positions = np.arange(len(k_values))
        for i in range(k_values.shape[1]):  # Iterate over K columns
            bottom = np.sum(k_values[:, :i], axis=1) if i > 0 else np.zeros(len(k_values))
            plt.bar(bar_positions, k_values[:, i], bottom=bottom, color=colors[i], label=f'K{i+1}')

        # Customize the plot
        plt.xlabel('Population selected')
        plt.ylabel('Proportion')
        plt.title('Probability Distribution for K1 to K5')
        plt.legend()

        # Hide x-ticks
        plt.xticks([])

        # Save the plot to a BytesIO buffer
        buf = BytesIO()
        plt.savefig(buf, format='png')
        plt.close()
        buf.seek(0)
        plot_url = base64.b64encode(buf.getvalue()).decode('utf8')
        buf.close()

        return render_template('admixture_results.html', plot_url=plot_url)


# Add a new route for downloading the pairwise FST matrix
@app.route('/download_fst_matrix')
def download_fst_matrix():
    # Set the file path to your pairwise_fst_matrix.txt
    file_path = "pairwise_fst_matrix.txt"
    # Provide the file for download
    return send_file(file_path, as_attachment=True, download_name='pairwise_fst_matrix.txt', mimetype='text/plain')


if __name__=="__main__":
    app.run(debug=True)
