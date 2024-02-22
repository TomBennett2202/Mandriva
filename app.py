from flask import Flask,render_template, redirect, request, url_for, send_file
import plotly.express as px
from plotly.offline import plot
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import numpy as np
import base64
from io import BytesIO

DATABASE = r'/your/file/path/mandriva_database.db' #keep the r as it should read the db

def query_db(query, args=(), one=False):
    con = sqlite3.connect(DATABASE)
    cur = con.cursor()
    cur.execute(query, args)
    rv = cur.fetchall()
    con.close()
    return (rv[0] if rv else None) if one else rv


app = Flask(__name__)

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/analysis')
def analysis():
    return render_template('analysis.html')

@app.route('/snp', methods=['GET', 'POST'])
def snp():
    if request.method == 'POST':
        selected_populations = request.form.getlist('pop[]')# Ensure the form has a 'population' input
        return redirect(url_for('result_snp', populations=','.join(selected_populations)))
    else:
        return render_template('snp.html')

@app.route('/result_snp/<populations>')
def result_snp(populations):
    population_list = populations.split(',')  # Split the populations back into a list
    placeholders = ','.join('?' for _ in population_list)  # Create placeholders for query parameters
    return render_template('result_snp.html')


@app.route('/cluster', methods=['GET', 'POST'])
def cluster():
    if request.method == 'POST':
        selected_populations = request.form.getlist('pop[]')  # Ensure the form has a 'population' input
        return redirect(url_for('show_plot', populations=','.join(selected_populations)))
    else:
        return render_template('cluster.html')
    
@app.route('/show_plot/<populations>')
def show_plot(populations):
    population_list = populations.split(',')  # Split the populations back into a list
    placeholders = ','.join('?' for _ in population_list)  # Create placeholders for query parameters

    query = f'SELECT Population, PC1, PC2 FROM RESULTS WHERE Population IN ({placeholders})'
    with sqlite3.connect(DATABASE) as con:  # Use context manager
        df = pd.read_sql_query(query, con, params=population_list)
    
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
    'SIB': '#4f2f4f'   # purple taupe
    }

    fig = px.scatter(df, x='PC1', y='PC2', color='Population', title='PCA Plot for Selected Populations (PC1: 52% Variance Explained, PC2: 16% Variance Explained)', color_discrete_map=color_map)
    plot_div = plot(fig, output_type='div', include_plotlyjs=False)

    return render_template('plot.html', plot_div=plot_div)


@app.route('/admixture', methods=['GET', 'POST'])
def admixture():
    if request.method == 'POST':
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
        plt.close()  # Close the plot to prevent it from displaying in the Jupyter Notebook or Python script
        buf.seek(0)
        plot_url = base64.b64encode(buf.getvalue()).decode('utf8')
        buf.close()

        return render_template('admixture_results.html', plot_url=plot_url)
    else:
        return render_template('admixture.html')


@app.route('/admixture_results')
def admixture_results():
    # Your code to present the results goes here
    return render_template('admixture_results.html')

@app.route('/download')
def download_file():
    p = "PCA image.png"
    return send_file(p,as_attachment=True)

if __name__=="__main__":
    app.run(debug=True)
