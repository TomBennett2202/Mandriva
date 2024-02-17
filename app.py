from flask import Flask,render_template, redirect, url_for, send_file

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/analysis')
def analysis():
    return render_template('analysis.html')

@app.route('/snp')
def snp():
    return render_template('snp.html')

@app.route('/cluster')
def cluster():
    return render_template('cluster.html')

@app.route('/admixture')
def admixture():
    return render_template('admixture.html')

@app.route('/result_analysis')
def result_analysis():
    return render_template('result_analysis.html')

@app.route('/result_snp')
def result_snp():
    return render_template('result_snp.html')

@app.route('/download')
def download_file():
    p = "PCA image.png"
    return send_file(p,as_attachment=True)

if __name__=="__main__":
    app.run(debug=True)
