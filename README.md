# **Mandriva Website**
Mandriva is a genomic analysis web application that displays SNP variation and population genetics on human chromosome 1 at the continental level, providing insights into human ancestry and clinical information. It offers clustering, admixture analysis, and statistical summaries with interactive visual representations to explore genetic relationships and differentiation.

## Installing Mandriva
### Data
* Please download the required database (mandriva_database.db) from [Mandriva](https://qmulprod-my.sharepoint.com/personal/bt23629_qmul_ac_uk/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fbt23629%5Fqmul%5Fac%5Fuk%2FDocuments%2FMandriva&ct=1709042648099&or=OWA%2DNT&cid=bf1c1404%2D4efb%2D5f3f%2D8046%2D524f80580ab0&ga=1) (QMUL email is required to access the data)
* If you want to visualize the database, a DB Browser for SQLite is required, which can be downloaded from [here](https://sqlitebrowser.org/dl/) 


### Instructions 

1. `git clone https://github.com/TomBennett2202/Mandriva.git`
2. `cd Mandriva`
3. `pip3 install -r requirements.txt`



* Please ensure that the mandriva_database.db and app.py are in the same folder while running the app.
If you wish to change the database directory, you should also change the *DATABASE* path in the application.



