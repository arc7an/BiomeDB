# BiomeDB

A Python module which uses functions from BioPython to process molecular biology data and communicate with external databases. To administrate Neo4j database py2neo module was used.

## Quick start
In Python shell:
````
from Data_base_classes_new_arch import *
````

This imports class, contained in module and additiinal modules such BioPython, py2neo, etc.
After that create an object to create and work with the base. The arguments are the link to the base and optional logger seting.

````
graph_database = Biome_db('http://localhost:7474/db/data/', logger_level = logging.INFO)
````

After that the initializator tries to connect to DB. Try to find whether there is any organism in the DB:

````
graph_database.find_organism()
````

The output is the list of strings with names of organisms:
````
Escherichia_coli_str._K-12_substr._MG1655
Out: [Node('http://localhost:7474/db/data/node/23489')]
````

If there is no organisms the list will be empty.

Let's add a new organism to the base. Before that you need to download a GenBank file form NCBI. For exapmle we are going to use Enterobacteria phage T7 genome (http://www.ncbi.nlm.nih.gov/nuccore/9627425).
````
graph_database.create_nodes('Enterobacteria_phage_T7.gb')
````

