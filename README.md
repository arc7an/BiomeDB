# BiomeDB

A Python module which uses functions from BioPython (ver. 1.63) to process molecular biology data and communicate with external databases. To administrate Neo4j (ver. 2.1.2) database py2neo (ver. 1.6.4) module was used. 

Links: [Neo4j](http://www.neo4j.org/download), [BioPython](http://biopython.org/wiki/Download), [py2neo](https://github.com/nigelsmall/py2neo)

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
## Adding new organism to the base
Let's add a new organism to the base. Before that you need to download a GenBank file form NCBI. For exapmle we are going to use [Enterobacteria phage T7 genome] (http://www.ncbi.nlm.nih.gov/nuccore/9627425).
````
graph_database.create_nodes('Enterobacteria_phage_T7.gb')
````

You may check now wheter organism is in the base:

````
graph_database.find_organism()
Escherichia_coli_str._K-12_substr._MG1655
Enterobacteria_phage_T7
[Node('http://localhost:7474/db/data/node/23489'),
 Node('http://localhost:7474/db/data/node/41232')]
````

Then you have to choose what organism to process further:

````
graph_database.set_organism('Enterobacteria_phage_T7')
Enterobacteria_phage_T7
````

## Ordering methods - NEXT and OVERLAP
Now you can use other methods to work with chosen organism. Let's make 'NEXT' and 'OVERLAP' relations between features (all elements of genome) of the T7:
````
graph_database.relation_next()
graph_database.relation_overlap()
````

If you want to set relation 'NEXT' and 'OVERLAP' between certain nodes with properties 'start' and 'end' you have to choose certain label:
````
graph_database.relation_next(genes)
graph_database.relation_overlap(DNA)
````

If you want to add a label to a node or a group of nodes you can use method 'add_label_to_nodes'. For example let's add for each 'Polypeptide' new label 'Enzyme':
````
graph_database.add_label_to_nodes('Polyppetide', 'Enzyme')
````
If you want to add a lavel for certain node you can use additional searching arguments. In this case you have set property key and value. In this example we are giong to add label 'Verified' to a node with label 'Term' and property 'name' with value 'DNA primase/helicase':
````
graph_database.add_label_to_nodes('Term', 'Verified', 'name', 'DNA primase/helicase')
````
## Editing nodes
There is also a method which allows to add new properties to nodes. Here you can add property and its value to a group of nodes with a certain label:
````
graph_database.add_property_to_nodes('Term', 'function', 'unknown')
````
You can also add property to a certain node, for example let's take the 'Term' node we used in previous method and change its name:
````
graph_database.add_property_to_nodes('Term', 'name', 'Example name', 'name', 'DNA primase/helicase')
````
## Homological proteins searching
There is a method that allow you to seek for homological proteins using both NCBI remote BLAST service and your local BLAST database if have install it. The arguments are offline (True/False) which define wether you want to you use remote NCBI BLAST or your local base, similar_quantity (int) which efine how many homological proteins you want to include to you base, e_value (float) - the BLAST expecte value parameter and db - a string containing path to your local BLAST database (if you choose online BLAST, set db = '' - empty string).
Example of online query:
````
graph_database.relation_similar(offline = False, similar_quantity = 3, e_value = 0.01, db='')
````
Example of offline query:
````
graph_database.relation_similar(offline = True, similar_quantity = 3, e_value = 0.01, db='/path/to/your/local/base')
````
I want to warn user that these processes are very slow - it takes 4-8 minutes to analyze one polypeptide. So you can break executing with CTRL+C and continue later. There is a veryfication for already analyzed polypeptides, so there will not be any repeats.
## Pattern seaching
There are two methods which allow to search simple patterns:

(x)-[r]->(y) and (a)-[r1]->(b)-[r2]->(c)

In the first case you need to use method. For example let's search for similar polypeptides:
````
out=graph_database.search_one_rel(label_out='Polypeptide', label_in='Polypeptide', relation='SIMILAR', return_number = -1)
````
"Out" contains list of tuples with  references (node, node, relation). You can fileter out put with return_number argument (default -1). When it equals -1, method returns the previous list, starting with 0 you can choose what element of the tuple you want to return. For example
````
out=graph_database.search_one_rel(label_out='Polypeptide', label_in='Polypeptide', relation='SIMILAR', return_number=2)
````
returns a list of relations only.

Now we want to find polypeptides that are a part of organism and we use the second method:
````
out=graph_database.search_two_rels('Polypeptide', 'Polypeptide', 'Organism', 'SIMILAR', 'PART_OF')
````
The output has the same structure as the previous one - it is a list of tuples. You can also filter output with return_number (default -1), it takes values from -1 to 4. For example let's get only homological polypeptides:
````
out=graph_database.search_two_rels(label_left='Polypeptide', label_middle='Polypeptide', label_right='Organism', relation_lm='SIMILAR', relation_mr='PART_OF', return_number=0)
````

Moreover these two methods allow to search using properties of the nodes, what can make queries more specific. All you need is just make input string argument as a Cypher query. For example we want to find a polypeptide which is encoded by a gene with a known locus_tag and that gene is a part of a certain organism:
````
out=graph_database.search_two_rels('Polypeptide', 'Gene{locus_tag:"T7p03"}', 'Organism{name:"Enterobacteria_phage_T7"}', 'ENCODES', 'PART_OF', 0)
````
The same structures can be used in search_one_rel method.
