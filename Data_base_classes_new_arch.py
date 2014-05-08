from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from py2neo import neo4j, node, rel, cypher
from time import ctime
import operator

#rm -rf /var/lib/neo4j/data/graph.db/* - purges DB

class Biome_db():

    def __init__(self, link):
        self.data_base = neo4j.GraphDatabaseService(link)
        try:
            self.sp = list(self.data_base.find('sp'))[0]
        except:
            self.sp, = self.data_base.create(node())
            self.sp.add_labels('sp')
        try:
            self.gb = list(self.data_base.find('gb'))[0]
        except:
            self.gb, = self.data_base.create(node())
            self.gb.add_labels('gb')
        try:
            self.pdb = list(self.data_base.find('pdb'))[0]
        except:
            self.pdb, = self.data_base.create(node())
            self.pdb.add_labels('pdb')

    def create_nodes(self, seq_name):
        """
        This function reads 'seq_name' in GenBank format, than it creates nodes with labels of organism
        name, sequence type and 'feature'.
        Each node contains 3 properties: 'start', 'end' and 'strand'.
        For each 'CDS'-node extra 'polypeptide'-node is created. 'CDS' and 'polypeptide' are linked by
        relationship with label 'ENCODED'. Each 'polypeptide'-node has property 'seq' which contains
        aminoacid sequence and label with organism name.
        For each 'polypeptide'-node 'XRef'-node is created. 'XRef'-nodes
        have properties 'uniprot', which contains the uniprot-link and 'qualifier' with value 'is_a'.
        'polypeptide' and 'XRef' are connected with relationship with label 'REF'.
        """
        rec = SeqIO.read(seq_name, 'genbank')
        organism = rec.features[0].qualifiers['organism'][0]
        organism = organism.replace(' ', '_')
        print organism
        organism_nodes = list(self.data_base.find(organism))
        if len(organism_nodes) == 0:
            self._create_features(rec, organism)
        else:
            print 'Organism is already in base'
        return organism

    def _cds2poly(self, feature, element, source):
        """
        Additional function for create_features and create_nodes.
        Makes 'polypeptide' nodes for 'CDS' nodes and relation between them
        """
        polypeptide, = self.data_base.create(node({'seq': feature.qualifiers['translation'][0]}))
        polypeptide.add_labels('polypeptide')
        #xref, = self.data_base.create(node({'uniprot': feature.qualifiers['db_xref'][1].split(':')[1],
        #                               'qualifier': 'is_a'}))
        #
        #xref.add_labels('XRef')
        #self.data_base.create(rel(xref, 'REF', polypeptide))
        self.data_base.create(rel(polypeptide,
                              ('XREF', {'id': feature.qualifiers['db_xref'][1].split(':')[1]}), self.sp))
        self.data_base.create(rel(element, 'ENCODE', polypeptide))
        self.data_base.create(rel(polypeptide, 'PART_OF', source))

    def _create_features(self, gb_record, organism):
        """
        Additional function for create_nodes.
        Makes nodes from input sequence elements with labels and properties.
        """
        source, = self.data_base.create(node(({'start': int(gb_record.features[0].location.start),
                                              'end': int(gb_record.features[0].location.end),
                                              'strand': gb_record.features[0].location.strand})))
        source.add_labels('source', organism)
        #gb_record.features[0]
        for feature in gb_record.features[1:]:
            element, = self.data_base.create(node({'start': int(feature.location.start),
                                              'end': int(feature.location.end),
                                              'strand': feature.location.strand}))
            element.add_labels('feature', str(feature.type))
            self.data_base.create(rel((element, 'PART_OF', source)))
            if feature.type == 'CDS':
                self._cds2poly(feature, element, source)

    def relation_next(self, organism):
        """
        The function creates relationships 'NEXT' in the 'data_base' between nodes with label 'label'.
        The nodes must have property 'start' as the function compare 2 'start' values
        to make the relationship.
        """
        #check_list = self._relation_checker('feature', organism, 'NEXT')
        check_list = self.search_two_rels('feature', 'feature', organism, 'NEXT', 'PART_OF')
        if len(check_list) == 0:
            ordered = self._feature_start_ordering(organism)
            for i in xrange(2, len(ordered)):
                self.data_base.create(rel(ordered[i][0], 'NEXT', ordered[i-1][0]))
        else:
            print 'NEXT relationships is already exist for ' + organism

    def _feature_start_ordering(self, organism):
        query_out = self.search_one_rel('feature', organism, 'PART_OF')
        ordering = {}
        for out in query_out:
            feature = out[0]
            ordering[feature] = feature.get_properties()['start']
        return sorted(ordering.iteritems(), key=operator.itemgetter(1))

    def relation_overlap(self, organism):
        """
        Function creates relationships 'OVERLAP' between features of organism with name 'organism_label'.
        If 'start' value of the node is located between 'start' and 'end' values of the previous nodes,
        the relationship is created.
        There is a check for property 'start'.
        """
        check_list = self.search_two_rels('feature', 'feature', organism, 'OVERLAP', 'PART_OF')
        if len(check_list) == 0:
            ordered = self._feature_start_ordering(organism)
            for i in xrange(2, len(ordered)):
                for j in xrange(i-1, -1, -1):
                    #<= or <?
                    if ordered[j][0]['end'] <= ordered[i][0]['end'] and ordered[j][0]['end'] >= ordered[i][0]['start']:
                        self.data_base.create(rel(ordered[i][0], 'OVERLAP', ordered[j][0]))
        else:
            print 'OVERLAP relationships is already exist for ' + organism


    def add_label_to_nodes(self, label, label_to_add):
        """
        The function adds new 'label_to_add' to the element with label 'label'.
        """
        elements = self.data_base.find(label)
        for element in elements:
            element.add_labels(label_to_add)

    def add_property_to_nodes(self, label, property_to_add, value_to_add):
        """
        Function adds new property 'property_to_add' with value 'value_to_add' to the element with label 'label'.
        """
        elements = self.data_base.find(label)
        for element in elements:
            properties = element.get_properties()
            properties[property_to_add] = value_to_add
            element.update_properties(properties)

    def gene_analyzer(self, organism):
        #make proper find nodes
        genes = list(self.data_base.find('gene'))
        cdss = list(self.data_base.find('CDS'))
        mrnas = list(self.data_base.find('mRNA'))
        for gene in genes:
            start = gene.get_properties()['start']
            end = gene.get_properties()['end']
            for cds in cdss:
                cds_start = cds.get_properties()['start']
                cds_end = cds.get_properties()['end']
                if start <= cds_start <= end and start <= cds_end <= end:
                    self.data_base.create(rel(cds, 'VERSION_OF', gene))
            for mrna in mrnas:
                mrna_start = mrna.get_properties()['start']
                mrna_end = mrna.get_properties()['end']
                if start <= mrna_start <= end and start <= mrna_end <= end:
                    self.data_base.create(rel(gene, 'ENCODE', mrna))

    def relation_similar(self, organism, offline = True, similar_quantity = 3, e_value = 0.01):
        """
        The function searches for nodes with labels 'polypeptide', then functions searches among
        incoming relationships for one with 'SIMILAR' label. If it is not found than 'polypeptide'
        'seq' value is BLASTed with blastp function. After that output XML is read.
        For the first 'similar_quantity' homological proteins nodes are made with 'SIMILAR' relationships.
        For each homological new 'polypeptide' several nodes with relation 'LINK' are created with property
        'qualifier' and value 'is_a'. Links labels are 'pdb', 'gb', and/or 'sp'.
        """
        polypeptides_all = self.data_base.find('polypeptide')
        #polypeptides_analyzed = list(self._search_similar_polypeptides(organism))
        polypeptides_analyzed = []
        blasted_searching = self.search_two_rels('polypeptide', 'polypeptide', organism, 'SIMILAR', 'PART_OF')
        for out in blasted_searching:
            polypeptides_analyzed.append(out[1])
        polypeptide_counter = 1
        for polypeptide in polypeptides_all:
            already_exists = False
            if polypeptide in polypeptides_analyzed:
                already_exists = True
            if already_exists == False:
                protein_txt = open('current_protein.txt', 'w')
                protein_txt.write(str(polypeptide.get_properties()['seq']))
                print ctime()
                protein_txt.close()
                print str(polypeptide.get_properties()['seq'])
                homolog_sequence, gb, pdb, sp = \
                    self._blaster('current_protein.txt', offline=True, similar_quantity=similar_quantity, e_value=e_value)
                print 'Blasted:' + str(polypeptide_counter)
                polypeptide_counter += 1
                for i in xrange(similar_quantity):
                    sim_polypeptide, = self.data_base.create(node({'seq': homolog_sequence}))
                    sim_polypeptide.add_labels('polypeptide')
                    self.data_base.create(rel(sim_polypeptide, 'SIMILAR', polypeptide))
                    if len(gb[i]) > 0:
                        self.data_base.create(rel(sim_polypeptide, ('XREF', {'id': gb[i]}), self.gb))
                    if len(pdb[i]) > 0:
                        self.data_base.create(rel(sim_polypeptide, ('XREF', {'id': pdb[i]}), self.pdb))
                    if len(sp[i]) > 0:
                        self.data_base.create(rel(sim_polypeptide, ('XREF', {'id': sp[i]}), self.sp))
            else:
                print 'Found! ' + str(polypeptide_counter)
                polypeptide_counter += 1

    def _search_similar_polypeptides(self, organism):
        session = cypher.Session()
        transaction = session.create_transaction()
        #query = 'MATCH (a:' + organism + ')-[r:PART_OF]-(b:' + label + ') RETURN b'
        query = 'match (poly_new:polypeptide)-[:`SIMILAR`]->(poly_org:polypeptide)-[:`PART_OF`]->(organism:' + organism + ') return poly_org'
        transaction.append(query)
        #print query
        transaction_out = transaction.execute()
        nodes = []
        for node in transaction_out[0]:
            #print node
            nodes.append(node[0])
        return nodes

    def _blaster(self, protein_sequence, offline = True, similar_quantity = 3, e_value = 0.01):
        if offline:
            cline = NcbiblastpCommandline(query = protein_sequence,
                                db = '/home/artem/BLAST_DB/nr', evalue = e_value,
                                outfmt = 5, out = 'prot_out.xml')
            cline()
            blast_write = open('prot_out.xml')
            blast_read = NCBIXML.read(blast_write)
        else:
            #add e_value for online
            blast_write = NCBIWWW.qblast('blastp', 'nr', protein_sequence)
            blast_read = NCBIXML.read(blast_write)
        homolog_seq, gb, pdb, sp = [], [], [], []
        #check alignments - they can be empty!
        for i in xrange(similar_quantity):
            gb.append(str(blast_read.alignments[i].title.partition('gb|')[-1].partition('|')[0]))
            pdb.append(str(blast_read.alignments[i].title.partition('pdb|')[-1].partition('|')[0]))
            sp.append(str(blast_read.alignments[i].title.partition('sp|')[-1].partition('|')[0]))
            homolog_seq.append(blast_read.alignments[i].hsps[0].match)
        #return blast_read.alignments[i].hsps[0].match, gb, pdb, sp
        return homolog_seq, gb, pdb, sp

    def search_one_rel(self, label_out, label_in, relation):
        session = cypher.Session()
        transaction = session.create_transaction()
        query = 'MATCH (node_out:' + label_out + ')' \
                '-[rel:`' + relation + '`]->' \
                ' (node_in:' + label_in + ') RETURN node_out, node_in, rel'
        transaction.append(query)
        transaction_out = transaction.execute()
        out = []
        for result in transaction_out[0]:
            out.append(result.values)
        return out

    def search_two_rels(self, label_left, label_middle, label_right, relation_lm, relation_mr):
        session = cypher.Session()
        transaction = session.create_transaction()
        query = 'MATCH (label_left:' + label_left + ')-[relation_lm:`' + relation_lm + '`]->(label_middle:' + label_middle + ')-[relation_mr:`' + relation_mr + '`]->(label_right:' + label_right + ') RETURN label_left, label_middle, label_right, relation_lm, relation_mr'
        transaction.append(query)
        transaction_out = transaction.execute()
        out = []
        for result in transaction_out[0]:
            out.append(result.values)
        return out

#Testing function
#def blast_standalone(protein_sequence = '/home/artem/work/reps/neo4j/prot.txt', e_value = 0.01, similar_quantity = 3):
#    cline=NcbiblastpCommandline(query=protein_sequence,
#                                db='/home/artem/BLAST_DB/nr', evalue=e_value,
#                                outfmt=5, out='prot_out.xml')
#    cline()
#    res = open('prot_out.xml')
#    rec = NCBIXML.read(res)
    #for record in rec.alignments[:similar_quantity]:
#seems to work. Check on more protein examples and make as class

