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
        self.sp = list(self.data_base.find('sp'))
        self.gb = list(self.data_base.find('gb'))
        self.pdb = list(self.data_base.find('pdb'))
        if len(self.sp) == 0:
            self.sp, = self.data_base.create(node())
            self.sp.add_labels('sp')
        if len(self.gb) == 0:
            self.gb, = self.data_base.create(node())
            self.gb.add_labels('gb')
        if len(self.pdb) == 0:
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
            print 'Organism is already in base.'
        return organism

    def _cds2poly(self, feature, element, source):
        """
        Additional function for create_features and create_nodes.
        Makes 'polypeptide' nodes for 'CDS' nodes and relation between them
        """
        polypeptide, = self.data_base.create(node({'seq': feature.qualifiers['translation'][0]}))
        polypeptide.add_labels('polypeptide')
        xref, = self.data_base.create(node({'uniprot': feature.qualifiers['db_xref'][1].split(':')[1],
                                       'qualifier': 'is_a'}))
        xref.add_labels('XRef')
        self.data_base.create(rel(xref, 'REF', polypeptide))
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
        gb_record.features[0]
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
        check_list = self._relation_checker('feature', organism, 'NEXT')
        if len(check_list) == 0:
            ordered = self._feature_start_ordering(organism)
            for i in xrange(2, len(ordered)):
                self.data_base.create(rel(ordered[i][0], 'NEXT', ordered[i-1][0]))
        else:
            print 'NEXT relationships is already exist for ' + organism
        #checking for already NEXT

    def _relation_checker(self, search_label, organism, relation):
        session = cypher.Session()
        transaction = session.create_transaction()
        query = 'MATCH (f1:`' + search_label + '`)-[:`' + relation + '`]->(f2:`'\
                  + search_label + '`)-[:`PART_OF`]->(org:`' + organism + '`) RETURN f2'
        transaction.append(query)
        transaction_out = transaction.execute()
        return transaction_out[0]

    def _feature_start_ordering(self, organism):
        #session = cypher.Session()
        #transaction = session.create_transaction()
        #query = 'MATCH (a:' + organism + ')-[r:PART_OF]-(b:feature) RETURN b'
        #transaction.append(query)
        #transaction_out = transaction.execute()
        #features = []
        #for element in transaction_out[0]:
        #    features.append(element[0])
        ##giva a to _feature_start_ordering
        #del transaction_out
        features = self._search_parts_of_organism(organism, 'feature')
        ordering = {}
        for feature in features:
            ordering[feature] = feature.get_properties()['start']
        return sorted(ordering.iteritems(), key=operator.itemgetter(1))

    def _search_parts_of_organism(self, organism, label):
        session = cypher.Session()
        transaction = session.create_transaction()
        query = 'MATCH (a:' + organism + ')-[r:PART_OF]-(b:' + label + ') RETURN b'
        transaction.append(query)
        transaction_out = transaction.execute()
        nodes = []
        for node in transaction_out[0]:
            #print node
            nodes.append(node[0])
        del transaction_out
        return nodes

    #def _feature_start_ordering(self, organism):
    #    elements = list(self.data_base.find(organism))
    #    features = []
    #    for i in elements:
    #        if 'feature' in i.get_labels():
    #            features.append(i)
    #    del elements
    #    #
    #    ordering = {}
    #    for feature in features:
    #        ordering[feature] = feature.get_properties()['start']
    #    return sorted(ordering.iteritems(), key=operator.itemgetter(1))

    def relation_overlap(self, organism):
        """
        Function creates relationships 'OVERLAP' between features of organism with name 'organism_label'.
        If 'start' value of the node is located between 'start' and 'end' values of the previous nodes,
        the relationship is created.
        There is a check for property 'start'.
        """
        check_list = self._relation_checker('feature', organism, 'OVERLAP')
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
        The function adds new 'label_to_add' to the element with known label 'label'.
        """
        elements = self.data_base.find(label)
        for element in elements:
            element.add_labels(label_to_add)

    def add_property_to_nodes(self, label, property_to_add, value_to_add):
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
            #print gene
            start = gene.get_properties()['start']
            end = gene.get_properties()['end']
            #print start, end
            for cds in cdss:
                cds_start = cds.get_properties()['start']
                cds_end = cds.get_properties()['end']
                #print cds_start, cds_end
                if start <= cds_start <= end and start <= cds_end <= end:
                    self.data_base.create(rel(cds, 'VERSION_OF', gene))
                    #print 'VERSION_OF'
            for mrna in mrnas:
                mrna_start = mrna.get_properties()['start']
                mrna_end = mrna.get_properties()['end']
                if start <= mrna_start <= end and start <= mrna_end <= end:
                    self.data_base.create(rel(gene, 'ENCODE', mrna))
                    #print 'ENCODE'

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
        polypeptides_analyzed = list(self._search_similar_polypeptides(organism))
        polypeptide_counter = 1
        #for polypeptide in polypeptides_all:
        #    incoming_rels = polypeptide.match_incoming()
        #    already_exists = False
        #    for incoming in incoming_rels:
        #        #cypher query match (c:polypeptide)-[:`SIMILAR`]->(a:polypeptide)-[:`PART_OF`]->(b:Enterobacteria_phage_T7) return a
        #        if (incoming.type == 'SIMILAR') and (organism in polypeptide.get_labels()):
        #            already_exists = True
        for polypeptide in polypeptides_all:
            already_exists = False
            if polypeptide in polypeptides_analyzed:
                already_exists =True
            if already_exists == False:
                protein_txt = open('current_protein.txt', 'w')
                protein_txt.write(str(polypeptide.get_properties()['seq']))
                print ctime()
                protein_txt.close()
                print str(polypeptide.get_properties()['seq'])
                protein_sequence, gb, pdb, sp = \
                    self._blaster('current_protein.txt', offline=True, similar_quantity=similar_quantity, e_value=e_value)
                #blast_write = NCBIWWW.qblast('blastp', 'nr', str(polypeptide.get_properties()['seq']))
                #blast_read = NCBIXML.read(blast_write)
                print 'Blasted ' + str(polypeptide_counter)
                polypeptide_counter += 1
                for i in xrange(similar_quantity):
                    sim_polypeptide, = self.data_base.create(node({'seq': protein_sequence}))
                    sim_polypeptide.add_labels('polypeptide')
                    self.data_base.create(rel(sim_polypeptide, 'SIMILAR', polypeptide))
                    if len(gb[i]) > 0:
                        gb_node, = self.data_base.create(node({'name': gb[i]}))
                        gb_node.add_labels('gb')
                        self.data_base.create(rel(gb_node, ('LINK', {'qualifier':'is_a'}), sim_polypeptide))
                    if len(pdb[i]) > 0:
                        pdb_node, = self.data_base.create(node({'name': pdb[i]}))
                        pdb_node.add_labels('pdb')
                        self.data_base.create(rel(pdb_node, ('LINK', {'qualifier':'is_a'}), sim_polypeptide))
                    if len(sp[i]) > 0:
                        sp_node, = self.data_base.create(node({'name': sp[i]}))
                        sp_node.add_labels('sp')
                        self.data_base.create(rel(sp_node, ('LINK', {'qualifier':'is_a'}), sim_polypeptide))
                    #print blast_read.alignments[i].title[:3] + blast_read.alignments[i].title[3:].split('|')[0] + '\n'
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
        del transaction_out
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
        gb, pdb, sp = [], [], []
        for i in xrange(similar_quantity):
            gb.append(str(blast_read.alignments[i].title.partition('gb|')[-1].partition('|')[0]))
            pdb.append(str(blast_read.alignments[i].title.partition('pdb|')[-1].partition('|')[0]))
            sp.append(str(blast_read.alignments[i].title.partition('sp|')[-1].partition('|')[0]))
        return blast_read.alignments[i].hsps[0].match, gb, pdb, sp

    #def _blast_ncbi(self, protein_sequence, similar_quantity):
    #    blast_write = NCBIWWW.qblast('blastp', 'nr', protein_sequence)
    #    blast_read = NCBIXML.read(blast_write)
    #    gb, pdb, sp = [], [], []
    #    for i in xrange(similar_quantity):
    #        gb.append(str(blast_read.alignments[i].title.partition('gb|')[-1].partition('|')[0]))
    #        pdb.append(str(blast_read.alignments[i].title.partition('pdb|')[-1].partition('|')[0]))
    #        sp.append(str(blast_read.alignments[i].title.partition('sp|')[-1].partition('|')[0]))
    #    return blast_read.alignments[i].hsps[0].match, gb, pdb, sp
    #
    #def _blast(self, protein_sequence, offline = True, similar_quantity = 3, e_value = 0.01):
    #    if offline:
    #        cline=NcbiblastpCommandline(query=protein_sequence,
    #                            db='/home/artem/BLAST_DB/nr', evalue=e_value,
    #                            outfmt=5, out='prot_out.xml')
    #        blast_write = open('prot_out.xml')
    #        blast_read = NCBIXML.read(blast_write)
    #    else:
    #        #add e_value for online
    #        blast_write = NCBIWWW.qblast('blastp', 'nr', protein_sequence)
    #        blast_read = NCBIXML.read(blast_write)
    #    gb, pdb, sp = [], [], []
    #    for i in xrange(similar_quantity):
    #        gb.append(str(blast_read.alignments[i].title.partition('gb|')[-1].partition('|')[0]))
    #        pdb.append(str(blast_read.alignments[i].title.partition('pdb|')[-1].partition('|')[0]))
    #        sp.append(str(blast_read.alignments[i].title.partition('sp|')[-1].partition('|')[0]))
    #    return blast_read.alignments[i].hsps[0].match, gb, pdb, sp

def blast_standalone(protein_sequence = '/home/artem/work/reps/neo4j/prot.txt', e_value = 0.01, similar_quantity = 3):
    cline=NcbiblastpCommandline(query=protein_sequence,
                                db='/home/artem/BLAST_DB/nr', evalue=e_value,
                                outfmt=5, out='prot_out.xml')
    cline()
    res = open('prot_out.xml')
    rec = NCBIXML.read(res)
    #for record in rec.alignments[:similar_quantity]:
#seems to work. Check on more protein examples and make as class

