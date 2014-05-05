from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from py2neo import neo4j, node, rel
import operator

def create_nodes(data_base, seq_name):
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
    organism_nodes = list(data_base.find(organism))
    if len(organism_nodes) == 0:
        _create_features(data_base, rec, organism)
    else:
        print 'Organism is already in base.'
    return organism

def _cds2poly(data_base, feature, organism, element):
    """
    Additional function for create_features and create_nodes.
    Makes 'polypeptide' nodes for 'CDS' nodes and relation between them
    """
    polypeptide, = data_base.create(node({'seq': feature.qualifiers['translation'][0]}))
    polypeptide.add_labels('polypeptide', organism)
    xref, = data_base.create(node({'uniprot': feature.qualifiers['db_xref'][1].split(':')[1],
                                   'qualifier': 'is_a'}))
    xref.add_labels('XRef')
    data_base.create(rel(xref, 'REF', polypeptide))
    data_base.create(rel(polypeptide, 'ENCODED', element))

def _create_features(data_base, gb_record, organism):
    """
    Additional function for create_nodes.
    Makes nodes from input sequence elements with labels and properties.
    """
    for feature in gb_record.features:
        element, = data_base.create(node({'start': int(feature.location.start),
                                          'end': int(feature.location.end),
                                          'strand': feature.location.strand}))
        element.add_labels(str(feature.type), gb_record.description.split(',')[0])
        element.add_labels('feature', organism)
        if feature.type == 'CDS':
            _cds2poly(data_base, feature, organism, element)


def relation_next(data_base, organism):
    """
    The function creates relationships 'NEXT' in the 'data_base' between nodes with label 'label'.
    The nodes must have property 'start' as the function compare 2 'start' values
    to make the relationship.
    """
    #elements = data_base.find(organism)
    #must find only features
    #elements = list(elements)
    #must not start from 'source'
    #start = elements[1]['start']
    #for i in xrange(2,len(elements)):
    #    if 'feature' in elements[i]:
    #        if elements[i]['start'] >= start:
    #            data_base.create(rel(elements[i], 'NEXT', elements[i-1]))
    #            start = elements[i]['start']
    ordered = _feature_start_ordering(data_base, organism)
    for i in xrange(2, len(ordered)):
        data_base.create(rel(ordered[i][0], 'NEXT', ordered[i-1][0]))
    #checking for already NEXT



def _feature_start_ordering(data_base, organism):
    elements = list(data_base.find(organism))
    features = []
    for i in elements:
        if 'feature' in i.get_labels():
            features.append(i)
    del elements
    ordering = {}
    for feature in features:
        ordering[feature] = feature.get_properties()['start']
    return sorted(ordering.iteritems(), key=operator.itemgetter(1))




def relation_overlap(data_base, organism):
    """
    Function creates relationships 'OVERLAP' between features of organism with name 'organism_label'.
    If 'start' value of the node is located between 'start' and 'end' values of the previous nodes,
    the relationship is created.
    There is a check for property 'start'.
    """
    #elements = data_base.find(organism)
    #element_list = []
    #positions = []
    #for element in elements:
    #    #Make better find!
    #    if element['start'] != None:
    #        element_list.append(element)
    #        positions.append([element['start'], element['end']])
    #positions = positions[1:]
    #for i in xrange(1, len(positions)):
    #    for j in xrange(i-1, -1, -1):
    #        if positions[j][1] <= positions[i][1] and positions[j][1] >= positions[i][0]:
    #            data_base.create(rel(element_list[i], 'OVERLAP', element_list[j]))

    ordered = _feature_start_ordering(data_base, organism)
    for i in xrange(2, len(ordered)):
        for j in xrange(i-1, -1, -1):
            #<= or <?
            if ordered[j][0]['end'] <= ordered[i][0]['end'] and ordered[j][0]['end'] >= ordered[i][0]['start']:
                data_base.create(rel(ordered[i][0], 'OVERLAP', ordered[j][0]))


def add_label_to_nodes(data_base, label, label_to_add):
    """
    The function adds new 'label_to_add' to the element with known label 'label'.
    """
    elements = data_base.find(label)
    for element in elements:
        element.add_labels(label_to_add)

def relation_similar(data_base, organism, similar_quantity = 3):
    """
    The function searches for nodes with labels 'polypeptide', then functions searches among
    incoming relationships for one with 'SIMILAR' label. If it is not found than 'polypeptide'
    'seq' value is BLASTed with blastp function. After that output XML is read.
    For the first 'similar_quantity' homological proteins nodes are made with 'SIMILAR' relationships.
    For each homological new 'polypeptide' several nodes with relation 'LINK' are created with property
    'qualifier' and value 'is_a'. Links labels are 'pdb', 'gb', and/or 'sp'.
    """
    polypeptides = data_base.find('polypeptide')
    polypeptide_counter = 1
    for polypeptide in polypeptides:
        incoming_rels = polypeptide.match_incoming()
        already_exists = False
        for incoming in incoming_rels:
            if (incoming.type == 'SIMILAR') and (organism in polypeptide.get_labels()):
                already_exists = True
        if already_exists == False:
            protein_sequence, gb, pdb, sp = _blast_ncbi(str(polypeptide.get_properties()['seq']), similar_quantity)
            #blast_write = NCBIWWW.qblast('blastp', 'nr', str(polypeptide.get_properties()['seq']))
            #blast_read = NCBIXML.read(blast_write)
            print 'Blasted ' + str(polypeptide_counter)
            polypeptide_counter += 1
            for i in xrange(similar_quantity):
                sim_polypeptide, = data_base.create(node({'seq': protein_sequence}))
                sim_polypeptide.add_labels('polypeptide')
                #gb = blast_read.alignments[i].title.partition('gb|')[-1].partition('|')[0]
                #pdb = blast_read.alignments[i].title.partition('pdb|')[-1].partition('|')[0]
                #sp = blast_read.alignments[i].title.partition('sp|')[-1].partition('|')[0]
                data_base.create(rel(sim_polypeptide, 'SIMILAR', polypeptide))
                if len(gb[i]) > 0:
                    gb_node, = data_base.create(node({'name': gb[i]}))
                    gb_node.add_labels('gb')
                    data_base.create(rel(gb_node, ('LINK', {'qualifier':'is_a'}), sim_polypeptide))
                if len(pdb[i]) > 0:
                    pdb_node, = data_base.create(node({'name': pdb[i]}))
                    pdb_node.add_labels('pdb')
                    data_base.create(rel(pdb_node, ('LINK', {'qualifier':'is_a'}), sim_polypeptide))
                if len(sp[i]) > 0:
                    sp_node, = data_base.create(node({'name': sp[i]}))
                    sp_node.add_labels('sp')
                    data_base.create(rel(sp_node, ('LINK', {'qualifier':'is_a'}), sim_polypeptide))
                #print blast_read.alignments[i].title[:3] + blast_read.alignments[i].title[3:].split('|')[0] + '\n'
        else:
            print 'Found! ' + str(polypeptide_counter)
            polypeptide_counter += 1

def _blast_ncbi(protein_sequence, similar_quantity):
    blast_write = NCBIWWW.qblast('blastp', 'nr', protein_sequence)
    blast_read = NCBIXML.read(blast_write)
    gb, pdb, sp = [], [], []
    for i in xrange(similar_quantity):
        gb.append(str(blast_read.alignments[i].title.partition('gb|')[-1].partition('|')[0]))
        pdb.append(str(blast_read.alignments[i].title.partition('pdb|')[-1].partition('|')[0]))
        sp.append(str(blast_read.alignments[i].title.partition('sp|')[-1].partition('|')[0]))
    return blast_read.alignments[i].hsps[0].match, gb, pdb, sp
