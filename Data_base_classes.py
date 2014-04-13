from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from py2neo import neo4j, node, rel

def create_nodes(seq_name, data_base):
    """
    This function reads 'seq_name' in GenBank format, than it creates nodes with labels of organism
    name, sequence type and 'feature'.
    Each node contains 3 properties: 'start', 'end' and 'strand'.
    For each 'CDS'-node extra 'polypeptide'-node is created. 'CDS' and 'polypeptide' are linked by
    relationship with label 'ENCODED'. Each 'polypeptide'-node has property 'seq' which contains
    aminoacid sequence and label with organism name.
    For each 'polypeptide'-node 'XRef'-node is created. 'XRef'-nodes
    have properties 'uniprot', which contains the uniprot-link and 'qualifier' with value 'is_a'.
    'polypeptide' and 'XRef' are connected with relationship with label 'REF'
    """
    rec = SeqIO.read(seq_name, 'genbank')
    organism = rec.features[0].qualifiers['organism'][0]
    for feature in rec.features:
        element, = data_base.create(node({'start': int(feature.location.start),
                                          'end': int(feature.location.end),
                                          'strand': feature.location.strand}))
        element.add_labels(str(feature.type), rec.description.split(',')[0])
        element.add_labels('feature', organism)
        if feature.type == 'CDS':
            polypeptide, = data_base.create(node({'seq': feature.qualifiers['translation'][0]}))
            polypeptide.add_labels('polypeptide', organism)
            xref, = data_base.create(node({'uniprot': feature.qualifiers['db_xref'][1].split(':')[1],
                                           'qualifier': 'is_a'}))
            xref.add_labels('XRef')
            data_base.create(rel(xref, 'REF', polypeptide))
            data_base.create(rel(polypeptide, 'ENCODED', element))

def relation_next(data_base, label):
    """
    The function creates relationships 'NEXT' in the 'data_base' between nodes with label 'label'.
    The nodes must have property 'start' as the function compare 2 'start' values
    to make the relationship.
    """
    elements = data_base.find(label)
    elements = list(elements)
    #must not start from 'source'
    start = elements[1]['start']
    for i in xrange(2,len(elements)):
        if elements[i]['start'] >= start:
            data_base.create(rel(elements[i], 'NEXT', elements[i-1]))
            start = elements[i]['start']


def relation_overlap(data_base, organism_label):
    """
    Function creates relationships 'OVERLAP' between features of organism with name 'organism_label'.
    If 'start' value of the node is located between 'start' and 'end' values of the previous nodes,
    the relationship is created.
    There is a check for property 'start'.
    """
    elements = data_base.find(organism_label)
    element_list = []
    positions = []
    for element in elements:
        #Make better find!
        if element['start'] != None:
            element_list.append(element)
            positions.append([element['start'], element['end']])
    positions = positions[1:]
    for i in xrange(1, len(positions)):
        for j in xrange(i-1, -1, -1):
            if positions[j][1] <= positions[i][1] and positions[j][1] >= positions[i][0]:
                data_base.create(rel(element_list[i], 'OVERLAP', element_list[j]))


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
        incomings = polypeptide.match_incoming()
        already_exists = False
        for incoming in incomings:
            if (incoming.type == 'SIMILAR') and (organism in polypeptide.get_labels()):
                already_exists = True
        if already_exists == False:
            blast_write = NCBIWWW.qblast('blastp', 'nr', str(polypeptide.get_properties()['seq']))
            blast_read = NCBIXML.read(blast_write)
            print 'Blasted ' + str(polypeptide_counter)
            polypeptide_counter += 1
            for i in xrange(similar_quantity):
                sim_polypeptide, = data_base.create(node({'seq': blast_read.alignments[i].hsps[0].match}))
                sim_polypeptide.add_labels('polypeptide')
                gb = blast_read.alignments[i].title.partition('gb|')[-1].partition('|')[0]
                pdb = blast_read.alignments[i].title.partition('pdb|')[-1].partition('|')[0]
                sp = blast_read.alignments[i].title.partition('sp|')[-1].partition('|')[0]
                data_base.create(rel(sim_polypeptide, 'SIMILAR', polypeptide))
                if len(gb) > 0:
                    gb_node, = data_base.create(node({'name':str(gb)}))
                    gb_node.add_labels('gb')
                    data_base.create(rel(gb_node, ('LINK', {'qualifier':'is_a'}), sim_polypeptide))
                if len(pdb) > 0:
                    pdb_node, = data_base.create(node({'name':str(pdb)}))
                    pdb_node.add_labels('pdb')
                    data_base.create(rel(pdb_node, ('LINK', {'qualifier':'is_a'}), sim_polypeptide))
                if len(sp) > 0:
                    sp_node, = data_base.create(node({'name':str(sp)}))
                    sp_node.add_labels('sp')
                    data_base.create(rel(sp_node, ('LINK', {'qualifier':'is_a'}), sim_polypeptide))
                print blast_read.alignments[i].title[:3] + blast_read.alignments[i].title[3:].split('|')[0] + '\n'
        else:
            print 'Found! ' + str(polypeptide_counter)
            polypeptide_counter += 1