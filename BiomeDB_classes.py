from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from py2neo import neo4j, node, rel, cypher
from time import ctime, time
import datetime
import warnings
import logging



class Biome_db():

    def __init__(self, link = 'http://localhost:7474/db/data/', logger_level = logging.INFO):
        logging.basicConfig(filename = 'BiomeDB.log', level = logger_level,
                            format = '%(asctime)s %(message)s - %(module)s.%(funcName)s', datefmt='%H:%M:%S-%d.%m.%y')
        self.logger = logging.getLogger(__name__)
        self.logger.info('Initialization of database')
        self.data_base = neo4j.GraphDatabaseService(link)
        self.date = ''
        self.current_organism = None
        self.current_organism_ref = None
        self.current_organism_name = None
        self.start_time = 0
        self.stop_time = 0
        self.cds_product_counter = 0
        self.circular_flag = False
        try:
            self._base_init()
            self.logger.info('Data bases nodes (SP, GB, PDB) are found/created.')
        except:
            self.logger.error('Could not connect to the base. The link is not correct or the server does not answer.')
            raise ValueError('The link is not correct or the server does not answer.')

    def _base_init(self):
        """
        Private method valuecheck suppose to be done by caller method. Never invoke directly!
        """
        try:
            self.sp = list(self.data_base.find('DB', 'name', 'SP'))[0]
            self.logger.info('Node DB with {name: SP} is already in base.')
        except:
            self.sp, = self.data_base.create(node({'name': 'SP'}))
            self.sp.add_labels('DB')
            self.logger.info('Node DB with {name: SP} created.')
        try:
            self.gb = list(self.data_base.find('DB', 'name', 'GB'))[0]
            self.logger.info('Node DB with {name: GB} is already in base.')
        except:
            self.gb, = self.data_base.create(node({'name': 'GB'}))
            self.gb.add_labels('DB')
            self.logger.info('Node DB with {name: GB} created.')
        try:
            self.pdb = list(self.data_base.find('DB', 'name', 'PDB'))[0]
            self.logger.info('Node DB with {name: PDB} is already in base.')
        except:
            self.pdb, = self.data_base.create(node({'name': 'PDB'}))
            self.pdb.add_labels('DB')
            self.logger.info('Node DB with {name: PDB} created.')

    def set_organism(self, organism):
        """
        Methods sets current organism user is going to work with. When a new organism is created it becomes current by default.
        """
        if not isinstance(organism, basestring):
            self.logger.error('TypeError: Argument must be a string.')
            raise TypeError('Argument must be a string.')

        organism = self._normalize(organism)
        self.current_organism_name = organism
        query_out = list(self.data_base.find('Organism', 'name', organism))
        if len(query_out) == 0:
            err_message = 'There is no node %s in the database or it does not have an "organism" label.' % organism
            self.logger.error(err_message)
            raise ValueError(err_message)

        if len(query_out) > 1:
            err_message = 'There are several nodes with name %s in the database.' % organism
            self.logger.error(err_message)
            raise Warning(err_message)

        self.current_organism = query_out[0]
        self.current_organism_ref = str(self.current_organism).split(' ')[0].split('(')[1]
        if self.current_organism.get_properties()['circular'] == 'True':
            self.circular_flag = True
        else:
            self.circular_flag = False
        log_message = '%s is chosen by user, link in the base: %s' % (self.current_organism_name, self.current_organism)
        self.logger.info(log_message)
        print self.current_organism_name

    def find_organism(self):
        """
        Method returns list of organisms contained in the base.
        """
        organism_list = list(self.data_base.find('Organism'))

        if len(organism_list) == 0:
            warn_message = 'There is no organism in the data base.'
            self.logger.warning(warn_message)
            warnings.warn(warn_message)

        for organism in organism_list:
            self.logger.info('Found organism(s) in the base.')
            print organism['name']
        return organism_list

    def _search_organism_one_rel(self, label, return_number = -1):
        if not isinstance(self.current_organism_name, basestring):
            err_message = 'The organism was not chosen. Set current organism before using methods.'
            self.logger.error(err_message)
            raise UserWarning(err_message)

        organism = 'Organism{name:"%s"}' % self.current_organism_name
        return self.search_one_rel(label, organism, 'PART_OF', return_number)

    def _search_organism_two_rels(self, first_node, second_node, relation, return_number = -1):
        if not isinstance(self.current_organism_name, basestring):
            err_message = 'The organism was not chosen. Set current organism before using methods.'
            self.logger.error(err_message)
            raise UserWarning(err_message)

        organism = 'Organism{name:"%s"}' % self.current_organism_name
        return self.search_two_rels(first_node, second_node, organism, relation, 'PART_OF', return_number)

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
        self.start_time = time()
        if not isinstance(seq_name, basestring):
            err_message = 'Argument must be a string.'
            self.logger.error(err_message)
            raise TypeError(err_message)

        try:
            rec = SeqIO.read(seq_name, 'genbank')
        except:
            err_message = 'The input file must be of a genbank format (.gb).'
            self.logger.error(err_message)
            raise ValueError(err_message)

        seq_file = open(seq_name, 'r')
        first_line = seq_file.readline()
        seq_file.close()
        if 'circular' in first_line:
            self.circular_flag = True
        organism = rec.features[0].qualifiers['organism'][0]
        organism = self._normalize(organism)
        print organism
        self._create_features(rec, organism)
        self.stop_time = time()
        print self.stop_time - self.start_time

    def _normalize(self, input_string):
        """
        Private method valuecheck suppose to be done by caller method. Never invoke directly!
        """
        if ' ' in input_string:
            return input_string.replace(' ', '_')
        else:
            return input_string

    def _cds2poly(self, cds, gene, gb_record):
        """
        Private method valuecheck suppose to be done by caller method. Never invoke directly!

        Additional function for create_features and create_nodes.
        Makes 'polypeptide' nodes for 'CDS' nodes and relation between them
        """
        try:
            polypeptide, = self.data_base.create(node({'seq': cds.qualifiers['translation'][0]}))
            id = self._find_sp_link(cds.qualifiers['db_xref'])
            xref = self._create_xref(id)
            self.data_base.create(rel(polypeptide, 'EVIDENCE', xref))
            self.data_base.create(rel(xref, 'LINK', self.sp))
        except:
            translate = str(gb_record.seq[cds.location.nofuzzy_start:cds.location.nofuzzy_end])
            translate = Seq(translate, generic_dna)
            translate = translate.translate()[:-1]
            polypeptide, = self.data_base.create(node({'seq': str(translate)}))

        polypeptide.add_labels('Polypeptide', 'BioEntity')
        self.data_base.create(rel(gene, 'ENCODES', polypeptide))
        self.data_base.create(rel(polypeptide, 'PART_OF', self.current_organism))
        try:
            self._create_term(cds.qualifiers['product'][0], polypeptide)
            self.cds_product_counter += 1
        except:
            pass

    def _find_sp_link(self, db_xref_list):
        for db_xref in db_xref_list:
            name, value = db_xref.split(':')
            if name == 'UniProtKB/Swiss-Prot':
                return value

    def _create_term(self, product, polypeptide):
        term, = self.data_base.create(node({'name': product}))
        term.add_labels('Term')
        self.data_base.create(rel(polypeptide, ('HAS_NAME', {'Source': 'GenBank'}), term))

    def _create_xref(self, id):
        """
        Private method valuecheck suppose to be done by caller method. Never invoke directly!
        """
        xref, = self.data_base.create(node({'id': id}))
        xref.add_labels('Reference')
        return xref

    def _create_features(self, gb_record, organism):
        """
        Additional function for create_nodes.
        Makes nodes from input sequence elements with labels and properties.
        !! NO UPDATE.
        """
        organism_node = list(self.data_base.find('Organism', 'name', organism))
        if len(organism_node) == 0:
            genome, = self.data_base.create(node(({'length': int(gb_record.features[0].location.end),
                                                  'name': organism,
                                                  'circular': str(self.circular_flag)})))
            genome.add_labels('Genome', 'Organism')
            log_message = 'Created node for organism %s' % organism
            self.logger.info(log_message)
        self.set_organism(organism)
        print 'Current organism: %s' % self.current_organism_name
        cdss, genes = [], []
        features = [feature for feature in gb_record.features if feature.type != 'mRNA']
        #upd check for created features
        print 'Creating features'
        i = 1
        for feature in features[1:]:
            if feature.type == 'CDS':
                cdss.append(feature)
            elif feature.type == 'gene':
                genes.append(feature)
            elif 'RNA' in feature.type or 'rna' in feature.type:
                if feature.type == 'rRNA':
                    feature.type = 'tRNA'
                elif feature.type == 'ncRNA':
                    feature.type = 'sRNA'
                self._rna_dna_maker(feature, 'RNA')
            else:
                self._rna_dna_maker(feature, 'DNA')
            i += 1
        log_message = '%d elements in record, %d CDSs and %d genes.' % (i, len(cdss), len(genes))
        self.logger.info(log_message)
        print log_message
        if len(cdss) != 0 or len(genes) != 0:
            print 'Creating genes! There are %d cdss and %d genes' % (len(cdss), len(genes))
            self._search_or_create_genes(genes, cdss, gb_record)

    def _rna_dna_maker(self, feature, feature_type):
        strand = self._strand2str(feature)
        feature_node, = self.data_base.create(node({'start': int(feature.location.start),
                                            'end': int(feature.location.end),
                                            'strand': strand}))
        if feature_type == 'DNA':
            feature_node.add_labels('Feature', str(feature.type).capitalize(), feature_type, 'BioEntity')
        elif feature_type == 'RNA':
            feature_node.add_labels('Feature', str(feature.type), feature_type, 'BioEntity')
        self.data_base.create(rel((feature_node, 'PART_OF', self.current_organism)))

    def _strand2str(self, feature):
        if feature.location.strand == 1:
            return 'forward'
        elif feature.location.strand == -1:
            return 'reverse'

    def _gene_cds_match(self, genes, cdss):
        pairs = []
        genes_copy = genes[:]
        for gene in genes:
            locus_tag = gene.qualifiers['locus_tag']
            for cds in cdss:
                if cds.qualifiers['locus_tag'] == locus_tag:
                    pairs.append([gene, cds])
                    cdss.pop(cdss.index(cds))
                    genes_copy.pop(genes_copy.index(gene))
        if len(cdss) != 0:
            warn_mess = 'CDSs without genes: %d!' % len(cdss)
            self.logger.warning(warn_mess)
            warnings.warn(warn_mess)
        log_message = 'pairs=%d, nocds genes=%d' % (len(pairs), len(genes_copy))
        self.logger.info(log_message)
        print log_message
        return pairs + genes_copy

    def _search_or_create_genes(self, genes, cdss, gb_record):
        pairs = self._gene_cds_match(genes, cdss)
        print 'Found gene-cds pairs.'
        session = cypher.Session()
        transaction = session.create_transaction()
        query = 'MATCH (g:Gene)--(o:Organism)' \
                ' WHERE o.name = "%s"' \
                ' RETURN g, g.start, g.end, g.strand' \
                % self.current_organism_name
        transaction.append(query)
        transaction_out = transaction.execute()
        genes_db_props = [result.values[1:] for result in transaction_out[0]]
        genes_db = [result.values[0] for result in transaction_out[0]]
        updated_genes = 0
        created_genes = 0
        self.cds_product_counter = 0
        for item in pairs:
            if isinstance(item, list):
                strand = self._strand2str(item[0])
                start = int(item[0].location.nofuzzy_start + 1)
                end = int(item[0].location.nofuzzy_end)
                locus_tag = item[0].qualifiers['locus_tag'][0]
                search_gene = (start, end, strand)
                if search_gene in genes_db_props:
                    #found_index = genes_db_props.index(search_gene)
                    gene = genes_db[genes_db_props.index(search_gene)]
                    gene.update_properties({'locus_tag': locus_tag})
                    updated_genes += 1
                else:
                    gene, = self.data_base.create(node({'start': start,
                                                        'end': end,
                                                        'strand': strand,
                                                        'locus_tag': locus_tag}))
                    gene.add_labels('Gene', 'DNA', 'BioEntity', 'Feature')
                    self.data_base.create(rel((gene, 'PART_OF', self.current_organism)))
                    created_genes += 1
                self._cds2poly(item[1], gene, gb_record)
            elif item.type == 'gene':
                print 'NoCDS gene!'
                gene, = self.data_base.create(node({'match (o:Genome), (f:Feature) create (o)<-[:PART_OF]-(f)start': item.location.nofuzzy_start + 1,
                                                    'end': item.location.nofuzzy_end,
                                                    'strand': self._strand2str(item),
                                                    'locus_tag': item.qualifiers['locus_tag'][0]}))
                gene.add_labels('Gene', 'DNA', 'BioEntity', 'Feature')
                created_genes += 1
            else:
                err_message ='Pairs must contain list of genbank records or a genbank record.'
                self.logger.error(err_message)
                raise TypeError(err_message)
        log_message = '%d genes updated, %d genes created, %d CDSs with product' % (updated_genes, created_genes, self.cds_product_counter)
        self.logger.info(log_message)

    def relation_next(self, type_of_node = 'Feature'):
        """
        The function creates relationships 'NEXT' in the 'data_base' between nodes with label 'type_of_node'.
        The nodes must have property 'start' as the function compare 2 'start' values
        to make the relationship.
        >>> db = Biome_db()
        >>> db.relation_next(8)
        Traceback (most recent call last):
            ...
            TypeError: Argument must be a string.
        >>> db.relation_next('8')
        Traceback (most recent call last):
        ...
            ValueError: 'There is no node 8 in the database or it does not have an "organism" label.'
        """
        self.start_time = time()
        if not isinstance(self.current_organism_name, basestring):
            err_message = 'The organism was not chosen. Set current organism before using methods.'
            self.logger.error(err_message)
            raise UserWarning(err_message)

        if len(list(self.data_base.find('Organism', 'name', self.current_organism_name))) == 0:
            err_message = 'There is no node %s in the database or it does not have an "Organism" label.' % self.current_organism_name
            self.logger.error(err_message)
            raise ValueError(err_message)

        for strand in ('forward', 'reverse'):
            features = self._feature_start_ordering(strand, type_of_node)
            for i in xrange(1, len(features)):
                self.data_base.create(rel(features[i-1], 'NEXT', features[i]))
            if self.circular_flag == True:
                self.data_base.create(rel(features[-1], 'NEXT', features[0]))
            log_message = 'Created %d relations for %s in %s strand in %s' % (len(features), type_of_node, strand, self.current_organism_name)
            self.logger.info(log_message)
        self.stop_time = time()
        print self.stop_time - self.start_time

    def _feature_start_ordering(self, strand, type_of_node):
        """
        Private method valuecheck suppose to be done by caller method. Never invoke directly!
        """
        session = cypher.Session()
        transaction = session.create_transaction()
        if strand == None:
            query = 'START organism = node(%s) MATCH (element:%s)-[r:`PART_OF`]->(organism) RETURN element ORDER BY element.start' % (self.current_organism_ref, type_of_node)
        else:
            query = 'START organism = node(%s) MATCH (element:%s)-[r:`PART_OF`]->(organism) WHERE element.strand="%s" OR element.strand="both" RETURN element ORDER BY element.start' % (self.current_organism_ref, type_of_node, strand)
        print query
        transaction.append(query)
        transaction_out = transaction.execute()
        out = [result.values[0] for result in transaction_out[0]]
        log_message = 'Amount of features to NEXT/OVERLAP %d' % len(out)
        self.logger.info(log_message)
        print log_message
        return out

    def relation_overlap(self, type_of_node='Feature'):
        """
        The function creates relationships 'OVERLAP' in the 'data_base' between nodes with label 'type_of_node'.
        The nodes must have property 'start' as the function compare 2 'start' values
        to make the relationship.
        """
        self.start_time = time()
        if not isinstance(self.current_organism_name, basestring):
            err_message = 'The organism was not chosen. Set current organism before using methods.'
            self.logger.error(err_message)
            raise UserWarning(err_message)

        if len(list(self.data_base.find('Organism', 'name', self.current_organism_name))) == 0:
            err_message = 'There is no node %s in the database or it does not have an "organism" label.' % self.current_organism_name
            self.logger.error(err_message)
            raise ValueError(err_message)

        self.start_time = time()
        features = self._feature_start_ordering(None, type_of_node)
        left_edge = [0]
        for i in xrange(2, len(features)):
            for j in xrange(i-1, max(left_edge)-1, -1):
                if features[j]['end'] >= features[i]['start']:
                    self.data_base.create(rel(features[i], 'OVERLAP', features[j]))
                else:
                    left_edge.append(j)
        log_message = 'Created %d relations for %s in %s' % (len(features), type_of_node, self.current_organism_name)
        self.logger.info(log_message)
        self.stop_time = time()
        print self.stop_time - self.start_time

    def add_label_to_nodes(self, label, label_to_add, search_property=None, search_value=None):
        """
        The function adds new 'label_to_add' to the element with label 'label'.
        """
        if not isinstance(label, basestring) or not isinstance(label_to_add, basestring):
            err_message = 'Arguments must be strings.'
            self.logger.error(err_message)
            raise TypeError(err_message)

        elements = list(self.data_base.find(label, search_property, search_value))
        if len(elements) == 0:
            err_message = 'There are no nodes with %s label, %s property and/or %s value in the database.' % (label, search_property, search_value)
            self.logger.error(err_message)
            raise UserWarning(err_message)

        for element in elements:
            element.add_labels(label_to_add)

    def add_property_to_nodes(self, label, property_to_add, value_to_add, search_property=None, search_value=None):
        """
        Function adds new property 'property_to_add' with value 'value_to_add' to the element with label 'label'.
        """
        if not isinstance(label, basestring) or not isinstance(property_to_add, basestring):
            err_message = 'Arguments "label" and "property_to_add" must be strings.'
            self.logger.error(err_message)
            raise TypeError(err_message)

        elements = list(self.data_base.find(label, search_property, search_value))
        if len(elements) == 0:
            err_message = 'There are no nodes with %s label, %s property and/or %s value in the database.' % (label, search_property, search_value)
            self.logger.error(err_message)
            raise UserWarning(err_message)
        for element in elements:
            properties = element.get_properties()
            properties[property_to_add] = value_to_add
            element.update_properties(properties)

    def relation_similar(self, db, offline = True, similar_quantity = 3, e_value = 0.01):
        """
        The function searches for nodes with labels 'polypeptide', then functions searches among
        incoming relationships for one with 'SIMILAR' label. If it is not found than 'polypeptide'
        'seq' value is BLASTed with blastp function. After that output XML is read.
        For the first 'similar_quantity' homological proteins nodes are made with 'SIMILAR' relationships.
        For each homological new 'polypeptide' several nodes with relation 'LINK' are created with property
        'qualifier' and value 'is_a'. Links labels are 'pdb', 'gb', and/or 'sp'.
        """
        if not isinstance(self.current_organism_name, basestring):
            err_message = 'Set current organism before using methods.'
            self.logger.error(err_message)
            raise UserWarning(err_message)

        if not isinstance(offline, bool):
            err_message = 'Argument must be a bool.'
            self.logger.error(err_message)
            raise TypeError(err_message)

        if not isinstance(similar_quantity, int):
            err_message = 'Argument must be a int.'
            self.logger.error(err_message)
            raise TypeError(err_message)

        if not isinstance(e_value, float):
            err_message = 'Argument must be a float.'
            self.logger.error(err_message)
            raise TypeError(err_message)

        if len(list(self.data_base.find('Organism', 'name', self.current_organism_name))) == 0:
            err_message = 'There is no node %s in the database or it does not have an "Organism" label.' % self.current_organism_name
            self.logger.error(err_message)
            raise ValueError(err_message)

        polypeptides_all = self._search_organism_one_rel('Polypeptide', 0)
        polypeptides_analyzed = self._search_organism_two_rels('Polypeptide', 'Polypeptide', 'SIMILAR', 1)
        polypeptide_counter = 0
        for polypeptide in polypeptides_all:
            if polypeptide in polypeptides_analyzed:
                polypeptide_counter += 1
                print 'Found! ' + str(polypeptide_counter)
                log_message = 'Polypeptide %s is already have similar polypeptides' % polypeptide
                self.logger.info(log_message)
            else:
                protein_txt = open('current_protein.txt', 'w')
                protein_txt.write(str(polypeptide.get_properties()['seq']))
                print ctime()
                protein_txt.close()
                print str(polypeptide.get_properties()['seq'])
                homolog_sequence, gb, pdb, sp, similar_quantity_new = \
                    self._blaster('current_protein.txt', offline, similar_quantity, e_value, db)
                polypeptide_counter += 1
                print 'Blasted:' + str(polypeptide_counter)
                log_message = '%d similar polypeptides are found for polypeptide %s' % (similar_quantity_new, polypeptide)
                print log_message
                self.logger.info(log_message)
                for i in xrange(similar_quantity_new):
                    sim_polypeptide, = self.data_base.create(node())
                    sim_polypeptide.add_labels('Polypeptide')
                    self.data_base.create(rel(sim_polypeptide, ('SIMILAR', {'Date': self.date}), polypeptide))
                    if len(gb[i]) > 0:
                        xref_gb = self._create_xref(gb[i])
                        self.data_base.create(rel(sim_polypeptide, 'REFERENCE_TO', xref_gb))
                        self.data_base.create(rel(xref_gb, 'LINK_TO', self.gb))
                    if len(pdb[i]) > 0:
                        xref_pdb = self._create_xref(pdb[i])
                        self.data_base.create(rel(sim_polypeptide, 'REFERENCE_TO', xref_pdb))
                        self.data_base.create(rel(xref_pdb, 'LINK_TO', self.pdb))
                    if len(sp[i]) > 0:
                        xref_sp = self._create_xref(sp[i])
                        self.data_base.create(rel(sim_polypeptide, 'REFERENCE_TO', xref_sp))
                        self.data_base.create(rel(xref_sp, 'LINK_TO', self.sp))

    def _blaster(self, protein_sequence, offline = True, similar_quantity = 3, e_value = 0.01, db = ''):
        """
        Private method valuecheck suppose to be done by caller method. Never invoke directly!
        """
        if offline == True or len(db) != 0:
            cline = NcbiblastpCommandline(query = protein_sequence,
                                db = db, evalue = e_value,
                                outfmt = 5, out = 'prot_out.xml')
            stdout_str, stderr_str = cline()
            blast_write = open('prot_out.xml')
            blast_read = NCBIXML.read(blast_write)
        else:
            seq_file = open(protein_sequence, 'r')
            protein_sequence_text = seq_file.read()
            seq_file.close()
            #add e_value for online
            blast_write = NCBIWWW.qblast('blastp', 'nr', protein_sequence_text, expect = e_value)
            blast_read = NCBIXML.read(blast_write)
        homolog_seq, gb, pdb, sp = [], [], [], []
        if similar_quantity < 0 or similar_quantity > len(blast_read.alignments):
            similar_quantity = len(blast_read.alignments)
        for i in xrange(similar_quantity):
            gb.append(str(blast_read.alignments[i].title.partition('gb|')[-1].partition('|')[0]))
            pdb.append(str(blast_read.alignments[i].title.partition('pdb|')[-1].partition('|')[0]))
            sp.append(str(blast_read.alignments[i].title.partition('sp|')[-1].partition('|')[0]))
            homolog_seq.append(blast_read.alignments[i].hsps[0].match)
        blast_write.close()
        return homolog_seq, gb, pdb, sp, similar_quantity

    def _current_organism_one_relation_search(self, node_label, relation):
        session = cypher.Session()
        transaction = session.create_transaction()
        query = 'START organism=node(%s) MATCH (organism)-[relation:%s]-(n:%s) RETURN relation' \
                % (self.current_organism_ref, relation, node_label)
        transaction.append(query)
        transaction_out = transaction.execute()
        return [result.values for result in transaction_out[0]]

    def search_one_rel(self, label_out, label_in, relation, return_number = -1):
        """
        Method makes search in the base for pattern (a)-[r]-(b) and returns list with links to nodes and relations [a,b,r].
        Argument 'return_number' filters output list, takes values:
        -1 (return all links), 0 (returns only (a)), 1 (returns only (b)), 2 (returns only (r)).
        """
        if not isinstance(label_in, basestring) or not isinstance(label_out, basestring)\
            or not isinstance(relation, basestring):
            raise TypeError('Arguments  label_out, label_in and relation must be a strings.')

        if not isinstance(return_number, int):
            raise TypeError('Argument must be int')

        if return_number > 2 or return_number < -1:
            raise ValueError('return_number can have values -1, 0, 1 or 2.')

        session = cypher.Session()
        transaction = session.create_transaction()
        if return_number == -1:
            query = 'MATCH (node_out:%s)' \
                    '-[rel:`%s`]-' \
                    ' (node_in:%s) RETURN node_out, node_in, rel' % (label_out, relation, label_in)
            log_message = 'Seqarching query: ' + query
            self.logger.info(log_message)
        if return_number == 0:
            query = 'MATCH (node_out:%s)' \
                    '-[rel:`%s`]-' \
                    ' (node_in:%s) RETURN DISTINCT node_out' % (label_out, relation, label_in)
            log_message = 'Seqarching query: ' + query
            self.logger.info(log_message)
        if return_number == 1:
            query = 'MATCH (node_out:%s)' \
                    '-[rel:`%s`]-' \
                    ' (node_in:%s) RETURN DISTINCT node_in' % (label_out, relation, label_in)
            log_message = 'Seqarching query: ' + query
            self.logger.info(log_message)
        if return_number == 2:
            query = 'MATCH (node_out:%s)' \
                    '-[rel:`%s`]-' \
                    ' (node_in:%s) RETURN DISTINCT rel' % (label_out, relation, label_in)
            log_message = 'Seqarching query: ' + query
            self.logger.info(log_message)
        transaction.append(query)
        transaction_out = transaction.execute()
        if len(transaction_out[0]) == 0:
            log_message = 'Nothing was found on your query'
            self.logger.warning(log_message)
            warnings.warn(log_message)
        if return_number == -1:
            return [result.values for result in transaction_out[0]]
        else:
            return [result.values[0] for result in transaction_out[0]]

    def search_two_rels(self, label_left, label_middle, label_right, relation_lm, relation_mr, return_number = -1):
        """
        Method makes search in the base for pattern (a)-[r1]-(b)-[r2]-(c) and returns list with links to nodes and relations [a,b,c, r1, r2].
        Argument 'return_number' filters output list, takes values:
        -1 (return all links), 0 (returns only (a)), 1 (returns only (b)), 2 (returns only (c)), 3 (returns only (r1)), 4(returns only (r2)).
        """
        if not isinstance(label_left, basestring) or not isinstance(label_middle, basestring)\
            or not isinstance(label_right, basestring) or not isinstance(relation_lm, basestring)\
            or not isinstance(relation_mr, basestring):
            raise TypeError('Arguments  label_left, label_middle, label_right, relation_lm and relation_mr must be a strings.')

        if not isinstance(return_number, int):
            raise TypeError('Argument must be int')

        if return_number > 4 or return_number < -1:
            raise ValueError('return_number can have values -1, 0, 1, 2, 3 or 4.')
        session = cypher.Session()
        transaction = session.create_transaction()
        if return_number == -1:
            query = 'MATCH (label_left:%s)-[relation_lm:`%s`]-(label_middle:%s)-[relation_mr:`%s`]-(label_right:%s) RETURN label_left, label_middle, label_right, relation_lm, relation_mr' \
                    % (label_left, relation_lm, label_middle, relation_mr, label_right)
            log_message = 'Seqarching query: ' + query
            self.logger.info(log_message)
        if return_number == 0:
            query = 'MATCH (label_left:%s)-[relation_lm:`%s`]-(label_middle:%s)-[relation_mr:`%s`]-(label_right:%s) RETURN DISTINCT label_left' \
                    % (label_left, relation_lm, label_middle, relation_mr, label_right)
            log_message = 'Seqarching query: ' + query
            self.logger.info(log_message)
        if return_number == 1:
            query = 'MATCH (label_left:%s)-[relation_lm:`%s`]-(label_middle:%s)-[relation_mr:`%s`]-(label_right:%s) RETURN DISTINCT label_middle' \
                    % (label_left, relation_lm, label_middle, relation_mr, label_right)
            log_message = 'Seqarching query: ' + query
            self.logger.info(log_message)
        if return_number == 2:
            query = 'MATCH (label_left:%s)-[relation_lm:`%s`]-(label_middle:%s)-[relation_mr:`%s`]-(label_right:%s) RETURN DISTINCT label_right' \
                    % (label_left, relation_lm, label_middle, relation_mr, label_right)
            log_message = 'Seqarching query: ' + query
            self.logger.info(log_message)
        if return_number == 3:
            query = 'MATCH (label_left:%s)-[relation_lm:`%s`]-(label_middle:%s)-[relation_mr:`%s`]-(label_right:%s) RETURN DISTINCT relation_lm' \
                    % (label_left, relation_lm, label_middle, relation_mr, label_right)
        if return_number == 4:
            query = 'MATCH (label_left:%s)-[relation_lm:`%s`]-(label_middle:%s)-[relation_mr:`%s`]-(label_right:%s) RETURN DISTINCT relation_mr' \
                    % (label_left, relation_lm, label_middle, relation_mr, label_right)
            log_message = 'Seqarching query: ' + query
            self.logger.info(log_message)
        transaction.append(query)
        transaction_out = transaction.execute()
        if len(transaction_out[0]) == 0:
            log_message = 'Nothing was found on your query: ' + query
            self.logger.warning(log_message)
            warnings.warn(log_message)
        if return_number == -1:
            return [result.values for result in transaction_out[0]]
        else:
            return [result.values[0] for result in transaction_out[0]]