import re
from collections import defaultdict
import logging

_TRANSLATION_TABLE = str.maketrans('ACTG', 'TGAC')
class Seq(str):
    def reverse_complement(self):
        return self[::-1].translate(_TRANSLATION_TABLE)

    def __getitem__(self, *args, **kwargs):
        return Seq(super().__getitem__(*args, **kwargs))

START_CODONS = ['TTG', 'CTG', 'ATG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']

logger = logging.getLogger('parseGFF3')

class RevalidationNecessary(Exception):
    pass

class Issue:
    def __init__(self, node):
        self.node = node

    def fix(self):
        raise NotImplementedError

    def _debug_log_msg(self, msg):
        logger.debug(f'{self.node.attributes["ID"]}: {msg}')

    def __message__(self):
        raise NotImplementedError

    def __str__(self):
        return f'Issue with {self.node.type} {self.node.attributes["ID"]} (line {self.node.line_nr}): {self.__message__()}'

class FeatureStartLeftOfSequenceRegionStart(Issue):
    def fix(self):
        self.node.start = self.node.sequence_region.start
        logger.debug('Truncating gene to the start of the sequence region.')
        raise RevalidationNecessary()

    def __message__(self):
        return f'Feature start ({self.node.start}) is to the left of the start of the containing sequence region ({self.node.sequence_region.name}: {self.node.sequence_region.start}'

class FeatureEndRightOfSequenceRegionEnd(Issue):
    def fix(self):
        self.node.end = self.node.sequence_region.end
        logger.debug('Truncating gene to the end of the sequence region.')
        raise RevalidationNecessary()

    def __message__(self):
        return f'Feature end ({self.node.end}) is to the right of the end of the containing sequence region ({self.node.sequence_region.name}: {self.node.sequence_region.end}'

class FeatureStartLeftOfParentIssue(Issue):
    def __init__(self, node, parent):
        super().__init__(node)
        self.parent = parent

    def fix(self):
        self.parent.start = self.node.start
        self._debug_log_msg('Extended start of parent feature.')
        raise RevalidationNecessary()

    def __message__(self):
        return f'Feature start ({self.node.start}) is to the left of one of its parents ({self.parent.start}).'

class FeatureEndRightOfParentIssue(Issue):
    def __init__(self, node, parent):
        super().__init__(node)
        self.parent = parent

    def fix(self):
        self.parent.end = self.node.end
        self._debug_log_msg('Extended end of parent feature.')
        raise RevalidationNecessary()

    def __message__(self):
        return f'Feature end ({self.node.end}) is to the right of one of its parents ({self.parent.end}).'

class CDSPhaseIncorrectIssue(Issue):
    def __init__(self, node, correct_phase):
        super().__init__(node)
        self.correct_phase = correct_phase

    def fix(self):
        logger.debug(f'')
        self.node.phase = self.correct_phase
        self._debug_log_msg('Corrected phase of CDS feature.')
        raise RevalidationNecessary()

    def __message__(self):
        if 'pseudogene' in self.node.parents[0].parents[0].attributes:
            can_be_ignored = f' (can be ignored - parent gene is pseudogene)'
        else:
            can_be_ignored = ''
        return f'CDS phase ({self.node.phase}) should be {self.correct_phase}{can_be_ignored}.'

class CodingSequenceStartIsNotStartCodon(Issue):
    def __message__(self):
        return f'Coding sequence does not start with a start codon.'

    def fix(self):
        self.node.parents[0].mark_as_pseudo()
        self._debug_log_msg('Marked gene as a pseudogene.')
        raise RevalidationNecessary()

class CodingSequenceEndIsNotStopCodon(Issue):
    def __message__(self):
        return f'Coding sequence does not end with a stop codon.'

    def fix(self):
        mRNA = self.node
        contig_seq = mRNA.sequence_region.sequence

        CDSs = mRNA.CDS_children()
        last_CDS = CDSs[-1]
        length = last_CDS.end - last_CDS.start + 1
        codon_aligned_length = 3*(length//3) + last_CDS.phase
        if mRNA.strand == '+':
            end = last_CDS.start + codon_aligned_length - 1
            adjacent_seq = str(contig_seq[end:min(len(contig_seq),end+30*3)])
        else:
            start = last_CDS.end - codon_aligned_length + 1
            adjacent_seq = str(contig_seq[max(0, start-30*3-1):start-1].reverse_complement())
        adjacent_codons = [adjacent_seq[i:i+3] for i in range(0, len(adjacent_seq), 3)]
        for stop_codon in STOP_CODONS:
            try:
                index = adjacent_codons.index(stop_codon)
            except ValueError:
                continue
            # Found a stop codon within 30 codons of the end - adjusting end.
            if mRNA.strand == '+':
                last_CDS.extend_coordinates(new_end=last_CDS.end + 3*index+3)
            else:
                last_CDS.extend_coordinates(new_start=last_CDS.start - 3*index-3)
            self._debug_log_msg('Extended coding sequence to include next stop codon.')
            raise RevalidationNecessary()
        else:
            # Check if we have an issue with internal stop codons first, if not mark as pseudo
            if not any([issubclass(type(issue), CodingSequenceContainsInternalStopCodons) for issue in self.node.issues]):
                self._debug_log_msg('No stop codon found, marking gene as pseudo.')
                mRNA.parents[0].mark_as_pseudo()
            else:
                self._debug_log_msg('No stop codon found, but holding off on marking gene as pseudo because there are internal stop codons present.')

class CodingSequenceContainsInternalStopCodons(Issue):
    def fix(self):
        coding_sequence = self.node.coding_sequence()
        stop_codon_presence = [str(coding_sequence[i:i+3]) in STOP_CODONS for i in range(0, len(coding_sequence)-3, 3)]
        if not any(stop_codon_presence):
            raise RuntimeError('Coding sequence is supposed to have an internal stop, but cannot find it.')
        first_stop_index = stop_codon_presence.index(True)

        old_length = len(coding_sequence)
        new_length = first_stop_index*3+3
        if new_length/old_length < 0.9:
            self._debug_log_msg('Marking as pseudo since length would be truncated more than 10% by first internal stop codon.')
            self.node.parents[0].mark_as_pseudo()
            raise RevalidationNecessary()

        CDSs = self.node.CDS_children()
        if len(CDSs) > 1:
            self._debug_log_msg('Cannot truncate coding sequence with multiple CDSs (needs to be implemented).')
            raise NotImplementedError
        CDS = CDSs[-1]
        self._debug_log_msg('Truncating coding sequence to first stop codon.')
        if self.node.strand == '+':
            CDS.end = CDS.start + 3*first_stop_index + 2
            raise RevalidationNecessary()
        else:
            CDS.start = CDS.end - 3*first_stop_index + 1
            raise RevalidationNecessary()

    def __message__(self):
        return f'Coding sequence contains internal stop codons.'

class ExonDoesNotContainCDSCompletely(Issue):
    def __init__(self, node, CDS):
        super().__init__(node)
        self.CDS = CDS

    def fix(self):
        exon = self.node
        CDS = self.CDS

        new_start = min(exon.start, CDS.start)
        new_end = max(exon.end, CDS.end)
        exon.extend_coordinates(new_start=new_start, new_end=new_end)
        self._debug_log_msg(f'Extended exon coordinates to contain CDS {CDS.attributes["ID"]}.')
        raise RevalidationNecessary()

    def __message__(self):
        return f'CDS extends beyond the range of the exon.'

class CDSOverlap(Issue):
    def __init__(self, node, other):
        super().__init__(node)
        self.other = other
    
    def __message__(self):
        return f'CDS overlaps another CDS {self.other.attributes["ID"]}.'    

class Node:
    def __init__(self, line_nr, sequence_region, source, entry_type, start, end, score, strand, phase, attributes, *args, **kwargs):
        start = int(start)
        end = int(end)
        if start > end:
            raise ValueError(f'Start needs to come before end at line nr {line_nr}.')
        if not issubclass(type(sequence_region), SequenceRegion):
            raise ValueError('sequence_region needs to be of type SequenceRegion.')
        if strand not in ['+', '-', '.']:
            raise ValueError(f'Invalid strand value at line nr {line_nr}')
        if phase not in ['.', '0', '1', '2']:
            raise ValueError(f'Invalid CDS phase value at line nr {line_nr}.')
        if not issubclass(type(attributes), str):
            raise ValueError('attributes needs to be a string.')
        if entry_type != self.type:
            raise ValueError('Called wrong node type!')

        self.issues = []
        self.line_nr = line_nr
        self.sequence_region = sequence_region
        self.source = source
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase if phase == '.' else int(phase)
        try:
            self.attributes = dict(item.split('=') for item in attributes.split(';'))
        except ValueError:
            raise ValueError(f'Invalid attributes entry on line nr {self.line_nr}.')

        self.children = []
        self.parents = []
        if 'Parent' in self.attributes:
            for pid in self.attributes['Parent'].split(','):
                try:
                    self.parents.append(self.sequence_region.node_registry[pid])
                except KeyError:
                    raise ValueError(f'Invalid parent ID at line nr {self.line_nr}')
            for p in self.parents:
                p.children.append(self)
                if self.strand != p.strand:
                    raise ValueError(f"Strand is different between parent and child at line nr {self.line_nr}.")
        self.sequence_region.add_node(self)

    @property
    def sequence(self):
        region_seq = self.sequence_region.sequence
        if region_seq is None:
            return None

        seq = region_seq[self.start-1:self.end]
        if self.strand == '-':
            seq = seq.reverse_complement()

        return seq

    def _validate(self):
        self.issues = []
        self._validate_start_end()
        self.validate()
        for child in self.children:
            child._validate()

    def _validate_start_end(self):
        for p in self.parents:
            if self.start < p.start:
                self.issues.append(FeatureStartLeftOfParentIssue(self, p))
            if self.end > p.end:
                self.issues.append(FeatureEndRightOfParentIssue(self, p))
        if self.start < self.sequence_region.start:
            self.issues.append(FeatureStartLeftOfSequenceRegionStart(self))
        if self.end > self.sequence_region.end:
            self.issues.append(FeatureEndRightOfSequenceRegionEnd(self))

    def delete(self):
        for child in list(self.children):
            child.delete()
        for p in self.parents:
            p.children.remove(self)
        del self.sequence_region.node_registry[self.attributes['ID']]

    def __str__(self):
        attr_str = ';'.join([f'{key}={value}' for key, value in self.attributes.items()])
        entries = [f'{self.sequence_region.name}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{attr_str}']
        for child in self.children:
            entries.append(str(child))
        return '\n'.join(entries)

    def apply_recursively(self, func):
        for child in self.children:
            child.apply_recursively(func)
        return func(self)

    def extend_coordinates(self, new_start=None, new_end=None):
        if new_start is not None:
            self.start = min(self.start, new_start)
        if new_end is not None:
            self.end = max(self.end, new_end)

class GeneNode(Node):
    type = 'gene'
    def validate(self):
        if self.type != 'gene':
            raise ValueError(f'GeneNode called with wrong type "{self.type}"')
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for genes.')
        if 'Parent' in self.attributes:
            raise ValueError('"Parent" attribute does not make sense for gene.')

    def mark_as_pseudo(self, reason='unknown'):
        self.attributes['pseudogene'] = reason
        for child in self.children:
            if child.type == 'mRNA':
                for grandchild in list(child.children):
                    if grandchild.type in ['CDS', 'three_prime_UTR', 'five_prime_UTR']:
                        grandchild.delete()

class MRNANode(Node):
    type = 'mRNA'
    def validate(self):
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for mRNA.')

        if len(self.parents) != 1:
            raise ValueError('mRNA needs exactly one parent')

        if self.parents[0].type != 'gene':
            raise ValueError('mRNA parent needs to be gene')

        # Validate phases
        CDSs = self.CDS_children()
        if len(CDSs) > 0:
            phases = self.calculate_CDS_phases()
            for p, CDS in zip(phases, CDSs):
                if p != CDS.phase:
                    self.issues.append(CDSPhaseIncorrectIssue(CDS, p))

        seq = self.coding_sequence()
        if seq is not None:
            first_codon = seq[:3]
            last_codon = seq[3*(len(seq)//3)-3:3*(len(seq)//3)]
            if first_codon not in START_CODONS:
                self.issues.append(CodingSequenceStartIsNotStartCodon(self))
            if last_codon not in STOP_CODONS:
                self.issues.append(CodingSequenceEndIsNotStopCodon(self))

            codons = [str(seq[i:i+3]) for i in range(0, len(seq)-3, 3)]
            for stop_codon in STOP_CODONS:
                if stop_codon in codons:
                    self.issues.append(CodingSequenceContainsInternalStopCodons(self))
                    break

        exons = [child for child in self.children if child.type == 'exon']
        for exon in exons:
            for CDS in CDSs:
                if (((CDS.start < exon.start) and (CDS.end > exon.start)) or
                    ((CDS.start < exon.end) and (CDS.end > exon.end))):
                    self.issues.append(ExonDoesNotContainCDSCompletely(exon, CDS))

    def calculate_CDS_phases(self):
        CDSs = self.CDS_children()
        if len(CDSs) == 0:
            return []
        phases = [0]
        for i, CDS in enumerate(CDSs[1:]):
            # Calculate phase needed. The phase of the first CDS feature is always p_1=0.
            # (i.e. the first codon of the first CDS always starts at the bp indicated by the start coordinate)
            # The subsequent phases p_i are calculated by taking the length of the preceding CDS L_{i-1}=e_{i-1}-s_{i-1}+1
            # and taking its modulo 3, L_{i-1} % 3, to calculate the number of "overhanging" bps, then
            # subtracting this from 3 to get the number of missing bps having to be filled in at the current CDS,
            # plus the preceding phase p_{i-1} and everything modulo 3 again. In summary:
            #  p_{i} = [ p_{i-1} + 3 - ( e_{i-1} - s_{i-1} + 1 ) % 3 ] % 3
            phases.append((phases[-1] + 3-((CDSs[i].end - CDSs[i].start + 1) % 3)) % 3)
        return phases

    def CDS_children(self):
        CDSs = sorted([child for child in self.children if child.type == 'CDS'], key=lambda c: c.start)
        if self.strand == '-':
            CDSs.reverse()
        return CDSs

    def coding_sequence(self, extend_start_bps=0, extend_end_bps=0):
        CDSs = self.CDS_children()
        if (len(CDSs) == 0) or (self.sequence_region.sequence is None):
            return None
        return Seq(''.join([str(CDS.sequence) for CDS in CDSs]))

    def extend_coordinates(self, new_start=None, new_end=None):
        super().extend_coordinates(new_start, new_end)
        self.parents[0].extend_coordinates(new_start, new_end)

class NcRNANode(Node):
    type = 'ncRNA'
    def validate(self):
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for ncRNA.')

        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('ncRNA parent needs to be gene')

class RRNANode(Node):
    type = 'rRNA'
    def validate(self):
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for rRNA.')

        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('rRNA parent needs to be gene')

class TRNANode(Node):
    type = 'tRNA'
    def validate(self):
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for tRNA.')

        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('tRNA parent needs to be gene')

class ExonNode(Node):
    type = 'exon'
    def validate(self):
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for exon.')

        for p in self.parents:
            if p.type not in ['mRNA', 'ncRNA', 'rRNA', 'tRNA']:
                raise ValueError('Exon parent needs to be an RNA')

    def extend_coordinates(self, new_start=None, new_end=None):
        super().extend_coordinates(new_start, new_end)
        self.parents[0].extend_coordinates(new_start, new_end)

class ThreePrimeUTRNode(Node):
    type = 'three_prime_UTR'
    def validate(self):
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for three_prime_UTR.')

        for p in self.parents:
            if p.type != 'mRNA':
                raise ValueError('three_prime_UTR parent needs to be mRNA')

class FivePrimeUTRNode(Node):
    type = 'five_prime_UTR'
    def validate(self):
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for five_prime_UTR.')

        for p in self.parents:
            if p.type != 'mRNA':
                raise ValueError('five_prime_UTR parent needs to be mRNA')

class CDSNode(Node):
    type = 'CDS'
    def validate(self):
        if self.phase not in [0, 1, 2]:
            raise ValueError('Phase needs to be 0, 1 or 2 for CDS.')

        if len(self.parents) != 1:
            raise ValueError('CDS needs exactly one parent')

        if self.parents[0].type != 'mRNA':
            raise ValueError('CDS parent needs to be mRNA')

    def extend_coordinates(self, new_start=None, new_end=None):
        super().extend_coordinates(new_start, new_end)
        self.parents[0].extend_coordinates(new_start, new_end)

class SLASNode(Node):
    type = 'SLAS'
    def validate(self):
        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('SLAS parent needs to be gene')

class PASNode(Node):
    type = 'PAS'
    def validate(self):
        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('PAS parent needs to be gene')

_gap_re = re.compile('NNN+')  # Gaps are at least 3 Ns
class SequenceRegion:
    def __init__(self, name, start, end, sequence=None):
        if (not issubclass(type(name), str)) or (name == '') or (' ' in name):
            raise ValueError('Name needs to be a valid string without spaces')
        if start > end:
            raise ValueError('Start needs to come before end.')
        if sequence is not None:
            if end - start + 1 != len(sequence):
                #logger.warning(f'Sequence region {name} specification length is different that the sequence length. Fixing sequence end.')
                end = start + len(sequence) - 1

        self.name = name
        self.start = start
        self.end = end
        self.sequence = Seq(sequence) if sequence is not None else None
        self.node_registry = {}

    def validate(self):
        for node in self.node_registry.values():
            if node.type == 'gene':
                node._validate()

        CDSs = sorted([entry for entry in self.node_registry.values() if entry.type == 'CDS'], key=lambda CDS: CDS.start)
        for i, CDS1 in enumerate(CDSs):
            for CDS2 in CDSs[i+1:]:
                if (CDS1.end >= CDS2.start):
                    CDS1.issues.append(CDSOverlap(CDS1, CDS2))
                    CDS2.issues.append(CDSOverlap(CDS2, CDS1))
    
    def add_node(self, node):
        if node.attributes['ID'] in self.node_registry:
            raise ValueError(f'Node with ID "{node.attributes["ID"]}" already exists in node registry!')
        if node.sequence_region != self:
            raise ValueError("Node doesn't belong to this sequence region!")
        #if (node.start < self.start) or (node.end > self.end):
        #    raise ValueError('Node is outside of sequence region')

        self.node_registry[node.attributes["ID"]] = node

    def trim_sequence(self, start, end, delete_overlapping_features=True):
        logger.debug(f'{self.name}: Trimming sequence region.')
        overlapping_features = [
            node for node in self.node_registry.values()
            if (((node.start <= start) and (node.end >= end)) or
                ((node.start >= start) and (node.end <= end)) or
                ((node.start <= start) and (node.end >= start)) or
                ((node.start <= end) and (node.end >= end)))]
        if (not delete_overlapping_features) and overlapping_features:
            raise NotImplementedError()
        overlapping_genes = {}
        def recursively_find_gene(node):
            if node.type == 'gene':
                overlapping_features.add(node)
            else:
                for p in node.parents[0]:
                    recursively_find_gene(p)
        for gene in overlapping_genes:
            gene.delete()
        
        trim_length = end - start + 2
        affected_features = [node for node in self.node_registry.values() if node.start > end]
        for feature in affected_features:
            feature.start -= trim_length
            feature.end -= trim_length
        
        if self.sequence is not None:
            seq = str(self.sequence)
            new_seq = seq[:start-1] + seq[end+1:]
            self.sequence = Seq(new_seq)

    def __str__(self):
        return f'##sequence-region\t{self.name}\t{self.start}\t{self.end}'

def _read_fasta_file(filename):
    with open(filename, 'r') as f:
        txt = f.read()
    return _parse_fasta(txt)

def _parse_fasta(txt):
    contigs = {}
    fasta_splits = re.split('>(.*)\n', txt)[1:]
    for i in range(0, len(fasta_splits), 2):
        name = fasta_splits[i].split()[0].strip()
        contigs[name] = Seq(fasta_splits[i+1].replace('\n', ''))
    return contigs

class GffFile:
    def __init__(self, gff_file, fasta_file=None, postpone_validation=False):
        if fasta_file is not None:
            self._contigs = _read_fasta_file(fasta_file)
        else:
            self._contigs = {}

        self._read_gff_file(gff_file)

        if not postpone_validation:
            self.validate()
    
    def _read_gff_file(self, filename):
        header_lines = []
        GFF_lines = []

        with open(filename, 'r') as f:
            lines = f.readlines()
        
            for line in f:
                if not line.startswith('##'):
                    break
                header_lines.append(line)
            f.seek(0)
            read_fasta = False
            for line_nr, line in enumerate(f):
                if line.startswith('##FASTA'):
                    read_fasta = True
                    break
                if line.startswith('#'):
                    continue
                GFF_lines.append((line_nr, line))
            if read_fasta:
                if len(self._contigs) > 0:
                    logger.warning('External FASTA file provided, but GFF contains FASTA section. Using the external data.')
                else:
                    self._contigs.update(_parse_fasta(f.read()))

        seqregs = {name: SequenceRegion(name, 1, len(sequence)+1, sequence) for name, sequence in self._contigs.items()}
        for line in header_lines:
            split = re.split('\s+', line)
            if split[0] == '##sequence-region':
                name, start, end = split[1:-1]
                if name in seqregs:
                    seqreg = seqregs[name]
                    if (int(start) != seqreg.start) or (int(end) != seqreg.end):
                        raise ValueError('FASTA and GFF disagree on start or end of sequence region')
                else:
                    seqregs[split[1]] = SequenceRegion(split[1], int(split[2]), int(split[3]), None)

        for line_nr, line in GFF_lines:
            splits = re.split('\s+', line[:-1])
            if len(splits) < 9:
                raise ValueError(f'Malformatted line on line nr {line_nr}.')
            elif len(splits) > 9:
                last_entry = ' '.join(splits[8:])
                splits = list(splits[:8]) + [last_entry]
            seq_id, source, entry_type, start, end, score, strand, phase, attributes = splits
            try:
                seqreg = seqregs[seq_id]
            except KeyError:
                raise ValueError(f'Unknown sequence region ID on line nr {line_nr}.')
            try:
                node = next(subclass for subclass in Node.__subclasses__() if subclass.type == entry_type)(line_nr+1, seqreg, *splits[1:])
            except:
                raise Exception(f'Exception raised on line nr {line_nr}.')

        self.sequence_regions = seqregs

    def validate(self):
        for seqreg in self.sequence_regions.values():
            seqreg.validate()

    def fix_issues(self):
        def fix_issues_recursively(node):
            for child in node.children:
                fix_issues_recursively(child)
            while node.issues:
                issue = node.issues[0]
                del node.issues[0]
                logger.debug(f'{issue.node.attributes["ID"]}: Fixing issue of type {type(issue).__name__}')
                try:
                    issue.fix()
                except NotImplementedError:
                    logger.debug(f'{issue.node.attributes["ID"]}: Unable to fix.')

        for seqreg in self.sequence_regions.values():
            for gene in [entry for entry in seqreg.node_registry.values() if entry.type == 'gene']:
                while True:
                    try:
                        fix_issues_recursively(gene)
                    except RevalidationNecessary:
                        logger.debug(f'{gene.attributes["ID"]}: Revalidating after fix')
                        gene._validate()
                        continue
                    break


    def gather_issues(self):
        issues = defaultdict(list)
        for seqreg in self.sequence_regions.values():
            seqreg.validate()
            for node in seqreg.node_registry.values():
                for issue in node.issues:
                    issues[type(issue).__name__].append(issue)
        return issues

    def validation_report(self):
        issues = self.gather_issues()
    
        report = 'Issue summary:\n'
        if len(issues) == 0:
            report += 'No issues found.\n'
        else:
            for issue_type in sorted(issues):
                report += f' * {issue_type}: {len(issues[issue_type])} found\n'

            for issue_type in sorted(issues):
                report += f'\n#### {issue_type} ####\n'
                for issue in issues[issue_type]:
                    report += str(issue) + '\n'

        return report

    def __getitem__(self, key):
        for seqreg in self.sequence_regions.values():
            try:
                return seqreg.node_registry[key]
            except KeyError:
                pass
        raise KeyError(f'No feature with ID "{key}" found.')

    def _sequence_region_header(self, skip_empty_sequences=True):
        header = ''
        for seqreg in self.sequence_regions.values():
            if skip_empty_sequences and (len(seqreg.node_registry) == 0):
                continue
            header += str(seqreg) + '\n'
        return header

    def __str__(self):
        header = '''##gff-version    3
##feature-ontology  so.obo
##attribute-ontology    gff3_attributes.obo
'''
        header += self._sequence_region_header()
        body = ''
        for seqreg in self.sequence_regions.values():
            for gene in (entry for entry in seqreg.node_registry.values() if entry.type == 'gene'):
                body += '###\n' + str(gene) + '\n'
        return header + body

    def save(self, filename, include_sequences=False):
        with open(filename, 'w') as f:
            f.write(str(self))
            if include_sequences:
                f.write('##FASTA\n')
                for seqreg in self.sequence_regions.values():
                    if (len(seqreg.node_registry) == 0) or (seqreg.sequence is None):
                        continue
                    seq = str(seqreg.sequence)
                    f.write(f'>{seqreg.name}\n')
                    for i in range(0, len(seq)+1, 80):
                        f.write(seq[i:i+80] + '\n')
    