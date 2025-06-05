from __future__ import annotations
import re
from collections import defaultdict
import bisect
import logging
import textwrap
import pathlib

from typing import Literal, Callable, Any, Pattern, TextIO, Type, TypeVar, cast

# Translation table helps to create the reverse complement
_TRANSLATION_TABLE: dict[int, int | None] = str.maketrans('ACTG', 'TGAC')
class Seq(str):
    """Sequence class to hold nucleotide sequence"""
    def reverse_complement(self) -> Seq:
        """Return the reverse complement of the given sequence"""
        return Seq(self[::-1].translate(_TRANSLATION_TABLE))

    def __getitem__(self, *args, **kwargs) -> Seq:
        return Seq(super().__getitem__(*args, **kwargs))

START_CODONS: list[str] = ['TTG', 'CTG', 'ATG']
"""List of start codons"""

STOP_CODONS: list[str] = ['TAA', 'TAG', 'TGA']
"""List of stop codons"""

CODON_TRANSLATION_TABLE: dict[str, str] = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "X",
        "TAG": "X",
        "TGT": "C",
        "TGC": "C",
        "TGA": "X",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
}
"""Table to translate between nucleotide triplets and amino acids."""
for stop_codon in STOP_CODONS:
    CODON_TRANSLATION_TABLE[stop_codon] = '*'

EXTRA_FEATURE_TYPES = {
    'protein_coding_gene': 'gene',
    'ncRNA_gene': 'gene',
}

logger: logging.Logger = logging.getLogger('parseGFF3')

class RevalidationNecessary(Exception):
    """Exception raised to trigger revalidation of the GFF"""
    pass

class Issue:
    """Base type for an issue found within the GFF file. Cannot be used as is, needs to be subclassed."""
    def __init__(self, node) -> None:
        self.node = node

    def fix(self) -> None:
        """Implement to fix the raised issue (if possible)."""
        raise NotImplementedError

    def _debug_log_msg(self, msg) -> None:
        logger.debug(f'{self.node.attributes["ID"]}: {msg}')

    def __message__(self) -> str:
        """Implement to display an informative issue message."""
        raise NotImplementedError

    def __str__(self) -> str:
        return f'Issue with {self.node.type} {self.node.attributes["ID"]} (line {self.node.line_nr}): {self.__message__()}'
    
class FeatureStartLeftOfSequenceRegionStart(Issue):
    """Invoked when a gene or gene feature has a start index smaller than the start of the sequence region it belongs to."""
    def fix(self) -> None:
        self.node.start = self.node.sequence_region.start
        logger.debug('Truncating gene to the start of the sequence region.')
        raise RevalidationNecessary()

    def __message__(self) -> str:
        return f'Feature start ({self.node.start}) is to the left of the start of the containing sequence region ({self.node.sequence_region.name}: {self.node.sequence_region.start}'

class FeatureEndRightOfSequenceRegionEnd(Issue):
    """Invoked when a gene or gene feature has an end index beyond the end of the sequence region it belongs to."""
    def fix(self) -> None:
        self.node.end = self.node.sequence_region.end
        logger.debug('Truncating gene to the end of the sequence region.')
        raise RevalidationNecessary()

    def __message__(self) -> str:
        return f'Feature end ({self.node.end}) is to the right of the end of the containing sequence region ({self.node.sequence_region.name}: {self.node.sequence_region.end}'

class FeatureStartLeftOfParentIssue(Issue):
    """Invoked when a gene or gene feature has a start index smaller than the start of its parent."""
    def __init__(self, node, parent) -> None:
        super().__init__(node)
        self.parent = parent

    def fix(self) -> None:
        self.parent.start = self.node.start
        self._debug_log_msg('Extended start of parent feature.')
        raise RevalidationNecessary()

    def __message__(self) -> str:
        return f'Feature start ({self.node.start}) is to the left of one of its parents ({self.parent.start}).'

class FeatureEndRightOfParentIssue(Issue):
    """Invoked when a gene or gene feature has an end index beyong the end of its parent."""
    def __init__(self, node, parent) -> None:
        super().__init__(node)
        self.parent = parent

    def fix(self) -> None:
        self.parent.end = self.node.end
        self._debug_log_msg('Extended end of parent feature.')
        raise RevalidationNecessary()

    def __message__(self) -> str:
        return f'Feature end ({self.node.end}) is to the right of one of its parents ({self.parent.end}).'

class CDSPhaseIncorrectIssue(Issue):
    """Invoked when the recorded CDS phase differs from the calculated phase."""
    def __init__(self, node, correct_phase) -> None:
        super().__init__(node)
        self.correct_phase = correct_phase

    def fix(self) -> None:
        logger.debug(f'')
        self.node.phase = self.correct_phase
        self._debug_log_msg('Corrected phase of CDS feature.')
        raise RevalidationNecessary()

    def __message__(self) -> str:
        if 'pseudogene' in self.node.parents[0].parents[0].attributes:
            can_be_ignored: str = f' (can be ignored - parent gene is pseudogene)'
        else:
            can_be_ignored = ''
        return f'CDS phase ({self.node.phase}) should be {self.correct_phase}{can_be_ignored}.'

class CodingSequenceStartIsNotStartCodon(Issue):
    """Invoked when a CDS does not start with a start codon.
    
    The fix is to mark the parent gene as a pseudogene."""
    def __message__(self) -> str:
        return f'Coding sequence does not start with a start codon.'

    def fix(self) ->  None:
        self.node.parents[0].mark_as_pseudo()
        self._debug_log_msg('Marked gene as a pseudogene.')
        raise RevalidationNecessary()

class CodingSequenceEndIsNotStopCodon(Issue):
    """Invoked when a CDS does not end with a stop codon.
    
    When fixing, the algorithm tries to find a stop codon within 30nt of the end of the CDS.
    If that cannot be found, the gene will be marked as a pseudogene, unless the CDS contains internal stop codons."""
    def __message__(self) -> str:
        return f'Coding sequence does not end with a stop codon.'

    def fix(self) -> None:
        mRNA = self.node
        contig_seq: Seq = mRNA.sequence_region.sequence

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
    """Invoked when a CDS contains internal stop codons.
    
    When fixing, the algorithm truncates the CDS to the first internal stop codon, unless the CDS length would shrink by more than 10%.
    In that case, the gene is marked as a pseudogene."""
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
    """Invoked when a CDS extends beyond its containing exon.
    
    When fixing, the exon's coordinates are extended to contain the CDS fully."""
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
    """Invoked when two CDSs are overlapping. No automatic fix available."""
    def __init__(self, node, other):
        super().__init__(node)
        self.other = other
    
    def __message__(self):
        return f'CDS overlaps another CDS {self.other.attributes["ID"]}.'    

class Node:
    """Base class of a GFF Node"""
    type: str = "__node__"

    def __init__(
            self,
            line_nr: int,
            sequence_region: SequenceRegion,
            source: str,
            entry_type: str,
            start: int | str,
            end: int | str,
            score: int | str,
            strand: str,
            phase: int | str,
            attributes: str,
            extra_feature_type: str | None = None,
            *args, **kwargs) -> None:
        """Initialise a GFF Node
        
        Arguments:
        line_nr -- the line number the `Node` is defined at in the GFF file
        sequence_region -- the `SequenceRegion` the `Node` belongs to
        source -- the source specifier from the line in the GFF file
        entry_type -- the type of `Node`
        start -- the start coordinate of the `Node`
        end -- the end coordinate of the `Node`
        score -- the score value of the `Node`
        phase -- the phase value of the `Node`
        """
        start: int = int(start)
        end: int = int(end)
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

        self.issues: list[Issue] = []
        self.line_nr = line_nr
        self.sequence_region = sequence_region
        self.source = source
        self.start: int = int(start)
        self.end: int = int(end)
        self.score = score
        self.strand = strand
        self.phase: str | int = phase if phase == '.' else int(phase)
        try:
            self.attributes: dict[str, str] = dict(item.split('=') for item in attributes.split(';'))
        except ValueError:
            raise ValueError(f'Invalid attributes entry on line nr {self.line_nr}.')

        # Allow for GFF-spec-subverting special feature types, such as "protein_coding_gene"
        self.extra_feature_type: str | None = extra_feature_type

        self.children: list[Node] = []
        self.parents: list[Node] = []
        # Assign parents
        if 'Parent' in self.attributes:
            for pid in self.attributes['Parent'].split(','):
                # This assumes that any parents are already read in - should be though because they need to be listed before any children.
                try:
                    self.parents.append(self.sequence_region.node_registry[pid])
                except KeyError:
                    raise ValueError(f'Invalid parent ID {pid} at line nr {self.line_nr}')
            for p in self.parents:
                p.children.append(self)
                if self.strand != p.strand:
                    raise ValueError(f"Strand is different between parent and child at line nr {self.line_nr}.")
        elif not self.toplevel:
            raise ValueError(f'Node has no parent but cannot be top level either on line {self.line_nr}')

        self.sequence_region.add_node(self)

    def add_parent(self, parent: Node):
        if self.strand != parent.strand:
            raise ValueError(f"Strand is different between parent and child.")
        self.parents.append(parent)
        parent.children.append(self)
        self.attributes['Parent'] = ','.join([p.attributes['ID'] for p in self.parents])

    @property
    def sequence(self) -> Seq | None:
        """The nucleotide sequence of the feature (if the sequence region's sequence is available, otherwise `None`)."""
        region_seq: Seq | None = self.sequence_region.sequence
        if region_seq is None:
            return None

        seq: Seq = region_seq[self.start-1:self.end]
        if self.strand == '-':
            seq = seq.reverse_complement()

        return seq
    
    def __len__(self) -> int:
        """Length (in nucleotides) of the feature."""
        return self.end - self.start + 1

    def _validate(self) -> None:
        # Run validation on this node and its children.
        self.issues = []
        self._validate_start_end()
        self.validate()
        for child in self.children:
            child._validate()

    def _validate_start_end(self) -> None:
        # Validate that the start and end of this feature is within the bounds of its parents and the sequence region
        for p in self.parents:
            if self.start < p.start:
                self.issues.append(FeatureStartLeftOfParentIssue(self, p))
            if self.end > p.end:
                self.issues.append(FeatureEndRightOfParentIssue(self, p))
        if self.start < self.sequence_region.start:
            self.issues.append(FeatureStartLeftOfSequenceRegionStart(self))
        if self.end > self.sequence_region.end:
            self.issues.append(FeatureEndRightOfSequenceRegionEnd(self))

    def delete(self) -> None:
        """Delete the feature node and any of its children."""
        for child in list(self.children):
            child.delete()
        for p in self.parents:
            p.children.remove(self)
        del self.sequence_region.node_registry[self.attributes['ID']]

    def __str__(self) -> str:
        attr_str: str = ';'.join([f'{key}={value}' for key, value in self.attributes.items()])
        if self.extra_feature_type is None:
            feature_type = self.type
        else:
            feature_type = self.extra_feature_type
        entries: list[str] = [f'{self.sequence_region.name}\t{self.source}\t{feature_type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{attr_str}']
        for child in self.children:
            entries.append(str(child))
        return '\n'.join(entries)
    
    def __repr__(self) -> str:
        return f'{self.type if self.extra_feature_type is None else self.extra_feature_type} {self.strand}[{self.start}, {self.end}]'

    def apply_recursively(self, func: Callable[[Node], Any]):
        """Apply the given function recursively to the node and its children.
        
        This is done depth-first, meaning any children are iterated over first.

        Arguments
        func -- function to apply (needs to take the `Node` as a single argument, return is ignored)"""
        for child in self.children:
            child.apply_recursively(func)
        return func(self)

    def extend_coordinates(self, new_start=None, new_end=None) -> None:
        """Extend the start and stop coordinates of the feature to the given coordinates.
        
        Crucially, this only allows extending, a feature cannot be shrunk this way (i.e. the new start cannot be to the right of the original).

        Arguments
        new_start -- new start coordinate (default `None`)
        new_end -- new end coordintate (default `None`)
        """
        if new_start is not None:
            self.start = min(self.start, new_start)
        if new_end is not None:
            self.end = max(self.end, new_end)

    def to_dict(self) -> dict[str, Any]:
        """Returns a `dict` containing the feature parameters."""
        return {
            'type': self.type,
            'start': self.start,
            'end': self.end,
            'phase': self.phase,
            'strand': self.strand,
            'score': self.score,
            'attributes': self.attributes,
            'sequence_region': self.sequence_region.name,
        }
    
    def children_of_type(self, node_type: type[NodeType]) -> list[NodeType]:
        children = cast(list[NodeType], self.children)
        return [
            node for node in children
            if node.type == node_type.type
        ]

    def has_child_of_type(self, type: str):
        for child in self.children:
            if child.type == type:
                return True
        return False

NodeType = TypeVar('NodeType', bound=Node)

class GenericNode(Node):
    """A catch-all node type used if `ignore_unknown_feature_types=True`."""
    type: str = '__generic__'
    # The node can be top level.
    toplevel: bool = True
    # This is only used if `ignore_unknown_feature_types` is enabled
    def __init__(self, line_nr, sequence_region, source, entry_type, start, end, score, strand, phase, attributes, *args, **kwargs) -> None:
        self.type = entry_type
        super().__init__(line_nr, sequence_region, source, entry_type, start, end, score, strand, phase, attributes, *args, **kwargs)

    def validate(self) -> None:
        """Stub function, cannot actually validate a generic node since we don't know anything about it."""
        pass

class GeneNode(Node):
    """Node type describing a gene feature."""
    type: str = 'gene'
    toplevel: bool = True
    def validate(self) -> None:
        """Validate the gene node."""
        if self.type != 'gene':
            raise ValueError(f'GeneNode called with wrong type "{self.type}"')
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for genes.')
        if 'Parent' in self.attributes:
            raise ValueError('"Parent" attribute does not make sense for gene.')

    def mark_as_pseudo(self, reason='unknown') -> None:
        """This marks the gene feature as a pseudo gene.
        
        It will delete any CDS, 3'UTR and 5'UTR features found in any mRNA child features.

        Arguments
        reason -- text describing the reason for marking as a pseudogene (default "unknown")
        """
        self.attributes['pseudogene'] = reason
        for child in self.children:
            if child.type == 'mRNA':
                for grandchild in list(child.children):
                    if grandchild.type in ['CDS', 'three_prime_UTR', 'five_prime_UTR']:
                        grandchild.delete()

class MRNANode(Node):
    """Node type describing an mRNA feature."""
    type: str = 'mRNA'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the mRNA node.
        
        This does basic consistence checks on the metadata, checks if the CDS starts with a start codon, ends with a stop codon and has no internal stop codons.
        It also checks if the exons cover the CDSs completely."""
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for mRNA.')

        if len(self.parents) != 1:
            raise ValueError('mRNA needs exactly one parent')

        if self.parents[0].type != 'gene':
            raise ValueError('mRNA parent needs to be gene')

        # Validate CDS child feature phases
        CDSs: list[Node] = self.CDS_children()
        if len(CDSs) > 0:
            phases = self.calculate_CDS_phases()
            for p, CDS in zip(phases, CDSs):
                if p != CDS.phase:
                    self.issues.append(CDSPhaseIncorrectIssue(CDS, p))

        # Validate start / stop codon presence and positon.
        # Only done when the coding sequence is available.
        seq: Seq | None = self.coding_sequence()
        if seq is not None:
            first_codon: Seq = seq[:3]
            last_codon: Seq = seq[3*(len(seq)//3)-3:3*(len(seq)//3)]
            if first_codon not in START_CODONS:
                self.issues.append(CodingSequenceStartIsNotStartCodon(self))
            if last_codon not in STOP_CODONS:
                self.issues.append(CodingSequenceEndIsNotStopCodon(self))

            # Check for internal stop codons
            codons: list[str] = [str(seq[i:i+3]) for i in range(0, len(seq)-3, 3)]
            for stop_codon in STOP_CODONS:
                if stop_codon in codons:
                    self.issues.append(CodingSequenceContainsInternalStopCodons(self))
                    break

        # Checks if the exon child features cover the contained CDSs completely.
        exons: list[Node] = [child for child in self.children if child.type == 'exon']
        for exon in exons:
            for CDS in CDSs:
                if (((CDS.start < exon.start) and (CDS.end > exon.start)) or
                    ((CDS.start < exon.end) and (CDS.end > exon.end))):
                    self.issues.append(ExonDoesNotContainCDSCompletely(exon, CDS))

    def calculate_CDS_phases(self) -> list[int]:
        """Calculates the phase values for all CDS child features."""
        CDSs: list[Node] = self.CDS_children()
        if len(CDSs) == 0:
            return []
        phases: list[int] = [0]
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

    def CDS_children(self) -> list[CDSNode]:
        """Returns a sorted list of CDS child features."""
        CDSs: list[Node] = sorted([child for child in self.children if child.type == 'CDS'], key=lambda c: c.start)
        if self.strand == '-':
            CDSs.reverse()
        return CDSs

    def coding_sequence(self, extend_start_bps=0, extend_end_bps=0) -> Seq | None:
        """Returns the nucleotide coding sequence covered by the contained CDSs."""
        CDSs: list[Node] = self.CDS_children()
        if (len(CDSs) == 0) or (self.sequence_region.sequence is None):
            return None
        return Seq(''.join([str(CDS.sequence) for CDS in CDSs]))
    
    def protein_sequence(self) -> str | None:
        """Returns the amino acid sequence corresponding to the mRNA's coding sequence."""
        coding_sequence: Seq | None = self.coding_sequence()
        if coding_sequence is None:
            return None
        protein_sequence: str = ""
        for codon in textwrap.wrap(coding_sequence, 3):
            try:
                protein_sequence += CODON_TRANSLATION_TABLE[codon]
            except KeyError:
                protein_sequence += 'X'
        return protein_sequence

    def extend_coordinates(self, new_start=None, new_end=None) -> None:
        """Extend the start and stop coordinates of the feature to the given coordinates.
        
        Crucially, this only allows extending, a feature cannot be shrunk this way (i.e. the new start cannot be to the right of the original).
        It also extends the parent's coordinates.

        Arguments
        new_start -- new start coordinate (default `None`)
        new_end -- new end coordintate (default `None`)
        """
        super().extend_coordinates(new_start, new_end)
        self.parents[0].extend_coordinates(new_start, new_end)

class NcRNANode(Node):
    """Node type describing a non-coding RNA feature."""
    type: str = 'ncRNA'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the ncRNA feature."""
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for ncRNA.')

        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('ncRNA parent needs to be gene')

class RRNANode(Node):
    """Node type describing a ribosomal RNA feature."""
    type: str = 'rRNA'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the rRNA feature."""
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for rRNA.')

        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('rRNA parent needs to be gene')

class SnRNANode(Node):
    """Node type describing a small nuclear RNA feature."""
    type: str = 'snRNA'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the snRNA feature."""
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for snRNA.')

        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('snRNA parent needs to be gene')
            
class SnoRNANode(Node):
    """Node type describing a small nucleolar RNA feature."""
    type: str = 'snoRNA'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the snoRNA feature."""
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for snoRNA.')

        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('snoRNA parent needs to be gene')

class ScRNANode(Node):
    """Node type describing a small conditional RNA feature."""
    type: str = 'scRNA'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the scRNA feature."""
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for scRNA.')

        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('scRNA parent needs to be gene')

class TRNANode(Node):
    """Node type describing a t-RNA feature."""
    type: str = 'tRNA'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the t-RNA feature."""
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for tRNA.')

        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('tRNA parent needs to be gene')

class ExonNode(Node):
    """Node type describing an exon."""
    type: str = 'exon'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the exon feature."""
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for exon.')

        for p in self.parents:
            if p.type not in ['mRNA', 'ncRNA', 'rRNA', 'tRNA']:
                raise ValueError('Exon parent needs to be an RNA')

    def extend_coordinates(self, new_start=None, new_end=None) -> None:
        """Extend the start and stop coordinates of the feature to the given coordinates.
        
        Crucially, this only allows extending, a feature cannot be shrunk this way (i.e. the new start cannot be to the right of the original).
        It also extends the parent's coordinates.

        Arguments
        new_start -- new start coordinate (default `None`)
        new_end -- new end coordintate (default `None`)
        """
        super().extend_coordinates(new_start, new_end)
        self.parents[0].extend_coordinates(new_start, new_end)

class ThreePrimeUTRNode(Node):
    """Node type describing a 3' untranslated region (UTR)."""
    type: str = 'three_prime_UTR'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the 3'UTR feature."""
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for three_prime_UTR.')

        for p in self.parents:
            if p.type != 'mRNA':
                raise ValueError('three_prime_UTR parent needs to be mRNA')

class FivePrimeUTRNode(Node):
    """Node type describing a 5' untranslated region (UTR)."""
    type: str = 'five_prime_UTR'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the 5'UTR feature."""
        if self.phase != '.':
            raise ValueError('Phase needs to be "." for five_prime_UTR.')

        for p in self.parents:
            if p.type != 'mRNA':
                raise ValueError('five_prime_UTR parent needs to be mRNA')

class CDSNode(Node):
    """Node type describing a coding sequence feature."""
    type: str = 'CDS'
    toplevel: bool = False
    def validate(self) -> None:
        """Validate the CDS feature."""
        if self.phase not in [0, 1, 2]:
            raise ValueError('Phase needs to be 0, 1 or 2 for CDS.')

        if len(self.parents) != 1:
            raise ValueError('CDS needs exactly one parent')

        if self.parents[0].type != 'mRNA':
            raise ValueError('CDS parent needs to be mRNA')

    def extend_coordinates(self, new_start=None, new_end=None) -> None:
        """Extend the start and stop coordinates of the feature to the given coordinates.
        
        Crucially, this only allows extending, a feature cannot be shrunk this way (i.e. the new start cannot be to the right of the original).
        It also extends the parent's coordinates.

        Arguments
        new_start -- new start coordinate (default `None`)
        new_end -- new end coordintate (default `None`)
        """
        super().extend_coordinates(new_start, new_end)
        self.parents[0].extend_coordinates(new_start, new_end)

class SLASNode(Node):
    """Node type describing a splice leader acceptor site."""
    type: str = 'SLAS'
    toplevel: bool = True
    def validate(self) -> None:
        """Validate a SLAS feature."""
        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('SLAS parent needs to be gene')
        if len(self) != 2:
            raise ValueError('SLAS feature needs to be of length 2')

class PASNode(Node):
    """Node type describing a poly-adenylation site."""
    type: str = 'PAS'
    toplevel: bool = True
    def validate(self) -> None:
        """Validate the PAS feature."""
        for p in self.parents:
            if p.type != 'gene':
                raise ValueError('PAS parent needs to be gene')
        if len(self) != 1:
            raise ValueError('PAS feature needs to be of length 1')

class STARTNode(Node):
    """Node type describing a start codon site."""
    type: str = 'START'
    toplevel: bool = True
    def validate(self) -> None:
        # TODO: Validate that this actually is the site of a start codon if sequence is present.
        pass

class STOPNode(Node):
    """Node type describing a stop codon site."""
    type: str = 'STOP'
    toplevel: bool = True
    def validate(self) -> None:
        # TODO: Validate that this actually is the site of a start codon if sequence is present.
        pass

# Regular expression to check for gaps.
_gap_re: Pattern[str] = re.compile('NNN+')  # Gaps are at least 3 Ns
class SequenceRegion:
    """Describes a sequence region (or contig)."""
    def __init__(self, name: str, start: int, end: int, sequence=None) -> None:
        """Initialize the sequence region.
        
        Arguments
        name -- Name of the contig (no spaces allowed)
        start -- Start coordinate of the contig.
        end -- End coordinate of the contig.
        sequence -- Nucleotide sequence of the contig."""
        if (not issubclass(type(name), str)) or (name == '') or (' ' in name):
            raise ValueError('Name needs to be a valid string without spaces')
        if start > end:
            raise ValueError('Start needs to come before end.')
        if sequence is not None:
            if end - start + 1 != len(sequence):
                #logger.warning(f'Sequence region {name} specification length is different that the sequence length. Fixing sequence end.')
                end = start + len(sequence) - 1

        self.name: str = name
        self.start: int = start
        self.end: int = end
        self.sequence: Seq | None = Seq(sequence) if sequence is not None else None
        self.node_registry: dict[str, Node] = {}

    def validate(self) -> None:
        """Validate the sequence region and all the feature nodes it contains. Checks for CDS feature overlaps."""
        for node in self.node_registry.values():
            if node.type == 'gene':
                node._validate()

        # Check for CDS overlaps
        CDSs: list = sorted([entry for entry in self.node_registry.values() if entry.type == 'CDS'], key=lambda CDS: CDS.start)
        for i, CDS1 in enumerate(CDSs):
            for CDS2 in CDSs[i+1:]:
                if (CDS1.end >= CDS2.start):
                    CDS1.issues.append(CDSOverlap(CDS1, CDS2))
                    CDS2.issues.append(CDSOverlap(CDS2, CDS1))
    
    def add_node(self, node: Node) -> None:
        """Add a given node."""
        if node.attributes['ID'] in self.node_registry:
            raise ValueError(f'Node with ID "{node.attributes["ID"]}" already exists in node registry!')
        if node.sequence_region != self:
            raise ValueError("Node doesn't belong to this sequence region!")
        #if (node.start < self.start) or (node.end > self.end):
        #    raise ValueError('Node is outside of sequence region')

        self.node_registry[node.attributes["ID"]] = node

    def trim_sequence(self, start, end, delete_overlapping_features=True) -> None:
        """Trim the sequence region to the given coordinates.
        
        Arguments
        start -- Trim to this start position
        end -- Trim to this end position
        delete_overlapping_features -- Allow deletion of features that overlap the new sequence region boundaries."""
        logger.debug(f'{self.name}: Trimming sequence region.')
        overlapping_features: list[Node] = [
            node for node in self.node_registry.values()
            if (((node.start <= start) and (node.end >= end)) or
                ((node.start >= start) and (node.end <= end)) or
                ((node.start <= start) and (node.end >= start)) or
                ((node.start <= end) and (node.end >= end)))]
        if (not delete_overlapping_features) and overlapping_features:
            raise NotImplementedError()
        overlapping_genes: list = []
        # BUG!!!!
        # This isn't currently being used to delete any genes, there is a bug!
        def recursively_find_gene(node) -> None:
            if node.type == 'gene':
                overlapping_features.add(node)
            else:
                for p in node.parents[0]:
                    recursively_find_gene(p)
        for gene in overlapping_genes:
            gene.delete()
        
        trim_length: int = end - start + 2
        affected_features: list[Node] = [node for node in self.node_registry.values() if node.start > end]
        for feature in affected_features:
            feature.start -= trim_length
            feature.end -= trim_length
        
        if self.sequence is not None:
            seq: str = str(self.sequence)
            new_seq: str = seq[:start-1] + seq[end+1:]
            self.sequence = Seq(new_seq)

    def closest_node_of_type(self, node: Node, node_types: list[str] | None, direction: str = 'both', strand: Literal['-'] | Literal['+'] | None = None ) -> list[Node]:
        '''Return the closest node of the specified type(s) in the specified direction.
        Please note that the direction is reversed if the given node's strand is '-'.
        '''
        type_filter_func: Callable[[Node], bool] = lambda x: True
        if node_types is not None:
            if not isinstance(node_types, list):
                node_types: list[str] = [node_types]
            if strand is not None:
                type_filter_func = lambda x: x.type in node_types and x.strand == strand
            else:
                type_filter_func = lambda x: x.type in node_types
        elif strand is not None:
            type_filter_func = lambda x: x.strand == strand
        
        # Get a list of all nodes of the given type(s) in the sequence region, sorted by the start cordinate.
        if node.strand == '+':
            nodes: list[Node] = sorted([feature for feature in self.node_registry.values() if type_filter_func(feature) and feature.start >= node.end], key=lambda x: x.start)
        elif node.strand == '-':
            nodes: list[Node] = sorted([feature for feature in self.node_registry.values() if type_filter_func(feature) and feature.end <= node.end], key=lambda x: x.start)

        if node.strand == '-':
            # NOTE: This works well if we don't have overlapping features. If we do, then the behavior is somewhat undefined!
            pass
        found_nodes = []
        if direction == 'forward' or direction == 'both':
            idx_fwd: int = bisect.bisect_right(nodes, node.end, key=lambda x: x.start)
            if idx_fwd < len(nodes):  # No node found otherwise
                found_nodes.append(nodes[idx_fwd])
        if direction == 'backward' or direction == 'both':
            idx_back: int = bisect.bisect_left(nodes[::-1], node.start, key=lambda x: x.end)
            if idx_back < len(nodes):  # No node found otherwise
                #print(len(nodes), idx_back)
                found_nodes.append(nodes[idx_back])
        return sorted(found_nodes, key=lambda x: abs(node.start - x.start))

    def __str__(self) -> str:
        return f'##sequence-region\t{self.name}\t{self.start}\t{self.end}'

    def nodes_of_type(self, node_type: type[NodeType]) -> list[NodeType]:
        return [
            node for node in self.node_registry.values()
            if node.type == node_type.type
        ]


def _read_fasta_file(filename: pathlib.Path | str) -> dict[str, Seq]:
    # Open and read the given FASTA file
    with open(filename, 'r') as f:
        txt: str = f.read()
    return _parse_fasta(txt)

def _parse_fasta(txt: str) -> dict[str, Seq]:
    # Parse a FASTA string into a dict with the keys as the contig names and the values as the sequences.
    contigs: dict[str, Seq] = {}
    fasta_splits: list[str] = re.split('>(.*)\n', txt)[1:]
    for i in range(0, len(fasta_splits), 2):
        name: str = fasta_splits[i].split()[0].strip()
        contigs[name] = Seq(fasta_splits[i+1].replace('\n', ''))
    return contigs

class GffFile:
    """Describes a GFF file."""
    def __init__(self, gff_file=None, fasta_file=None, postpone_validation=True, ignore_unknown_feature_types=False) -> None:
        """Initialize a GFF file object.
        
        If `gff_file` is `None`, creates an empty GFF file, that can optionally be populated with sequences from `fasta_file`.

        Arguments
        gff_file -- File name of the GFF file (default `None`)
        fasta_file -- File name of a FASTA file containing the contig sequences (default `None`)
        postpone_validation --  Whether to postpone running GFF validation, will run validation directly after loading if `False` (default `True`).
        ignore_unknown_feature_types -- Raise an exception if an unknown feature type is encountered if `False` (default `False`).
        """
        if fasta_file is not None:
            self._contigs: dict[str, Seq] = _read_fasta_file(fasta_file)
        else:
            self._contigs = {}

        if gff_file is not None:
            self._read_gff_file(gff_file, ignore_unknown_feature_types)
        else:
            self.sequence_regions: dict[str, SequenceRegion] = self._generate_seqregs()

        if not postpone_validation:
            self.validate()
    
    def _generate_seqregs(self) -> dict[str, SequenceRegion]:
        # Generate empty sequence regions from the given contigs.
        return {name: SequenceRegion(name, 1, len(sequence)+1, sequence) for name, sequence in self._contigs.items()}

    def _read_gff_file(self, filename, ignore_unknown_feature_types) -> None:
        # Parse GFF file
        header_lines: list[str] = []
        GFF_lines: list[str] = []

        # Read GFF file line by line
        with open(filename, 'r') as f:
            f: TextIO
            # Read in header (until we encounter a line not starting with "##")
            for line in f:
                if not line.startswith('##'):
                    break
                header_lines.append(line)

            # Return to the start of the file.
            f.seek(0)
            read_fasta: bool = False
            for line_nr, line in enumerate(f):
                # If we encounter a line starting with "##FASTA", stop reading GFF lines and treat the rest as a FASTA file.
                if line.startswith('##FASTA'):
                    read_fasta = True
                    break
                # Ignore comments
                if line.startswith('#'):
                    continue
                # Record line
                GFF_lines.append((line_nr, line))
            if read_fasta:
                if len(self._contigs) > 0:
                    # If we already read in an external FASTA, use that, but warn the user.
                    logger.warning('External FASTA file provided, but GFF contains FASTA section. Using the external data.')
                else:
                    # If the file has a FASTA portion at the end, use that.
                    self._contigs.update(_parse_fasta(f.read()))

        # New sequence regions
        seqregs: dict[str, SequenceRegion] = self._generate_seqregs()
        # The header contains info on the sequence regions (among other things, but we only use the sequence region info)
        # This is formatted this way: "##sequence-region <name> <start> <stop>"
        for line in header_lines:
            split: list[str] = re.split('\s+', line)
            if split[0] == '##sequence-region':
                name: str
                start: str
                end: str
                name, start, end = split[1:-1]
                if name in seqregs:
                    seqreg: SequenceRegion = seqregs[name]
                    if (int(start) != seqreg.start) or (int(end) != seqreg.end):
                        ValueError(f'FASTA and GFF disagree on start or end of sequence region {name} - start {start} vs {seqreg.start} and end {end} vs {seqreg.end}')
                else:
                    seqregs[split[1]] = SequenceRegion(split[1], int(split[2]), int(split[3]), None)

        # Parse the body of the GFF
        for line_nr, line in GFF_lines:
            # Each line is formatted this way:
            # <sequence ID> <source> <type> <start> <end> <score> <strand> <phase> <attributes>
            # We split the line on whitespace, and make sure that we have nine components
            splits: list[str] = re.split('\s+', line[:-1])
            if len(splits) < 9:
                raise ValueError(f'Malformatted line on line nr {line_nr}.')
            elif len(splits) > 9:
                # We need to treat the last entry in a special way, because the attributes can contain whitespace.
                last_entry: str = ' '.join(splits[8:])
                splits = list(splits[:8]) + [last_entry]
            
            _: str
            seq_id: str
            entry_type: str
            seq_id, _, entry_type, start, end, _, _, _, _ = splits

            try:
                seqreg = seqregs[seq_id]
            except KeyError:
                raise ValueError(f'Unknown sequence region ID on line nr {line_nr}.')
            
            # There are GFF files with non-spec feature types, such as "protein_coding_gene".
            # We treat them as a gene feature, but remember the non-standard type for saving.
            # TODO: We should treat this as a validation / fix issue.
            if entry_type in EXTRA_FEATURE_TYPES:
                extra_feature_type = entry_type
                entry_type = EXTRA_FEATURE_TYPES[entry_type]
            else:
                extra_feature_type = None
            
            try:
                # Generate a new node for the GFF entry.
                # This raises a StopIteration exception if it finds a node type it doesn't know.
                node: Node = next(subclass for subclass in Node.__subclasses__() if subclass.type == entry_type)(line_nr+1, seqreg, splits[1], entry_type, *splits[3:], extra_feature_type=extra_feature_type)
            except StopIteration as e:
                # Produce a GenericNode if instructed to ignore unknown feature types, otherwise raise an Exception.
                if ignore_unknown_feature_types:
                    node: Node = GenericNode(line_nr+1, seqreg, *splits[1:])
                else:
                    raise ValueError(f'Unknown feature type {entry_type} on line nr {line_nr}.')
            except Exception as e:
                raise Exception(f'Exception raised on line nr {line_nr}.')

        self.sequence_regions = seqregs

    def validate(self) -> None:
        """Validate the GFF file, its sequence regions and the contained feature nodes."""
        for seqreg in self.sequence_regions.values():
            seqreg.validate()

    def fix_issues(self) -> None:
        """Recursively fix any issues found during validation that can be fixed automatically.
        
        This needs to be run after `GFFFile.validate`, otherwise no issues have been populated.
        """
        def fix_issues_recursively(node) -> None:
            # Issues are fixed depth-first.
            for child in node.children:
                fix_issues_recursively(child)
            while node.issues:
                issue: Issue = node.issues[0]
                del node.issues[0]
                logger.debug(f'{issue.node.attributes["ID"]}: Fixing issue of type {type(issue).__name__}')
                try:
                    issue.fix()
                except NotImplementedError:
                    logger.debug(f'{issue.node.attributes["ID"]}: Unable to fix.')

        # Go through all sequence regions and run fixes.
        for seqreg in self.sequence_regions.values():
            for gene in [entry for entry in seqreg.node_registry.values() if entry.type == 'gene']:
                while True:
                    try:
                        fix_issues_recursively(gene)
                    except RevalidationNecessary:
                        # On a gene-level, revalidation after a fix might be necessary.
                        logger.debug(f'{gene.attributes["ID"]}: Revalidating after fix')
                        gene._validate()
                        continue
                    break

    def gather_issues(self) -> dict[str, list[Issue]]:
        """Collect issues from all feature nodes. Returns a dict with the issue type as key and a list of occurrences as the value."""
        issues: defaultdict[str, list[Issue]] = defaultdict(list)
        for seqreg in self.sequence_regions.values():
            seqreg.validate()
            for node in seqreg.node_registry.values():
                for issue in node.issues:
                    issues[type(issue).__name__].append(issue)
        return dict(issues)

    def validation_report(self) -> str:
        """Return a report on issues found as a string."""
        issues: dict[str, list[Issue]] = self.gather_issues()
    
        report: str = 'Issue summary:\n'
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
    
    def merge(self, gff: GffFile):
        """Merge another GFF file into this one."""
        srnames_this: set[str] = set(self.sequence_regions)
        srnames_other: set[str] = set(gff.sequence_regions)
        if srnames_this.intersection(srnames_other):
            raise ValueError('Unable to merge because the given GffFile object contains sequence region names that overlap with existing sequence region names.')
        self._contigs.update(gff._contigs)
        self.sequence_regions.update(gff.sequence_regions)

    def __getitem__(self, key: str) -> Node:
        """Return the feature node with the given ID."""
        for seqreg in self.sequence_regions.values():
            try:
                return seqreg.node_registry[key]
            except KeyError:
                pass
        raise KeyError(f'No feature with ID "{key}" found.')

    def _sequence_region_header(self, skip_empty_sequences=True) -> str:
        # Construct the header lines for all sequence regions
        header: str = ''
        for seqreg in self.sequence_regions.values():
            if skip_empty_sequences and (len(seqreg.node_registry) == 0):
                continue
            header += str(seqreg) + '\n'
        return header

    def __str__(self) -> str:
        """Format and return the GFF file as a string."""
        header: str = '''##gff-version    3
##feature-ontology  so.obo
##attribute-ontology    gff3_attributes.obo
'''
        header += self._sequence_region_header()
        body: str = ''
        for _, seqreg in sorted(self.sequence_regions.items(), key=lambda x: x[0]):
            for feature in sorted((entry for entry in seqreg.node_registry.values() if entry.toplevel and not entry.parents), key=lambda x: int(x.start)):
                body += '###\n' + str(feature) + '\n'
        return header + body

    def save(self, filename, include_sequences=False) -> None:
        """Save the GFF file.
        
        Arguments
        filename -- Name of the output file
        include_sequences -- Whether to append the sequences as a FASTA portion (default `False)
        """
        with open(filename, 'w') as f:
            f.write(str(self))
            if include_sequences:
                f.write('##FASTA\n')
                for seqreg in self.sequence_regions.values():
                    if (len(seqreg.node_registry) == 0) or (seqreg.sequence is None):
                        continue
                    seq: str = str(seqreg.sequence)
                    f.write(f'>{seqreg.name}\n')
                    for i in range(0, len(seq)+1, 80):
                        f.write(seq[i:i+80] + '\n')
