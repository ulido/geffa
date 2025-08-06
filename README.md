# Geffa
Parse, manipulate and save GFF3 genome annotation files with python code.

## Installation

Installation is quick and straightforward using `pip`:
```bash
python3 -m pip install git+https://github.com/ulido/geffa
```

## Quickstart guide

### Load a GFF file:
```python
from geffa import GffFile

gff = GffFile("annotations.gff")
```
This will load the annotations contained in `annotations.gff` into a [`GffFile`][geffa.geffa.GffFile] object. If it is a combined GFF/FASTA file, it will also load the sequences for all contigs. Alternatively, to load sequences from an external FASTA file `sequences.fasta` you can do:
```python
gff = GffFile("annotations.gff", fasta_file="file.fasta")
```

### Print number of genes per sequence region:
The [`SequenceRegion`][geffa.geffa.SequenceRegion] object can be used to access information and GFF nodes on the sequence region / contig level. All sequence regions of a GFF file are stored within the [GffFile.sequence_regions][geffa.geffa.GffFile] structure:
```python
for sequence_region in gff.sequence_regions.values():
    gene_features = [feature for feature in sequence_region.node_registry.values() if feature.type == 'gene']
    print(sequence_region.name, len(gene_features))
```

### Output nucleic sequences
As an example, here is how to output all UTR sequences in the GFF. This will obviously only work if we have sequence information loaded (either from an included FASTA part or an external FASTA file).
```python
# First write the UTR sequences into a dict
sequences = {}
for sequence_region in gff.sequence_regions.values():
    for UTR in (feature for feature in sequence_region.node_registry.value() if feature.type in ['three_prime_UTR', 'five_prime_UTR']):
        # The returned sequence is already reverse complemented for negative strand features
        sequence = UTR.sequence
        if sequence is not None:
            sequences[UTR.attributes['ID']] = sequence

# Output the sequences into a FASTA file (we could have done this in one step, but it is more instructive like this)
with open('UTR_sequences.fa', 'w') as output:
    for UTR_id, sequence in sequences.items():
        output.write(f">{UTR_id}\n{sequence}\n")
```

### Output amino acid sequences
This snippet will output the amino acid sequences of all genes containing an mRNA feature into a FASTA file. Again, we need to have sequence information loaded for this.
```python
# First write the amino acid sequences into a dict
aa_sequences = {}
for sequence_region in gff.sequence_regions.values():
    for mRNA in (feature for feature in sequence_region.node_registry.value() if feature.type == 'mRNA'):
        aa_sequence = mRNA.protein_sequence()
        if aa_sequence is not None:
            aa_sequences[mRNA.parents[0].attributes['ID']] = aa_sequence

# Output the amino acid sequences into a FASTA file
with open('aa_sequences.fa', 'w') as output:
    for gene_id, aa_sequence in aa_sequences.items():
        output.write(f">{gene_id}\n{aa_sequence}\n")
```