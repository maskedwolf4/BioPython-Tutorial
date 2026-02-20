# Biopython: Complete Tutorial — Beginner to Intermediate

> **Prerequisites:** You know Python basics — variables, loops, functions, file I/O, and ideally a tiny bit of biology.  
> **Goal:** Go from "what is Biopython?" to confidently parsing sequences, querying databases, running alignments, and building real projects.

---

## Table of Contents

1. [Setup & Installation](#1-setup--installation)
2. [Core Concept: The Seq Object](#2-core-concept-the-seq-object)
3. [SeqRecord: Sequences with Metadata](#3-seqrecord-sequences-with-metadata)
4. [SeqIO: Reading & Writing Files](#4-seqio-reading--writing-files)
5. [Working with GenBank Files](#5-working-with-genbank-files)
6. [Querying NCBI with Entrez](#6-querying-ncbi-with-entrez)
7. [BLAST: Finding Similar Sequences](#7-blast-finding-similar-sequences)
8. [Pairwise Sequence Alignment](#8-pairwise-sequence-alignment)
9. [Multiple Sequence Alignment](#9-multiple-sequence-alignment)
10. [Protein Structure with Bio.PDB](#10-protein-structure-with-biopdb)
11. [Phylogenetic Trees with Bio.Phylo](#11-phylogenetic-trees-with-biophylo)
12. [Sequence Motifs](#12-sequence-motifs)
13. [Project Ideas](#13-project-ideas)

---

## 1. Setup & Installation

```bash
# Install Biopython
pip install biopython

# Verify it worked
python -c "import Bio; print(Bio.__version__)"
```

Biopython depends on a few optional tools for certain features:

| Tool | What it's for | Install |
|------|--------------|---------|
| NumPy | Some structure/alignment calculations | `pip install numpy` |
| Matplotlib | Plotting phylogenetic trees | `pip install matplotlib` |
| BLAST+ | Local BLAST searches | Download from NCBI |
| ClustalW / MUSCLE | Multiple alignment | Download separately |

For most beginner/intermediate work, just `pip install biopython` is enough.

---

## 2. Core Concept: The Seq Object

The `Seq` object is the heart of Biopython. Think of it like a Python string with superpowers — it knows it's biological data.

### Creating Sequences

```python
from Bio.Seq import Seq

# DNA sequence
dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

# RNA sequence
rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")

# Protein sequence
protein = Seq("MAIVMGR*KGAR*")

# Seq behaves like a string
print(len(dna))          # 39
print(dna[0:6])          # ATGGCC  (slicing works)
print(dna.count("ATG"))  # count occurrences
print("ATG" in dna)      # True
```

### DNA Operations

```python
from Bio.Seq import Seq

dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

# Complement (A↔T, G↔C)
print(dna.complement())
# TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC

# Reverse complement (5'→3' of opposite strand)
print(dna.reverse_complement())
# CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT

# Transcription: DNA → mRNA (T becomes U)
mrna = dna.transcribe()
print(mrna)
# AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG

# Back-transcription: mRNA → DNA
back = mrna.back_transcribe()
print(back)  # same as original dna
```

### Translation: DNA/RNA → Protein

```python
from Bio.Seq import Seq

dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

# Full translation (includes stop codon as *)
print(dna.translate())
# MAIVMGR*KGAR*

# Stop at first stop codon (most common use)
print(dna.translate(to_stop=True))
# MAIVMGR

# Use a different codon table (e.g., mitochondrial)
print(dna.translate(table="Vertebrate Mitochondrial", to_stop=True))

# Translate RNA directly
rna = dna.transcribe()
print(rna.translate(to_stop=True))
# MAIVMGR  (same result)
```

### GC Content

```python
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

gc = gc_fraction(dna)
print(f"GC content: {gc:.1%}")   # GC content: 53.8%
```

### Joining Sequences

```python
from Bio.Seq import Seq

seq1 = Seq("ATGCGT")
seq2 = Seq("TTTAAA")

combined = seq1 + seq2
print(combined)  # ATGCGTTTTAAA
```

---

## 3. SeqRecord: Sequences with Metadata

In real bioinformatics, a sequence is never just letters — it has an ID, a name, a description, and often annotations (features like genes, promoters, etc.). The `SeqRecord` wraps all of this together.

```python
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Create a SeqRecord manually
record = SeqRecord(
    Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
    id="gene_001",
    name="MyGene",
    description="Hypothetical protein, Homo sapiens",
)

print(record.id)          # gene_001
print(record.name)        # MyGene
print(record.description) # Hypothetical protein, Homo sapiens
print(record.seq)         # ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
print(len(record))        # 39
```

### Annotations

```python
record.annotations["organism"] = "Homo sapiens"
record.annotations["molecule_type"] = "DNA"
record.annotations["date"] = "01-JAN-2024"

print(record.annotations)
# {'organism': 'Homo sapiens', 'molecule_type': 'DNA', 'date': '01-JAN-2024'}
```

### Sequence Features (SeqFeature)

Features describe regions of a sequence — like "this region is a CDS (coding sequence)" or "this is a promoter."

```python
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

record = SeqRecord(Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"), id="gene_001")

# Add a CDS feature from position 0 to 21
cds_feature = SeqFeature(
    FeatureLocation(0, 21),      # start, end (0-based, end exclusive)
    type="CDS",
    qualifiers={"gene": ["myGene"], "product": ["hypothetical protein"]}
)

record.features.append(cds_feature)

# Access features
for feature in record.features:
    print(feature.type, feature.location)
    print(feature.qualifiers)
    
    # Extract the feature sequence
    feature_seq = feature.extract(record.seq)
    print(feature_seq)
```

---

## 4. SeqIO: Reading & Writing Files

`SeqIO` is Biopython's universal file reader/writer. Think of it like `pandas.read_csv()` but for biological sequence formats.

### Supported Formats

| Format | Description | Use |
|--------|-------------|-----|
| `fasta` | Simple sequence + header | Most common |
| `fastq` | Sequences + quality scores | NGS reads |
| `genbank` or `gb` | Rich annotations | NCBI downloads |
| `embl` | European format | EBI databases |
| `clustal` | Alignment format | After alignment |
| `phylip` | Phylogenetics | Tree building |
| `pdb-seqres` | Protein structure sequences | PDB files |

### Reading a FASTA File

```python
from Bio import SeqIO

# SeqIO.parse() returns an iterator — memory-efficient for large files
for record in SeqIO.parse("sequences.fasta", "fasta"):
    print(f"ID: {record.id}")
    print(f"Length: {len(record.seq)}")
    print(f"First 50 bases: {record.seq[:50]}")
    print("---")
```

### Reading a Single Sequence

```python
from Bio import SeqIO

# SeqIO.read() for files with exactly ONE sequence
record = SeqIO.read("single_sequence.fasta", "fasta")
print(record.id)
```

### Load All into a List or Dictionary

```python
from Bio import SeqIO

# Load into a list (for small files)
records = list(SeqIO.parse("sequences.fasta", "fasta"))
print(f"Total sequences: {len(records)}")
print(records[0].id)

# Load into a dictionary (fast lookup by ID)
record_dict = SeqIO.to_dict(SeqIO.parse("sequences.fasta", "fasta"))
my_record = record_dict["gene_001"]  # instant lookup
```

### Writing Sequences

```python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Write a single record
record = SeqRecord(Seq("ATGGCCATTGTAATG"), id="my_gene", description="test")
SeqIO.write(record, "output.fasta", "fasta")

# Write multiple records
records = [
    SeqRecord(Seq("ATGGCC"), id="gene_1", description=""),
    SeqRecord(Seq("TTTAAA"), id="gene_2", description=""),
]
SeqIO.write(records, "output_multi.fasta", "fasta")

print(f"Wrote {len(records)} sequences")
```

### Converting Between Formats

```python
from Bio import SeqIO

# Convert FASTA to GenBank
count = SeqIO.convert("sequences.fasta", "fasta", "sequences.gb", "genbank")
print(f"Converted {count} records")

# Convert GenBank to FASTA
count = SeqIO.convert("sequences.gb", "genbank", "sequences.fasta", "fasta")
```

### Filtering Sequences

```python
from Bio import SeqIO

# Keep only sequences longer than 500 bp
long_sequences = [
    record for record in SeqIO.parse("sequences.fasta", "fasta")
    if len(record.seq) > 500
]

SeqIO.write(long_sequences, "long_sequences.fasta", "fasta")
print(f"Kept {len(long_sequences)} sequences > 500 bp")
```

---

## 5. Working with GenBank Files

GenBank format is the richest format — it stores the sequence plus all known annotations (genes, exons, regulatory regions, taxonomy, references, etc.).

```python
from Bio import SeqIO

record = SeqIO.read("sequence.gb", "genbank")

# Basic info
print(record.id)           # Accession number e.g. NM_007294
print(record.name)         # Locus name
print(record.description)  # Full description

# Organism
print(record.annotations["organism"])
print(record.annotations["taxonomy"])   # Full lineage

# All features
for feature in record.features:
    if feature.type == "CDS":  # coding sequences only
        gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
        product = feature.qualifiers.get("product", ["unknown"])[0]
        location = feature.location
        print(f"Gene: {gene_name} | Product: {product} | Location: {location}")

# Extract a specific feature's sequence
for feature in record.features:
    if feature.type == "CDS":
        cds_seq = feature.extract(record.seq)
        protein = cds_seq.translate(to_stop=True)
        print(f"Protein: {protein[:30]}...")
        break
```

---

## 6. Querying NCBI with Entrez

`Bio.Entrez` lets you search and download data from all NCBI databases (GenBank, PubMed, Taxonomy, Protein, etc.) directly from Python.

### The Rules

> **Always set `Entrez.email`** — NCBI requires it and will block you without it.  
> **Add delays** between requests if downloading many records (use `time.sleep(0.5)`).

### Searching a Database

```python
from Bio import Entrez

Entrez.email = "you@example.com"  # REQUIRED

# Search for BRCA1 gene sequences in humans
handle = Entrez.esearch(
    db="nucleotide",
    term="BRCA1[Gene Name] AND Homo sapiens[Organism] AND mRNA[Filter]",
    retmax=10  # limit results
)
results = Entrez.read(handle)
handle.close()

print(f"Total results: {results['Count']}")
print(f"IDs: {results['IdList']}")
```

### Fetching Records

```python
from Bio import Entrez, SeqIO
import time

Entrez.email = "you@example.com"

# Search
handle = Entrez.esearch(db="nucleotide", term="BRCA1[Gene] AND human[Organism]", retmax=5)
results = Entrez.read(handle)
handle.close()

ids = results["IdList"]

# Fetch each record
for ncbi_id in ids:
    handle = Entrez.efetch(
        db="nucleotide",
        id=ncbi_id,
        rettype="gb",      # GenBank format
        retmode="text"
    )
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    print(f"{record.id}: {record.description[:60]}")
    print(f"  Length: {len(record.seq)} bp")
    print(f"  Organism: {record.annotations.get('organism', 'N/A')}")
    
    time.sleep(0.5)  # be polite to NCBI's servers
```

### Getting Database Info

```python
from Bio import Entrez

Entrez.email = "you@example.com"

# List all available databases
handle = Entrez.einfo()
databases = Entrez.read(handle)
print(databases["DbList"])

# Get details about a specific database
handle = Entrez.einfo(db="nucleotide")
info = Entrez.read(handle)
print(info["DbInfo"]["Description"])
```

### Fetching PubMed Articles

```python
from Bio import Entrez

Entrez.email = "you@example.com"

# Search PubMed for CRISPR papers from 2023
handle = Entrez.esearch(db="pubmed", term="CRISPR Cas9[Title] AND 2023[Year]", retmax=5)
results = Entrez.read(handle)
handle.close()

# Fetch summaries
ids = results["IdList"]
handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
print(handle.read())
handle.close()
```

---

## 7. BLAST: Finding Similar Sequences

BLAST (Basic Local Alignment Search Tool) finds sequences in databases that are similar to your query sequence. You can run it remotely through NCBI or locally if you install BLAST+.

### Remote BLAST (via NCBI)

```python
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

# Load a query sequence
record = SeqIO.read("query.fasta", "fasta")

print("Running BLAST... (this takes 30-60 seconds)")

# Run BLASTn against the nucleotide database
result_handle = NCBIWWW.qblast(
    program="blastn",    # blastn, blastp, blastx, tblastn, tblastx
    database="nt",       # nt, nr, refseq_rna, refseq_protein, swissprot
    sequence=record.seq,
    hitlist_size=10      # return top 10 hits
)

# Save results to file (good practice — avoid re-running)
with open("blast_results.xml", "w") as out:
    out.write(result_handle.read())
result_handle.close()

print("BLAST done! Parsing results...")
```

### Parsing BLAST Results

```python
from Bio.Blast import NCBIXML

# Parse saved results
with open("blast_results.xml") as result_handle:
    blast_records = list(NCBIXML.parse(result_handle))

blast_record = blast_records[0]

print(f"Query: {blast_record.query}")
print(f"Database: {blast_record.database}")
print(f"Number of hits: {len(blast_record.alignments)}\n")

for alignment in blast_record.alignments[:5]:  # top 5
    hsp = alignment.hsps[0]  # best HSP (High-Scoring Pair)
    
    print(f"Hit: {alignment.title[:60]}")
    print(f"  Length: {alignment.length} bp")
    print(f"  Score: {hsp.score}")
    print(f"  E-value: {hsp.expect:.2e}")
    print(f"  Identity: {hsp.identities}/{hsp.align_length} ({hsp.identities/hsp.align_length:.1%})")
    print(f"  Query range: {hsp.query_start}–{hsp.query_end}")
    print(f"  Subject range: {hsp.sbjct_start}–{hsp.sbjct_end}")
    print()
```

### Understanding E-value

| E-value | Interpretation |
|---------|---------------|
| < 1e-100 | Identical sequences |
| 1e-50 to 1e-100 | Near-perfect match |
| 1e-10 to 1e-50 | Strong similarity |
| 1e-3 to 1e-10 | Moderate match |
| > 0.01 | Likely noise — treat with caution |

### BLAST Program Reference

| Program | Query | Database | Use Case |
|---------|-------|----------|----------|
| `blastn` | DNA | DNA | Find similar DNA sequences |
| `blastp` | Protein | Protein | Find similar proteins |
| `blastx` | DNA | Protein | Translate DNA, search proteins |
| `tblastn` | Protein | DNA | Search translated DNA |
| `tblastx` | DNA | DNA | Translated comparison |

---

## 8. Pairwise Sequence Alignment

Alignment finds the best way to line up two sequences to maximize similarity, accounting for insertions and deletions (gaps).

### Global vs Local Alignment

- **Global (Needleman-Wunsch):** Aligns the entire length of both sequences — best when sequences are similar in length and overall structure.
- **Local (Smith-Waterman):** Finds the best-matching region — best when looking for a shared domain within otherwise different sequences.

```python
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

seq1 = "ATCGTTACG"
seq2 = "ATCGTCG"

# Global alignment
print("=== Global Alignment ===")
alignments = pairwise2.align.globalms(
    seq1, seq2,
    2,   # match score
    -1,  # mismatch penalty
    -2,  # gap open penalty
    -0.5 # gap extend penalty
)

for alignment in alignments[:2]:
    print(format_alignment(*alignment))

# Local alignment
print("=== Local Alignment ===")
alignments = pairwise2.align.localms(
    seq1, seq2,
    2, -1, -2, -0.5
)

for alignment in alignments[:2]:
    print(format_alignment(*alignment))
```

### Understanding Alignment Scoring

```python
# pairwise2.align.XXYY where:
# XX = substitution type: 
#   "global" = global alignment
#   "local"  = local alignment
# YY = scoring method:
#   "ms" = match score + mismatch penalty
#   "ds" = use a substitution matrix (BLOSUM, PAM)
#   "dc" = match, mismatch, gap open, gap extend

# For protein alignment, use a substitution matrix
from Bio import pairwise2
from Bio.Align import substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")

protein1 = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPNEECLFLERLEENHEGKQWWNAR"
protein2 = "MAEGEITAFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPNEECLFLERLEENHEGKQWWNAR"

alignments = pairwise2.align.globalds(
    protein1, protein2,
    blosum62,   # substitution matrix
    -10,        # gap open
    -0.5        # gap extend
)

print(pairwise2.format_alignment(*alignments[0]))
```

### Available Substitution Matrices

```python
from Bio.Align import substitution_matrices

# See all available matrices
print(substitution_matrices.load())
# BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90
# PAM30, PAM70, PAM250, BLOSUM45, ...

# BLOSUM62 = standard for protein alignment
# BLOSUM80 = for closely related proteins
# BLOSUM45 = for distantly related proteins
```

---

## 9. Multiple Sequence Alignment

When you have many sequences, you align them all together. Biopython wraps external tools (MUSCLE, ClustalW) and lets you read/write alignment files.

### Reading an Alignment File

```python
from Bio import AlignIO

# Read a multiple sequence alignment in Clustal format
alignment = AlignIO.read("alignment.aln", "clustal")

print(f"Number of sequences: {len(alignment)}")
print(f"Alignment length: {alignment.get_alignment_length()}")
print()

for record in alignment:
    print(f"{record.id:20s} {record.seq[:60]}")
```

### Supported Alignment Formats

```python
# FASTA alignment
alignment = AlignIO.read("alignment.fasta", "fasta")

# Phylip format (used in phylogenetics)
alignment = AlignIO.read("alignment.phy", "phylip")

# Stockholm format (used by Pfam/Rfam)
alignment = AlignIO.read("alignment.sto", "stockholm")
```

### Calculating Alignment Statistics

```python
from Bio import AlignIO
import numpy as np

alignment = AlignIO.read("alignment.aln", "clustal")
aln_len = alignment.get_alignment_length()
n_seqs = len(alignment)

# Conservation at each position
conservation = []
for i in range(aln_len):
    column = alignment[:, i]  # all residues at position i
    most_common = max(set(column), key=column.count)
    freq = column.count(most_common) / len(column)
    conservation.append(freq)

print(f"Average conservation: {sum(conservation)/len(conservation):.2%}")
print(f"Fully conserved positions: {sum(1 for c in conservation if c == 1.0)}")
```

### Running MUSCLE via Python (if installed)

```python
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import subprocess

# Assumes MUSCLE is installed and in PATH
muscle_cline = MuscleCommandline(
    input="sequences.fasta",
    out="aligned.fasta"
)

# Run it
stdout, stderr = muscle_cline()

# Read the result
alignment = AlignIO.read("aligned.fasta", "fasta")
print(f"Aligned {len(alignment)} sequences, length {alignment.get_alignment_length()}")
```

---

## 10. Protein Structure with Bio.PDB

`Bio.PDB` lets you parse PDB files (3D structures of proteins/nucleic acids from the Protein Data Bank).

### The SMCRA Hierarchy

```
Structure
└── Model(s)         (different conformations, usually just Model 0)
    └── Chain(s)     (A, B, C... — different polypeptide chains)
        └── Residue(s)  (amino acids or nucleotides)
            └── Atom(s)  (N, CA, C, O, CB...)
```

### Parsing a PDB File

```python
from Bio.PDB import PDBParser

# Parse a PDB file (QUIET suppresses warnings)
parser = PDBParser(QUIET=True)
structure = parser.get_structure("my_protein", "1abc.pdb")

# Navigate the hierarchy
model = structure[0]           # first (usually only) model
chain_a = model["A"]           # chain A

# Iterate through residues
for residue in chain_a:
    res_name = residue.resname  # e.g. "ALA", "GLY", "TRP"
    res_id = residue.id[1]     # residue number
    print(f"  Residue: {res_name} {res_id}")
```

### Downloading a PDB File Automatically

```python
from Bio.PDB import PDBList

pdbl = PDBList()

# Download 1ABC from RCSB PDB
pdbl.retrieve_pdb_file("1ABC", pdir=".", file_type="pdb")
# Downloads as pdb1abc.ent — rename to 1abc.pdb if needed
```

### Extracting Atoms and Coordinates

```python
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "1abc.pdb")

# Get all alpha-carbon atoms (one per residue)
ca_atoms = []
for atom in structure.get_atoms():
    if atom.name == "CA":
        ca_atoms.append(atom)

print(f"Number of residues (CA atoms): {len(ca_atoms)}")

# Get coordinates of first CA atom
first_ca = ca_atoms[0]
x, y, z = first_ca.get_vector()
print(f"Position of first CA: ({x:.2f}, {y:.2f}, {z:.2f}) Å")
```

### Calculating Distances

```python
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "1abc.pdb")
model = structure[0]

# Get two specific atoms
chain_a = model["A"]
residue_10 = chain_a[10]
residue_20 = chain_a[20]

ca_10 = residue_10["CA"]
ca_20 = residue_20["CA"]

# Distance between alpha carbons
distance = ca_10 - ca_20  # Bio.PDB overloads '-' for distance
print(f"Distance between residue 10 and 20: {distance:.2f} Å")
```

### Getting Sequence from Structure

```python
from Bio.PDB import PDBParser, PPBuilder

parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "1abc.pdb")

# Build polypeptides (continuous segments)
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())
```

---

## 11. Phylogenetic Trees with Bio.Phylo

Phylogenetics studies evolutionary relationships. Biopython can read, write, and visualize trees.

### Reading a Tree

```python
from Bio import Phylo

# Read a Newick format tree
tree = Phylo.read("tree.nwk", "newick")

# Print as ASCII art
Phylo.draw_ascii(tree)
```

### A Newick Tree Looks Like This

```
((Human:0.1, Chimp:0.12):0.3, (Mouse:0.4, Rat:0.38):0.35, Dog:0.6);
```

### Traversing a Tree

```python
from Bio import Phylo

tree = Phylo.read("tree.nwk", "newick")

# All terminals (leaves/taxa)
for clade in tree.get_terminals():
    print(f"Leaf: {clade.name}, Branch length: {clade.branch_length}")

# All internal nodes
for clade in tree.get_nonterminals():
    print(f"Internal node, confidence: {clade.confidence}")

# Find a specific taxon
human = tree.find_any("Human")
print(f"Human branch length: {human.branch_length}")
```

### Tree Statistics

```python
from Bio import Phylo

tree = Phylo.read("tree.nwk", "newick")

# Total branch length
total_length = tree.total_branch_length()
print(f"Total tree length: {total_length:.3f}")

# Distance between two taxa
dist = tree.distance("Human", "Chimp")
print(f"Human-Chimp distance: {dist:.3f}")

# Number of taxa
print(f"Number of taxa: {tree.count_terminals()}")
```

### Plotting with Matplotlib

```python
from Bio import Phylo
import matplotlib.pyplot as plt

tree = Phylo.read("tree.nwk", "newick")

fig, ax = plt.subplots(figsize=(10, 6))
Phylo.draw(tree, axes=ax)
plt.title("Phylogenetic Tree")
plt.tight_layout()
plt.savefig("phylo_tree.png", dpi=150)
plt.show()
```

---

## 12. Sequence Motifs

Motifs are short, recurring sequence patterns — like transcription factor binding sites. `Bio.motifs` handles their creation and searching.

### Creating a Motif from Instances

```python
from Bio import motifs
from Bio.Seq import Seq

# Known binding site instances
instances = [
    Seq("TACAA"),
    Seq("TACGC"),
    Seq("TACAC"),
    Seq("TACCC"),
    Seq("AACCC"),
]

# Create motif
m = motifs.create(instances)

print(m.consensus)     # Most common base at each position
print(m.counts)        # Count matrix
```

### Position-Specific Scoring Matrix (PSSM)

```python
from Bio import motifs
from Bio.Seq import Seq

instances = [Seq("TACAA"), Seq("TACGC"), Seq("TACAC"), Seq("AACCC")]
m = motifs.create(instances)

# Background frequencies (assume equal for simplicity)
background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

# Build PSSM (Position-Specific Scoring Matrix)
pssm = m.counts.normalize(pseudocounts=0.5).log_odds(background)

# Search a sequence for the motif
test_seq = Seq("TACAAGGCTGCATACGCATAACA")
for position, score in pssm.search(test_seq, threshold=3.0):
    print(f"Match at position {position}: score {score:.2f}")
```

---

## 13. Project Ideas

The best way to learn Biopython is to build something real. Here are projects organized by difficulty, with the specific skills each one teaches.

---

### 🟢 Beginner Projects

> Focus: `Seq`, `SeqIO`, basic file parsing, simple scripts

---

#### Project 1: FASTA Statistics Reporter

**What it does:** Takes a FASTA file as input and prints a summary report — number of sequences, min/max/average length, GC content distribution.

**Skills:** `SeqIO.parse()`, `Seq`, `gc_fraction`, file I/O, basic statistics

```python
# Starter skeleton
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def analyze_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    lengths = [len(r.seq) for r in records]
    gc_values = [gc_fraction(r.seq) for r in records]
    
    print(f"Sequences: {len(records)}")
    print(f"Shortest: {min(lengths)} bp")
    print(f"Longest: {max(lengths)} bp")
    print(f"Average GC: {sum(gc_values)/len(gc_values):.1%}")
    # TODO: add histogram of lengths

analyze_fasta("sequences.fasta")
```

**Stretch goals:**
- Output as a CSV report
- Plot a histogram of sequence lengths with matplotlib
- Filter sequences by GC content range and save to a new file

---

#### Project 2: DNA Toolkit

**What it does:** A simple command-line tool that takes a DNA sequence and outputs: complement, reverse complement, mRNA, protein, GC content, and ORF (open reading frame) finder.

**Skills:** `Seq` object, transcription, translation, string operations, argparse

**Stretch goals:**
- Find all 6 reading frames
- Highlight start/stop codons in output
- Accept both raw sequences and FASTA files

---

#### Project 3: Codon Usage Calculator

**What it does:** Reads a set of CDS (coding sequences) in FASTA format and calculates codon usage frequencies for each amino acid.

**Skills:** `SeqIO`, `Seq.translate()`, string parsing, `Counter` from Python stdlib, data output

**Stretch goals:**
- Compare codon usage between two organisms
- Visualize with a bar chart
- Calculate Codon Adaptation Index (CAI) — a more advanced metric

---

#### Project 4: FASTQ Quality Filter

**What it does:** Reads a FASTQ file (NGS reads with quality scores) and filters out reads with average quality below a threshold, saving the passing reads to a new file.

**Skills:** `SeqIO` with fastq format, understanding Phred quality scores, filtering logic

```python
# Key concept: quality scores are stored in record.letter_annotations["phred_quality"]
from Bio import SeqIO

def filter_by_quality(input_file, output_file, min_avg_quality=20):
    passing = []
    for record in SeqIO.parse(input_file, "fastq"):
        quals = record.letter_annotations["phred_quality"]
        avg_qual = sum(quals) / len(quals)
        if avg_qual >= min_avg_quality:
            passing.append(record)
    
    SeqIO.write(passing, output_file, "fastq")
    return len(passing)
```

**Stretch goals:**
- Also filter by minimum length
- Generate a quality report showing the distribution of scores
- Implement adapter trimming

---

#### Project 5: GenBank Feature Extractor

**What it does:** Downloads a GenBank record (by accession number) and extracts all CDS features — gene name, product, location, and sequence — into a CSV file.

**Skills:** `Bio.Entrez`, `SeqIO`, GenBank format, `SeqFeature`, CSV writing

**Stretch goals:**
- Extract features by type (gene, rRNA, tRNA, etc.)
- Translate all CDS and save as a FASTA protein file
- Handle multi-chromosome records

---

### 🟡 Intermediate Projects

> Focus: BLAST, Entrez, alignment, combining tools, real analysis pipelines

---

#### Project 6: Gene Family Finder

**What it does:** Takes a query gene sequence, BLASTs it against NCBI, downloads the top 10 hits, aligns them, and reports percent identity between all pairs.

**Skills:** BLAST, Entrez, `pairwise2`, alignment scoring, combining multiple modules

**Pipeline:**
1. BLAST query → get top hit accession numbers
2. Use Entrez to download full sequences
3. Run pairwise alignments between all pairs
4. Build a similarity matrix and print it

**Stretch goals:**
- Build a phylogenetic tree from the pairwise distances
- Identify conserved domains across the family
- Generate a multiple sequence alignment with MUSCLE

---

#### Project 7: Protein Domain Analyzer

**What it does:** Parses a set of protein structures from PDB, extracts sequences, runs a multiple sequence alignment, and identifies conserved positions.

**Skills:** `Bio.PDB`, `PPBuilder`, `AlignIO`, alignment analysis, `Bio.PDB.DSSP` (secondary structure)

**Stretch goals:**
- Map conserved positions back onto the structure
- Calculate solvent accessibility of conserved residues
- Visualize with PyMOL (external tool integration)

---

#### Project 8: Taxonomic Sequence Downloader

**What it does:** Given a taxonomic group (e.g., "Betacoronavirus"), fetches all available complete genome sequences from NCBI, filters by length and date, and saves them with metadata.

**Skills:** Entrez (esearch + efetch + elink), handling large result sets, rate limiting, metadata parsing

**Pipeline:**
```python
# Search by taxonomy
handle = Entrez.esearch(
    db="nucleotide",
    term="Betacoronavirus[Organism] AND complete genome[Title]",
    retmax=500
)
```

**Stretch goals:**
- Add a progress bar with tqdm
- Store results in a SQLite database
- Build a date-filtered subset for downstream analysis

---

#### Project 9: Promoter Motif Scanner

**What it does:** For a set of co-expressed genes (provided as a list of IDs), fetches their upstream promoter sequences from NCBI, and uses motif discovery to find shared binding sites.

**Skills:** Entrez, `SeqFeature` location math, `Bio.motifs`, PSSM scoring

**Pipeline:**
1. Fetch gene records → find transcription start sites
2. Extract 1000 bp upstream
3. Build a motif from known TFBS
4. Score all promoters with PSSM

**Stretch goals:**
- Integrate with JASPAR database motifs (downloadable as JASPAR format)
- Correct for background nucleotide frequencies
- Output a sequence logo (using logomaker library)

---

#### Project 10: Phylogenetic Tree Builder

**What it does:** Given a gene name, fetches sequences from multiple species, aligns them, infers a distance matrix, and builds + displays a neighbor-joining phylogenetic tree.

**Skills:** Entrez, `AlignIO`, `Bio.Phylo`, `Bio.Phylo.TreeConstruction`, tree visualization

**Pipeline:**
```python
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO, Phylo

# Step 1: fetch sequences from NCBI (multiple species)
# Step 2: run MUSCLE alignment → read with AlignIO
alignment = AlignIO.read("aligned.fasta", "fasta")

# Step 3: calculate distance matrix
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)

# Step 4: build tree
constructor = DistanceTreeConstructor(calculator, "nj")  # neighbor-joining
tree = constructor.build_tree(alignment)

# Step 5: visualize
Phylo.draw_ascii(tree)
```

**Stretch goals:**
- Add bootstrap support values
- Compare NJ tree vs UPGMA tree
- Annotate leaves with organism metadata

---

### 🔴 Advanced Projects

> Focus: Integration, scale, automation, combining Biopython with ML/data science

---

#### Project 11: Genome Annotation Pipeline

**What it does:** Automated pipeline that takes a raw genome FASTA, predicts ORFs, BLASTs predicted proteins, and generates a GFF3 annotation file.

**Skills:** ORF finding, BLAST (local), file format writing, pipeline architecture, subprocess module

**Components:**
- ORF predictor (all 6 reading frames, filter by length)
- Local BLASTp against SwissProt
- GFF3 writer (can use gffutils or write custom)
- Summary statistics report

---

#### Project 12: Structural Similarity Comparator

**What it does:** Downloads a set of PDB structures, calculates pairwise RMSD (structural similarity), and clusters them into structural families.

**Skills:** `Bio.PDB`, `Superimposer`, RMSD calculation, clustering (scipy/sklearn), visualization

```python
from Bio.PDB import Superimposer

# Superimpose structure B onto A
sup = Superimposer()
sup.set_atoms(atoms_fixed, atoms_moving)  # corresponding atoms
sup.apply(structure_b.get_atoms())

print(f"RMSD: {sup.rms:.3f} Å")
```

---

#### Project 13: Variant Effect Predictor

**What it does:** Given a list of SNPs (chromosome position + variant allele), fetches the reference sequence, determines what protein change each SNP causes (synonymous, missense, nonsense), and outputs a structured report.

**Skills:** Entrez, coordinate math, `Seq.translate()`, variant annotation, pandas for output

---

#### Project 14: Metagenomics Taxonomy Classifier

**What it does:** Reads FASTQ reads from a metagenomic sample, BLASTs a subset against NCBI nt, fetches taxonomy info for each hit, and generates a taxonomic composition profile.

**Skills:** BLAST, Entrez taxonomy, result aggregation, visualization (pie charts, Krona-style plots)

---

#### Project 15: Machine Learning Sequence Classifier

**What it does:** Extracts sequence features (k-mer frequencies, GC content, codon usage) from labeled sequences and trains a classifier (random forest or neural network) to predict gene function or organism.

**Skills:** Biopython feature extraction + scikit-learn or PyTorch, model evaluation, cross-validation

```python
# Feature extraction example
from collections import Counter
from Bio.Seq import Seq

def kmer_features(seq, k=3):
    """Extract k-mer frequency vector from a sequence."""
    seq_str = str(seq).upper()
    total = len(seq_str) - k + 1
    counts = Counter(seq_str[i:i+k] for i in range(total))
    return {kmer: count / total for kmer, count in counts.items()}
```

---

## Quick Reference Cheatsheet

### Key Imports

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Entrez, AlignIO, Phylo
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, PDBList, PPBuilder, Superimposer
from Bio.SeqUtils import gc_fraction
from Bio.motifs import create as create_motif
```

### Common Patterns

```python
# Read FASTA
records = list(SeqIO.parse("file.fasta", "fasta"))

# Filter sequences
long = [r for r in records if len(r.seq) > 1000]

# Write sequences
SeqIO.write(long, "long.fasta", "fasta")

# Translate DNA
protein = Seq("ATGGCC...").translate(to_stop=True)

# NCBI search + fetch
Entrez.email = "you@example.com"
handle = Entrez.esearch(db="nucleotide", term="BRCA1[Gene] AND human")
ids = Entrez.read(handle)["IdList"]
record = SeqIO.read(Entrez.efetch(db="nucleotide", id=ids[0], rettype="gb"), "genbank")

# BLAST
result = NCBIWWW.qblast("blastn", "nt", sequence)

# Pairwise alignment
alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -2, -0.5)

# PDB structure
structure = PDBParser(QUIET=True).get_structure("name", "file.pdb")
```

---

## Where to Go Next

| Resource | Link |
|----------|------|
| Official Biopython Tutorial | https://biopython.org/DIST/docs/tutorial/Tutorial.html |
| Biopython Cookbook | https://biopython.org/wiki/Category:Cookbook |
| API Documentation | https://biopython.org/docs/latest/api/ |
| GitHub Examples | https://github.com/biopython/biopython/tree/master/Doc/examples |
| Biostars Community | https://www.biostars.org |
| Rosalind (practice problems) | https://rosalind.info/problems/list-view/ |

> **Tip:** Rosalind.info is the best site for practicing bioinformatics problem-solving step by step — many problems are perfectly matched to Biopython.

---

*Tutorial version: February 2026 | Biopython 1.83+*
