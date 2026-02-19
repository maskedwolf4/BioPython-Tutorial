#Basic
# Your first 5 lines of Biopython
# from Bio.Seq import Seq

# dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
# print(dna.translate()) # → protein sequence - MAIVMGR*KGAR*
# print(dna.complement()) # → complementary strand - TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC
# print(dna.reverse_complement()) # CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT

# Sequence manipulation
# from Bio.Seq import Seq
# from Bio.SeqUtils import gc_fraction

# dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

# # Transcription: DNA → mRNA
# mrna = dna.transcribe() # AUGGCCAUUGUAAUGGGCCGCUGAAAGGGU...

# # Translation: mRNA → protein
# protein = dna.translate() # MAIVMGR*KGAR*
# protein_clean = dna.translate(to_stop=True) # Stops at first stop codon

# # GC Content (0-1 scale)
# gc = gc_fraction(dna) # e.g. 0.5384...
# print(f"Protein: {protein_clean}") # Protein: MAIVMGR
# print(f"GC content: {gc:.1%}") # GC content: 56.4%


# Read a FASTA file

# from Bio import SeqIO

# Parse a FASTA file — returns an iterator of SeqRecord objects
# for record in SeqIO.parse("sequences.fasta", "fasta"):
# print(record.id) # sequence identifier
# print(record.description) # full header
# print(record.seq) # the sequence itself
# print(len(record.seq)) # length

# # Convert FASTA → GenBank format in one line
# SeqIO.convert("sequences.fasta", "fasta", "output.gb", "genbank")

# # Load all into a dict (indexed by ID)
# records = SeqIO.to_dict(SeqIO.parse("sequences.fasta", "fasta"))
# my_gene = records["gene_id_123"]

# BLAST Search
# from Bio.Blast import NCBIWWW, NCBIXML
# from Bio import SeqIO

# Read a sequence from file
# record = SeqIO.read("query.fasta", "fasta")

# # Run BLASTn against NCBI nt database (takes a few seconds)
# result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

# # Parse results
# blast_records = NCBIXML.parse(result_handle)
# for blast_record in blast_records:
# for alignment in blast_record.alignments[:5]: # top 5
# hsp = alignment.hsps[0]
# print(f"Hit: {alignment.title[:50]}")
# print(f" E-value: {hsp.expect}, Score: {hsp.score}")

# Query NCBI
# from Bio.Blast import NCBIWWW, NCBIXML
# from Bio import SeqIO

# # Read a sequence from file
# record = SeqIO.read("query.fasta", "fasta")

# # Run BLASTn against NCBI nt database (takes a few seconds)
# result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

# # Parse results
# blast_records = NCBIXML.parse(result_handle)
# for blast_record in blast_records:
# for alignment in blast_record.alignments[:5]: # top 5
# hsp = alignment.hsps[0]
# print(f"Hit: {alignment.title[:50]}")
# print(f" E-value: {hsp.expect}, Score: {hsp.score}")

#3D Protein
# from Bio.PDB import PDBParser, PPBuilder

# parser = PDBParser(QUIET=True)
# structure = parser.get_structure("protein", "1abc.pdb")

# # Hierarchy: Structure → Model → Chain → Residue → Atom
# model = structure[0] # first model
# chain_a = model["A"] # chain A

# for residue in chain_a.get_residues():
# print(residue.resname, residue.id[1]) # e.g. ALA 23

# # Get all CA (alpha carbon) atom positions
# ca_atoms = [atom for atom in structure.get_atoms() if atom.name == "CA"]
# print(f"{len(ca_atoms)} residues in structure")