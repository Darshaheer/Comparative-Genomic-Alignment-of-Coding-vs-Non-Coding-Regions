# Loading Libraries
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
from collections import Counter

# ClustalW Path
Clustal_path = r"Your Path to the File.exe"

#Entrez Email
Entrez.email = "your email"

#Sequences to Download
gene_list = [
    "NM_001126117.2", "NM_001126115.2", "NM_001126118.2",
    "NM_001126112.3", "NM_001126116.2", "NM_001126113.3",
    "NM_001126114.3", "NM_001276698.3", "NM_001276697.3",
    "NM_001276696.3"
]

# ----------------------------------1. Sequence Downloading ------------------------------------------
with Entrez.efetch(db= "nucleotide", id= gene_list, rettype="gb", retmode="text") as file:
    records = list(SeqIO.parse(file, "genbank"))

SeqIO.write(records, "My_sequences.gb", "genbank")
print(f"{len(records)} sequences have been saved to the file")

# ----------------------------------2. Separating CDS & Non-Coding Seq  ------------------------------------------
cds_seq = []
non_cseq = []

for idx, record in enumerate(records):
    full_seq = record.seq

    for feature in record.features:
        if feature.type == "CDS":
            cds = feature.extract(record.seq)
            cds_seq.append(cds)

    non_seq = full_seq
    for feature in record.features:
        if feature.type == "CDS":
            non_seq = non_seq[:feature.location.start] + non_seq[feature.location.end:]
    non_cseq.append(non_seq)


with open("CDS_Sequence.fasta", "w") as rec:
        for i, seq in enumerate(cds_seq):
            rec.write(f">CDS_{i + 1}\n{seq}\n")



with open("Non-Coding_Seq.fasta", "w") as recs:
    for i, seq in enumerate(non_cseq):
        recs.write(f">Non-Coding_{i + 1}\n{seq}\n")

print("CDS and Non-Coding regions have been saved")

# ----------------------------------3. Multiple Sequence Alignment ------------------------------------------
clustalw_cline = ClustalwCommandline(Clustal_path, infile="CDS_Sequence.fasta")
stdout, stderr = clustalw_cline()
print("CDS ALignment Complete")
cds_alignment = AlignIO.read("CDS_Sequence.aln", "clustal")
print("CDS Alignment Length:", cds_alignment.get_alignment_length())

clustalw_cline = ClustalwCommandline(Clustal_path, infile="Non-Coding_Seq.fasta")
stdout, stderr = clustalw_cline()
print("Non Coding Sequences Alignment Complete")
ncs_alignment = AlignIO.read("Non-Coding_Seq.aln", "clustal")
print("CDS Alignment Length:", ncs_alignment.get_alignment_length())

AlignIO.write(cds_alignment, "CDS_Aligned.fasta", "fasta")
AlignIO.write(ncs_alignment, "NonCoding_Aligned.fasta", "fasta")

print("Both CDS and Non-coding alignments saved for further analysis!")

# ----------------------------------4. GC Content Comparison ------------------------------------------
cds_records = list(SeqIO.parse("CDS_Aligned.fasta", "fasta"))
ncs_records = list(SeqIO.parse("NonCoding_Aligned.fasta", "fasta"))
cds_gc = [gc_fraction(rec.seq) * 100 for rec in cds_records]
ncs_gc = [gc_fraction(rec.seq) * 100 for rec in ncs_records]

print("CDS GC%:", cds_gc)
print("Non-coding GC%:", ncs_gc)

# Visualization
plt.figure(figsize=(10,5))
plt.hist(cds_gc, bins=10, alpha=0.6, label="CDS", color="blue", edgecolor="black")
plt.hist(ncs_gc, bins=10, alpha=0.6, label="Non-coding", color="orange", edgecolor="black")
plt.title("GC Content Distribution (CDS vs Non-coding)")
plt.xlabel("GC%")
plt.ylabel("Frequency")
plt.legend()
plt.tight_layout()
plt.show()

# Saving Results
with open("GC_summary.txt", "w") as f:
    f.write("GC Content Comparison\n")
    f.write("CDS Sequences:\n")
    for i, gc in enumerate(cds_gc, 1):
        f.write(f"CDS_{i}: {gc:.2f}%\n")
    f.write("\nNon-coding Sequences:\n")
    for i, gc in enumerate(ncs_gc, 1):
        f.write(f"NonCoding_{i}: {gc:.2f}%\n")

print("GC summary saved to GC_summary.txt")

# ----------------------------------5. Highly Conserved Regions ------------------------------------------
cds_alignment = AlignIO.read("CDS_Aligned.fasta", "fasta")
ncs_alignment = AlignIO.read("NonCoding_Aligned.fasta", "fasta")

def conservation_score(alignment):
    scores = []
    aln_length = alignment.get_alignment_length()
    for i in range(aln_length):  
        column = alignment[:, i]  
        most_common = max(set(column), key=column.count)
        score = column.count(most_common) / len(column)  
        scores.append(score)
    return scores

cds_scores = conservation_score(cds_alignment)
ncs_scores = conservation_score(ncs_alignment)

# Visualization
plt.figure(figsize=(12,6))
plt.plot(cds_scores, label="CDS", color="blue")
plt.plot(ncs_scores, label="Non-coding", color="orange")
plt.title("Conservation Across Sequences (CDS vs Non-coding)")
plt.xlabel("Position in Alignment")
plt.ylabel("Conservation Score (0–1)")
plt.legend()
plt.tight_layout()
plt.show()

# Saving Results
with open("Conservation_summary.txt", "w") as f:
    f.write("Conservation Scores\n\n")
    f.write("CDS:\n")
    f.write(", ".join([f"{s:.2f}" for s in cds_scores]) + "\n\n")
    f.write("Non-coding:\n")
    f.write(", ".join([f"{s:.2f}" for s in ncs_scores]) + "\n")


print("Conservation summary saved to Conservation_summary.txt")
