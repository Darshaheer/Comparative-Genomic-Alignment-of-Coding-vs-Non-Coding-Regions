# Comparative-Genomic-Alignment-of-Coding-vs-Non-Coding-Regions
This project performs comparative genomic analysis of coding and non-coding regions. Sequences are retrieved from NCBI, aligned using ClustalW, and analyzed for GC content distribution and conservation patterns. The pipeline highlights functional constraints on coding regions versus variability in non-coding sequences.

## 🚀 Features  

- Automated sequence retrieval from **NCBI Entrez**  
- Extraction and separation of **CDS vs non-coding regions**  
- **Multiple sequence alignment** using ClustalW2  
- **GC content analysis** and distribution visualization  
- **Conservation scoring** across aligned sequences  
- Results exported to text reports and plots  

---

## 📂 Repository Structure  

├── README.md
├── LICENSE
├── .gitingore
├── Results # Folder with plots (GC, Conservation), Conservation_summary.txt, GC_summary.txt, CDS_Sequence.fasta, Non-Coding_Seq.fasta, CDS_Aligned.fasta, NonCoding_Aligned.fasta, My_sequences.gb
└── Comparative_Analysis.py # Main Python pipeline script

---

## 🛠️ Requirements  

- Python 3.13.5  
- [Biopython](https://biopython.org/)  
- [Matplotlib](https://matplotlib.org/)  
- [ClustalW2](http://www.clustal.org/clustal2/) installed locally  
- Internet connection (for NCBI Entrez sequence download)  

Install dependencies with:  

```bash
pip install biopython matplotlib
📖 Usage
Clone the repository
```
bash
Copy code
git clone https://github.com/your-username/comparative-alignment.git
cd comparative-alignment
Set the ClustalW2 path inside the script (pipeline.py):

python
Copy code
Clustal_path = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
Run the pipeline

bash
Copy code
python pipeline.py

## Outputs generated

Aligned CDS and non-coding sequences (.fasta, .aln)

GC content histograms (.png)

Conservation plots (.png)

Summary reports (GC_summary.txt, Conservation_summary.txt)

## Example Results
GC Content Distribution
Compares base composition between CDS (blue) and non-coding regions (orange).

Conservation Plot
Shows evolutionary conservation along sequence alignments, highlighting stronger conservation in CDS regions compared to non-coding.

## Biological Insight
Coding regions (CDS): Show higher conservation and distinct GC content due to selective pressure to maintain protein structure and function.

Non-coding regions: Display greater variability, reflecting less evolutionary constraint.

This analysis can help in identifying functional regions, studying evolutionary divergence, and distinguishing coding vs regulatory roles.

## Notes
Replace the Gene Id's with your desired gene Id's.

## LICENSE
MIT License This project is open-source. You are free to use and modify it for research or educational purposes.

