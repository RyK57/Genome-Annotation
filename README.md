# Genome-Annotation
Python Streamlit program to analyze and visualize genome sequences 
This code is a tool for annotating genome sequences in FASTA format. It identifies regions with high GC content or repetitive sequences and returns a list of annotations, which are plotted on the sequence using Bokeh library.

The main functions in this code are:

get_protein_type(coding_seq): Determines the protein type (e.g. 'Human' or 'Bacterial') based on the start codon of the coding sequence.
annotate_genome(sequence): Annotates the genome sequence by identifying regions with high GC content or repetitive sequences.
plot_annotations(sequence_len, annotations): Plots the annotations on the sequence.
read_fasta_file(contents): Reads a FASTA file contents and returns a dictionary containing the sequences, coding sequences, non-coding sequences, protein types, and annotations.
The code uses the following Python libraries:

streamlit: A framework for building web applications with Python.
matplotlib: A plotting library for creating static, animated, and interactive visualizations in Python.
Bio.Seq: A module from Biopython that provides classes for representing biological sequences.
Bio.Data: A module from Biopython that provides data dictionaries for various biological data types.
bokeh: A library for creating interactive visualizations in web browsers.
To use this code, you can run it in a Python environment that has these libraries installed. The code can be used to analyze any genome sequence in FASTA format.
