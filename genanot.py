import streamlit as st
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Data import CodonTable
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool

def get_protein_type(coding_seq):
    """
    Determines the protein type (e.g. 'Human' or 'Bacterial') based on the start codon of the coding sequence.
    """
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    start_codons = standard_table.start_codons
    codon = str(coding_seq[:3])
    if codon in start_codons:
        aa = standard_table.forward_table[codon]
        if aa in ["M", "I"]:
            return "Human"
        else:
            return "Bacterial"
    else:
        return "Unknown"

def annotate_genome(sequence):
    """
    Annotates the genome sequence by identifying regions with high GC content or repetitive sequences.
    Returns a list of annotations.
    """
    annotations = []
    gc_threshold = 0.5
    repeat_threshold = 10
    for i in range(0, len(sequence), 100):
        gc_content = sequence[i:i+100].count("G") + sequence[i:i+100].count("C")
        gc_ratio = gc_content / 100
        if gc_ratio > gc_threshold:
            annotations.append({"type": "GC", "start": i, "end": i+100})
        repeat_count = 1
        for j in range(i+1, min(i+100, len(sequence))):
            if sequence[j] == sequence[j-1]:
                repeat_count += 1
            else:
                if repeat_count >= repeat_threshold:
                    annotations.append({"type": "Repeat", "start": j-repeat_count, "end": j})
                repeat_count = 1
    return annotations

def plot_annotations(sequence_len, annotations):
    """
    Plots the annotations on the sequence.
    """
    plot_width = 800
    plot_height = 50
    plot = figure(plot_width=plot_width, plot_height=plot_height, tools="pan,wheel_zoom,box_zoom,reset")
    plot.toolbar.logo = None
    plot.toolbar_location = None
    plot.yaxis.visible = False
    plot.x_range.start = 0
    plot.x_range.end = sequence_len
    source = ColumnDataSource(data=dict(start=[], end=[], type=[]))
    hover = HoverTool(tooltips=[("Type", "@type"), ("Range", "@start - @end")])
    plot.add_tools(hover)
    plot.quad(left="start", right="end", bottom=0, top=1, source=source, color="red")
    for i, annotation in enumerate(annotations):
        source.stream(dict(start=[annotation["start"]], end=[annotation["end"]], type=[annotation["type"]]), rollover=1)
    return plot

def read_fasta_file(contents):
    """
    Reads a FASTA file contents and returns a dictionary containing the sequences.
    """
    sequences = {}
    sequence = ""
    name = ""
    for line in contents.split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if name:
                # annotate the previous sequence and add to dictionary
                coding_seq = Seq(sequence).translate(to_stop=True)
                noncoding_seq = Seq(sequence).translate(to_stop=False)
                annotations = annotate_genome(sequence)
                sequences[name] = {
                    "sequence": sequence,
                    "coding": str(coding_seq),
                    "non-coding": str(noncoding_seq),
                    "coding_type": get_protein_type(coding_seq),
                    "non-coding_type": "N/A",
                    "annotations": annotations
                }
                # plot the annotations
                fig, ax = plt.subplots()
                ax.set_title(f"Annotations for {name}")
                ax.set_xlabel("Sequence Position")
                ax.set_ylim([-1, len(annotations)])
                for i, ann in enumerate(annotations):
                    if ann["type"] == "GC":
                        color = "green"
                    else:
                        color = "red"
                    ax.axvspan(ann["start"], ann["end"], facecolor=color, alpha=0.2)
                    ax.text(ann["start"], i, ann["type"], fontsize=10, ha="left", va="center")
                st.pyplot(fig)
                sequence = ""
            name = line[1:]
        else:
            sequence += line
    # annotate the last sequence and add to dictionary
    coding_seq = Seq(sequence).translate(to_stop=True)
    noncoding_seq = Seq(sequence).translate(to_stop=False)
    annotations = annotate_genome(sequence)
    sequences[name] = {
        "sequence": sequence,
        "coding": str(coding_seq),
        "non-coding": str(noncoding_seq),
        "coding_type": get_protein_type(coding_seq),
        "non-coding_type": "N/A",
        "annotations": annotations
    }
    # plot the annotations
    fig, ax = plt.subplots()
    ax.set_title(f"Annotations for {name}")
    ax.set_xlabel("Sequence Position")
    ax.set_ylim([-1, len(annotations)])
    for i, ann in enumerate(annotations):
        if ann["type"] == "GC":
            color = "green"
        else:
            color = "red"
        ax.axvspan(ann["start"], ann["end"], facecolor=color, alpha=0.2)
        ax.text(ann["start"], i, ann["type"], fontsize=10, ha="left", va="center")
    st.pyplot(fig)
    return sequences


def main():
    st.title("Genome Annotation and Data Visualization")
    st.title("-Rithvik Sabnekar")
    st.write("Upload a text file containing FASTA sequences:")

    file = st.file_uploader("Upload file", type=["txt"])
    st.markdown(f'<a href="https://drive.google.com/file/d/1coCSpNrDI599WICcsCj5cs6C8U9Hb46g/view?usp=sharing" download>Example Input file (Carpodacus Mexicanus)</a>', unsafe_allow_html=True)

    
    if file is not None:
        contents = file.read().decode("utf-8")
        sequences = read_fasta_file(contents)
        st.write("FASTA Sequences:")
        for name, sequence_info in sequences.items():
            st.write("Annotations:")
            st.write(sequence_info["annotations"])

main()
