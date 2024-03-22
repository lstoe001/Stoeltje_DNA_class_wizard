#!/usr/bin/env python
# coding: utf-8

###HELPFUL VARIABLES


standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


### class seq
class seq:

    def __init__(self, name, organism, sequence, type):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    def info(self):
        print(f"Name: {self.name}")
        print(f"Organism: {self.organism}")
        print(f"Sequence: {self.sequence}")
        print(f"Type: {self.type}")

    def length(self):
        return len(self.sequence)

    def fasta_out(self):
        filename = "{}.fa".format(self.name)
        with open(filename, "w") as f:
            f.write(
                ">"
                + self.name
                + "_"
                + self.organism
                + "_"
                + self.type
                + "\n"
                + self.sequence
            )
            f.close()


### PROTEIN CLASS
class protein(seq):
    aa_mol_weights = {
        "A": 89.09,
        "C": 121.15,
        "D": 133.1,
        "E": 147.13,
        "F": 165.19,
        "G": 75.07,
        "H": 155.16,
        "I": 131.17,
        "K": 146.19,
        "L": 131.17,
        "M": 149.21,
        "N": 132.12,
        "P": 115.13,
        "Q": 146.15,
        "R": 174.2,
        "S": 105.09,
        "T": 119.12,
        "V": 117.15,
        "W": 204.23,
        "X": 0,
        "Y": 181.19,
    }

    def __init__(self, name, organism, sequence, type, size):
        self.size = size
        super().__init__(name, organism, sequence, type)

    def fasta_out(self):
        filename = "{}.fa".format(self.name)
        with open(filename, "w") as f:
            f.write(
                ">"
                + self.name
                + "_"
                + self.organism
                + "_"
                + self.type
                + "_"
                + self.size
                + "\n"
                + self.sequence
            )
            f.close()

    def mol_weight(self):
        total_weight = sum(self.aa_mol_weights.get(aa, 0) for aa in self.sequence)
        print(total_weight)


### CLASS NUCLEOTIDE
class nucleotide(seq):
    def __init__(self, name, organism, sequence, type):

        super().__init__(name, organism, sequence, type)

    def gc_content(self):
        gc_count = self.sequence.count("G") + self.sequence.count("C")
        total_count = len(self.sequence)
        gc_percentage = (gc_count / total_count) * 100
        print("GC content percentage:", gc_percentage)


### class DNA
class DNA(nucleotide):
    def __init__(self, name, sequence, organism, type):
        super().__init__(name, organism, sequence, type)

    def transcribe(self):
        transcribed_sequence = self.sequence.replace("T", "U")
        print(transcribed_sequence)

    def reverse_complement(self):
        complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
        reverse_complement_seq = "".join(
            complement_dict.get(base, base) for base in reversed(self.sequence)
        )
        return reverse_complement_seq

    def six_frames(self):
        frames = []
        # Forward strand frames
        for i in range(3):
            frames.append(self.sequence[i:])
        # Reverse complement frames
        reverse_seq = self.reverse_complement()
        for i in range(3):
            frames.append(reverse_seq[i:])
        print(frames)


### class RNA
class RNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def start(self):
        start_codon = "AUG"
        start_index = self.sequence.find(start_codon)
        if start_index != -1:
            print("Start codon found at index:", start_index)
        else:
            print("Start codon not found.")

    def translate(self):
        start_codon = "AUG"
        start_index = self.sequence.find(start_codon)
        if start_index != -1:
            codons = [
                self.sequence[i : i + 3]
                for i in range(start_index, len(self.sequence), 3)
            ]
            amino_acids = []
            for codon in codons:
                if codon in standard_code:
                    amino_acid = standard_code[codon]
                    if amino_acid == "*":
                        break
                    else:
                        amino_acids.append(amino_acid)
                else:
                    print("Unknown codon:", codon)
            protein_sequence = "".join(amino_acids)
            return protein_sequence
        else:
            print("Start codon not found. Cannot translate.")


### TEST

uidA = DNA(
    name="uidA",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
    organism="Bacteria",
    type="DNA",
)


uidA.fasta_out()


uidA.six_frames()
print(uidA.reverse_complement())


uidA.transcribe()


uidA_RNA = RNA(
    name="uidA_RNA",
    sequence="CGCAUGUUACGUCUGUAGAAACCCCAACCCGUGAAAUCAAAAAA",
    organism="Bacteria",
    type="RNA",
)


uidA_RNA.fasta_out()


uidA_RNA.translate()


uidA_protein = protein(
    name="uidA_protein", sequence="MLRL", organism="Bacteria", type="protein", size="17"
)


uidA_protein.fasta_out()


uidA_protein.mol_weight()
