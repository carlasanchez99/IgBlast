import re
import collections
import csv
import pandas as pd  # Import pandas for data manipulation
from Bio.Seq import Seq
import os
import subprocess
import xlsxwriter
from Bio.Align import PairwiseAligner
import numpy as np



class QueryParser:
    """Parses out the basic query information."""

    def __init__(self):
        self.regex = re.compile(r'^Query= ((\S| )*)')

    def parse(self, line, out_d):
        matches = re.match(self.regex, line)
        if matches:
            # Store only the sequence ID in the output dictionary
            out_d['mAb name'] = matches.group(1)  # Store the query
            return out_d
        else:
            return False
        
class AlignmentLine:
    """Helper class for parsing each line of the alignments result."""
    
    def __init__(self, line):
        self.id = None
        self.unique_key = None
        self.start = None
        self.end = None
        self.is_query = False
        self.left = ''
        self.type = ''
        self.gene_type = ''
        self.percent_identity = ''
        self.fraction = ''
        self.width = None
        self.span = None
        self.is_header = False
        self.is_translation = False
        self.appended_count = 0
        
        self.hit_regex = re.compile(r'([VDJ])\s+(\S*)\s+(\S*)\s+(\S*)\s+([0-9]+)\s+(\S*)\s+([0-9]+)')
        self.query_regex = re.compile(r'(\S*Query\S*)\s+([0-9]+)\s+(\S*)\s+([0-9]+)')
        self.header_regex = re.compile(r'[<\->]')
        
        # Parse the line using read_line
        self.read_line(line)

    def read_line(self, line):
        # Query Line
        matches = re.search(self.query_regex, line)
        if matches:
            self.id = matches.group(1)
            self.start = int(matches.group(2))
            self.left = self.id + '  ' + str(self.start)
            self.al_string = matches.group(3)
            self.end = int(matches.group(4))
            self.is_query = True
            self.span = matches.span(3)
            self.width = self.span[1] - self.span[0]
            self.unique_key = self.id + '-' + str(self.end)
            return

        # Header Line
        matches = re.search(self.header_regex, line)
        if matches:
            self.id = 'header'
            self.unique_key = self.id
            self.is_header = True
            self.al_string = line
            self.start = ''
            self.end = ''
            return

        # If we made it here we must be a translation
        self.id = 'translation'
        self.unique_key = self.id
        self.is_translation = True
        self.al_string = line
        self.start = ''
        self.end = ''

class SignificantAlignmentParser():
    """Parses out information about each hit and alignment for significant matches."""

    def __init__(self):
        self.required = True
        self.failure_regex = re.compile(r'^.+No hits found.+$')
        self.trigger_regex = re.compile(r'^Sequences producing significant alignments')
        self.hit_regex = re.compile(r'(.*?)[ ]+([0-9.\-e]+)[ ]+([0-9.\-e]+)')  # Captures gene, bit_score, e_value
        self.halt_regex = re.compile(r'^Domain classification requested')
        self.triggered = False
    def parse(self, line, out_d):
        # Check for the halt condition
        if re.match(self.halt_regex, line):
            self.triggered = False
            return out_d
        # When triggered, capture the first hit's e-value
        if self.triggered:
            matches = re.match(self.hit_regex, line)
            if matches:
                # Only set e-value if it hasn't been set already
                if out_d.get('E-value') is None:  # Check if 'E-value' is None
                    out_d['E-value'] = float(matches.group(3))  # Store e_value directly
                if out_d.get('RefGene') is None:  # Check if 'E-value' is None
                    out_d['RefGene'] = (matches.group(1))  # Store e_value directly
                return out_d  # Return immediately after capturing the first e-value

        # Start parsing significant alignments
        if re.match(self.trigger_regex, line):
            self.triggered = True  # Set triggered state to True for subsequent lines
            return out_d
        return False
    
dictionary_aa = {
        'A': 1.8,   # Alanine
        'R': -4.5,  # Arginine
        'N': -3.5,  # Asparagine
        'D': -3.5,  # Aspartic acid
        'C': 2.5,   # Cysteine
        'E': -3.5,  # Glutamic acid
        'Q': -3.5,  # Glutamine
        'G': -0.4,  # Glycine
        'H': -3.2,  # Histidine
        'I': 4.5,   # Isoleucine
        'L': 3.8,   # Leucine
        'K': -3.9,  # Lysine
        'M': 1.9,   # Methionine
        'F': 2.8,   # Phenylalanine
        'P': -1.6,  # Proline
        'S': -0.8,  # Serine
        'T': -0.7,  # Threonine
        'W': -0.9,  # Tryptophan
        'Y': -1.3,  # Tyrosine
        'V': 4.2    # Valine
    }



# VDJSummaryParser class to parse the VDJ summary details
class VDJSummaryParser:
    """Parses out VDJ summary information."""

    def __init__(self):
        self.regex = re.compile(r'^V-\(D\)-J rearrangement summary for query sequence \((.*)\)\.')
        self.triggered = False
        self.fields = None
    def parse(self, line, out_d):
        if self.triggered:
            item_d = dict(zip(self.fields, line.strip().split('\t')))
            for key, val in item_d.items():
                stripped_key = key.strip()
                # Set families and e-values based on gene matches
                if stripped_key == 'Top V gene match':
                    v_gene_matches = re.findall(r'V(\d+D?-?\d+)', val)
                    out_d['V+allele'] = val
                    if v_gene_matches:
                        # Eliminar duplicados con set y ordenar para mantener el orden original
                        unique_v_genes = sorted(set(v_gene_matches), key=v_gene_matches.index)
                        new_vkl_value = '/'.join(unique_v_genes)  # Genera la nueva coincidencia de Vkl
                        # Si 'Vkl' no está en out_d o su valor es diferente al nuevo, entonces lo agregamos
                        if 'V' not in out_d or out_d['V'] != new_vkl_value:
                            out_d['V'] = new_vkl_value
                if stripped_key == 'Top D gene match':
                    out_d['D+allele'] = val

                    if '-' in val and '*' in val and '/' not in val:
                        d_gene_matches = re.findall(r'D(\d+D?-?\d+)', val)
                    elif '/' in val and '*' in val:
                        val = re.sub(r',\s+', ',', val)
                        # Expresión regular ajustada para extraer las partes relevantes
                        matches = re.findall(r'I[GDJH][A-Za-z]*(\d+)/([A-Za-z0-9\-]+)\*\d+', val)

                        # Verificar y formatear las coincidencias
                        if matches:
                            d_gene_matches = [f"{match[0]}{match[1]}" for match in matches]
                    else:
                        d_gene_matches="NA"
                    if d_gene_matches != "NA" :
                        # Eliminar duplicados con set y ordenar para mantener el orden original
                        unique_d_genes = sorted(set(d_gene_matches), key=d_gene_matches.index)
                        new_dh_value = '/'.join(unique_d_genes)  # Genera la nueva coincidencia de Vkl
                        # Si 'Vkl' no está en out_d o su valor es diferente al nuevo, entonces lo agregamos
                        if 'DH' not in out_d or out_d['DH'] != new_dh_value:
                            out_d['DH'] = new_dh_value
            # Manejo de los diferentes valores de Top J gene match
                elif stripped_key == 'Top J gene match':
                # Buscar todas las coincidencias que sigan el patrón J con números
                    j_gene_matches = re.findall(r'J(\d+)', val)
                    out_d['J+allele'] = val
                    if j_gene_matches:
                        # Eliminar duplicados con set y ordenar para mantener el orden original
                        unique_j_genes = sorted(set(j_gene_matches), key=j_gene_matches.index)
                        new_jkl_value = '/'.join(unique_j_genes)  # Genera la nueva coincidencia de Vkl
                        # Si 'Vkl' no está en out_d o su valor es diferente al nuevo, entonces lo agregamos
                        if 'J' not in out_d or out_d['J'] != new_jkl_value:
                            out_d['J'] = new_jkl_value
                   
            # Set default values if matches are not found
                elif stripped_key == 'Top C gene match':
        # Tomar solo el primer valor antes de la coma
                    first_c_gene = val.split(',')[0].strip()  # Solo el primer valor, quitando espacios extra
                    out_d['Subtype'] = first_c_gene
                elif stripped_key == 'Chain type':
        # Tomar solo el primer valor antes de la coma
                    out_d['VK/VL/VH'] = val
                elif stripped_key == 'V-J frame':
                    # Verificar si está en frame o out of frame
                    if "In-frame" in val:
                        out_d['In frame'] = 'YES'
                    elif "Out-of-frame" in val:
                        out_d['In frame'] = 'NO'
                    else:
                        out_d['In frame'] = 'NO'
                else:
                    value = val.split(',')[0]
                    out_d[stripped_key] = value
            self.triggered = False
            return out_d

        matches = re.match(self.regex, line)
        if matches:
            self.fields = matches.group(1).split(',')
            self.triggered = True
            return out_d
        return False

    @staticmethod
    def set_family(prop, value, output):
        try:
            star_index = value.index('*')
            output[prop] = value[0:star_index]
        except ValueError:
            output[prop] = value

class Cregion_ext5:
    """Parses out VDJ summary information, including nucleotide sequences."""
    def __init__(self):
        self.regex = re.compile(r'^V-\(D\)-J rearrangement summary for query sequence \((.*)\)\.')
        self.query_regex = re.compile(r'lcl\|Query_\d+_reversed')  # Detects Query line
        self.query_sequence = ""
        self.ref_sequence = ""
        self.triggered = False
        self.fields = None
        self.capture_next_line = False  # Flag to capture the very next line

    def parse(self, line, out_d):
        # Parse summary fields if triggered
        if self.triggered:
            item_d = dict(zip(self.fields, line.strip().split('\t')))
            if ' Top C gene match' in item_d:
                out_d['Subtype'] = item_d[' Top C gene match'].split(',')[0].strip()
            self.triggered = False
            return out_d

        # Check for summary pattern and set fields
        match = re.match(self.regex, line)
        if match:
            self.fields = match.group(1).split(',')
            self.triggered = True
            return out_d
        # Capture nucleotides from the Query line
        return False


# SubRegionParser class to parse the sub-region sequence details
class SubRegionParser:
    """Parses out Subregion information from the 'Sub-region sequence details' section,
       specifically capturing only the CDR3 translation and its length."""
    # Dictionary of hydropathy values for each amino acid, as floats
    def __init__(self):
        # Regex to detect the 'Sub-region sequence details' header
        self.regex = re.compile(r'^Sub-region sequence details \((.*)\)')
        self.triggered = False
        self.fields = None

    def parse(self, line, out_d):
        if self.triggered:
            # Split the line by tab characters and map it to fields
            region = dict(zip(self.fields, line.strip().split('\t')))
            region_type = region['type']  # For example, 'CDR3'
            
            # Only process if the region type is 'CDR3'
            if region_type == 'CDR3':
                # We are only interested in the 'translation' field
                if ' translation' in region:
                    translation = region[' translation']
                    out_d['CDR3 AA'] = translation
                    out_d['Length'] = len(translation)  # Count the number of letters in the translation
                    # Count specific basic amino acids R, K, H
                    out_d['(+)' + 'CDR3'] = sum(1 for aa in translation if aa in {'R', 'K', 'H'})
                    # Count specific acidic amino acids D, E
                    out_d['(-)' + 'CDR3'] = sum(1 for aa in translation if aa in {'D', 'E'})
                    # Calculate difference
                    out_d['TOTAL (+/-) CDR3'] = out_d['(+)' + 'CDR3'] - out_d['(-)' + 'CDR3']
                    # Calculate GRAVY score 
                    total_hydropathy = sum(dictionary_aa.get(aa, 0) for aa in translation)
                    out_d['GRAVY CDR3'] = total_hydropathy  # Total hydropathy value
                    if len(translation) > 0:
                        out_d['GRAVY CDR3'] = total_hydropathy / len(translation)
                    else:
                        out_d['GRAVY CDR3'] = 0.0
                    start_position = int(region.get(' start', 0))  # Default to 0 if not present
                    end_position = int(region.get(' end', 0))      # Default to 0 if not present
                    out_d['CDR3 Start'] = start_position
                    out_d['CDR3 End'] = end_position
                    out_d["CDR3 Nuc"] = region.get('nucleotide sequence', 0)
                    


            # After parsing, reset the triggered state
            self.triggered = False
            return out_d

        # Check if the line matches the subregion details header
        matches = re.match(self.regex, line)
        if matches:
            # Extract the fields from the header
            self.fields = matches.group(1).split(',')
            self.fields.insert(0, 'type')  # Add 'type' for the first column (e.g., CDR3)
            self.triggered = True  # Set triggered state for parsing the next line
            return out_d

        return False

# AlignmentSummaryParser class to parse the alignment summary details
class AlignmentSummaryParser:
    def __init__(self):
        # Regex to match the header of the alignment summary section
        self.regex = re.compile(r'^Alignment summary between query and top germline V gene hit \((.*)\)')
        self.triggered = False
        self.fields = None
        self.mismatches_summary = collections.defaultdict(float)  # Store mismatches by region type
        self.gaps_summary = collections.defaultdict(float)  # Store gaps by region type
        self.positions_summary = {}  # New dictionary to store from/to positions

        self.total_mismatches = 0.0  # To accumulate total mismatches (total mutations)
        self.total_length = 0.0  # To accumulate total length for SHM calculation

    def parse(self, line, output):
        """Parses each line for alignment summary, calculates mismatches, SHM%, and stores them in a dictionary."""
        if self.triggered:
            # Create a dictionary from the fields and corresponding values in the line
            alignment = dict(zip(self.fields, line.strip().split('\t')))
            alignment_type = alignment['type']  # e.g., FR1, CDR1, CDR3 (V gene only)

            # Check if this is the Total line
            if 'Total' in alignment_type:
                self.triggered = False

                # Store mismatches by region and total mismatches in the output
                output['NCL MUT'] = dict(self.mismatches_summary)  # Convert to a regular dict for output
                output['Insertions/Deletions'] = dict(self.gaps_summary)  # Convert to a regular dict for output
                output['Positions'] = dict(self.positions_summary)  # Add positions to output
                output['Mut'] = self.total_mismatches  # Store total mismatches (total mutations)

                # Calculate SHM (somatic hypermutation percentage)
                shm = self.total_mismatches / self.total_length if self.total_length > 0 else 0.0
                output['SHM%'] = float(shm * 100)  # Add SHM percentage to the output

                # Return final output dictionary with mismatches, total mutations, SHM%, and positions
                return output

            # Extract mismatches, gaps, and length for the current region
            mismatches = float(alignment.get(' mismatches', 0))  # Default to 0 if not present
            gaps = float(alignment.get(' gaps', 0))  # Default to 0 if not present
            length = float(alignment.get(' length', 0))  # Default to 0 if not present

            # Extract from and to positions
            from_position = int(alignment.get('from', 0))  # Get from position
            to_position = int(alignment.get(' to', 0))  # Get to position
            
            # Store adjusted positions for the current alignment region
            self.positions_summary[alignment_type] = {
                'from':from_position,
                'to': to_position
            }
            # Store mismatches and gaps by region type (e.g., FR1, CDR1, etc.)
            self.mismatches_summary[alignment_type] += mismatches
            self.gaps_summary[alignment_type] += gaps
            # Accumulate total mismatches and length for SHM calculation
            self.total_mismatches += mismatches
            self.total_length += length

            return output

        # Check if the line matches the header to start processing
        matches = re.match(self.regex, line)
        if matches:
            # Extract the fields from the header
            self.fields = matches.group(1).split(',')
            self.fields.insert(0, 'type')  # Add 'type' at the start
            self.triggered = True  # Begin parsing
            return output

        return False

class Nuc:
    def __init__(self):
        self.trigger_regex = re.compile(r'^Alignments')
        self.halt_regex = re.compile(r'^Lambda')
        self.query_regex = re.compile(r'lcl\|Query_\d+_reversed')  # Detecting Query line
        self.ref_regex = None  # Set dynamically based on RefGene value
        self.query_sequence = ""
        self.ref_sequence = ""
        self.triggered = False
        self.first_query_number = None
        self.first_ref_number = None


    def reset_vars(self):
        """Reset internal variables to prepare for a new analysis."""
        self.query_sequence = ""
        self.ref_sequence = ""
        self.triggered = False
        self.first_query_number = None
        self.first_ref_number = None


    def set_refgene_regex(self, ref_gene):
        """Define regex to detect RefGene."""
        self.ref_regex = re.compile(re.escape(ref_gene)) if ref_gene else None

    def parse(self, line):
        """Parse a line to extract query and ref gene sequences."""
        # Start capturing if "Alignments" is found
        if re.match(self.trigger_regex, line):
            self.triggered = True
            return

        if not self.triggered:
            return  # Skip until the alignment block

        # Stop processing on "Lambda" line
        if re.match(self.halt_regex, line):
            return self.finish()

        # Extract query sequence nucleotides
        if re.search(self.query_regex, line):
            nucleotides = ''.join(re.findall(r'[ACGT\-]+', line))
            self.query_sequence += nucleotides
            if self.first_query_number is None:
                query_number_match = re.search(r'Query_\d+_reversed\s+(\d+)', line)
                if query_number_match:
                    self.first_query_number = query_number_match.group(1)

        # Capture ref gene sequence if line matches ref_regex
        elif self.ref_regex and re.search(self.ref_regex, line):
            dash_count = self.ref_regex.pattern.count('-')
            ref_nucleotides = ''.join(re.findall(r'[-ACGT]+', line))
            ref_nucleotides = ref_nucleotides.replace('G', '', 1)
            ref_nucleotides = ref_nucleotides.replace('-', '', dash_count)
            self.ref_sequence += ref_nucleotides
            if self.first_ref_number is None:
                # Buscamos el primer número antes de la secuencia de nucleótidos (antes de ACGT)
                ref_number_match = re.search(r'\s+(\d+)\s+[ACGT-]', line)  # Buscamos un número seguido de nucleótidos
                
                if ref_number_match:
                    self.first_ref_number = ref_number_match.group(1)  # Guardamos el primer número encontrado
    # Imprimir el primer número encontrado
    def finish(self):
        """Return query and ref sequences."""
        return {"query_sequence": self.query_sequence, "ref_sequence": self.ref_sequence, "FR1 first": self.first_query_number, "FR1 first ref": self.first_ref_number}



class AlignmentParser2():
    """Parses out Alignment summary information.
    This parser is responsible for extracting translation sequences from alignment data."""
    def __init__(self):
        self.trigger_regex = re.compile('^Alignments')
        self.halt_regex = re.compile('^Lambda')
        self.reset_vars()

    def reset_vars(self):
        """Resets the variables for a new parsing session."""
        self.triggered = False
        self.alignment_lines = []
        self.translation_string = None  # This will store the translation string

    def parse(self, line, out_d, previous_line_whitespace):
        """Parses a line of alignment data."""
        if re.match(self.trigger_regex, line):
            self.triggered = True
            return out_d

        if not self.triggered:
            return False

        if re.match(self.halt_regex, line):
            return self.finish(out_d)
        alignment_line = AlignmentLine(line)  # Assuming AlignmentLine class is defined
        self.alignment_lines.append(alignment_line)
        # Check if the line is a translation line
        if alignment_line.is_translation:  
            self.translation_string = alignment_line.al_string.strip()            
        return out_d

    def finish(self, out_d):
        # Initialize the 'Alignments' key as a list to store alignment strings
        out_d['All seq'] = []
        out_d['Indices'] = []  # Nueva clave para almacenar los índices
        previous_index = None  # Para almacenar el índice de la alineación anterior
        amino_acid_lines = []
        capturing_amino_acid = False  # Marca cuando estamos en una secuencia de aminoácidos
        amino_acid_sequence = None  # Para almacenar la secuencia de aminoácidos encontrada
        start_index = None 
        for index, alignment in enumerate(self.alignment_lines):
            stripped_alignment = alignment.al_string.strip()
            if (all(len(word) == 1 for word in stripped_alignment) and not re.search(r'[\d%]', stripped_alignment) and "-" not in stripped_alignment and ">" not in stripped_alignment):
                amino_acid_lines.append(stripped_alignment)
        filtered_amino_acid_lines = [amino_acid_lines[i] for i in range(1, len(amino_acid_lines)) if amino_acid_lines[i - 1] == ""]
        filtered_amino_acid_lines = [line.replace(" ", "") for line in filtered_amino_acid_lines if line]
        out_d['All seq']="".join(filtered_amino_acid_lines)
        # Unir todas las cadenas de alineación en una sola cadena sin espacios
        total_hydropathy = sum(dictionary_aa.get(aa, 0) for aa in out_d['All seq'])
        if len(out_d['All seq']) > 0:
            out_d['GRAVY All seq'] = total_hydropathy / len(out_d['All seq'])
        else:
            out_d['GRAVY All seq'] = 0.0
        out_d['(+)' + 'All seq'] = sum(1 for aa in out_d['All seq'] if aa in {'R', 'K', 'H'})
        # Count specific acidic amino acids D, E
        out_d['(-)' + 'All seq'] = sum(1 for aa in out_d['All seq'] if aa in {'D', 'E'})
        # Calculate difference
        out_d['TOTAL (+/-) All seq'] = out_d['(+)' + 'All seq'] - out_d['(-)' + 'All seq']
        # Restablecer el estado del analizador
        self.reset_vars()
        return out_d


class AlignmentSummaryParserProtein:
    def __init__(self):
        # Regex to match the header of the alignment summary section
        self.regex = re.compile(r'^Alignment summary between query and top germline V gene hit \((.*)\)')
        self.triggered = False
        self.fields = None
        self.mismatches_summary = collections.defaultdict(float)  # Store mismatches by region type
        self.total_mismatches_excl_cdr3 = 0.0  # To accumulate total mismatches excluding CDR3
        self.total_length = 0.0  # To accumulate total length for SHM calculation

    def parse(self, line, output):
        """Parses each line for alignment summary, calculates mismatches, and stores them in a dictionary."""
        if self.triggered:
            # Create a dictionary from the fields and corresponding values in the line
            alignment = dict(zip(self.fields, line.strip().split('\t')))
            alignment_type = alignment['type']  # e.g., FR1, CDR1, CDR3 (V gene only)
            # Check if this is the Total line
            if 'Total' in alignment_type:
                self.triggered = False
                # Store mismatches by region and total mismatches in the output
                output['AA MUT_region'] = dict(self.mismatches_summary)  # Convert to a regular dict for output
                output['AA MUT'] = self.total_mismatches_excl_cdr3  # Store total mismatches excluding CDR3
                # Calculate SHM (somatic hypermutation percentage)
                shm = self.total_mismatches_excl_cdr3 / self.total_length if self.total_length > 0 else 0.0
                output['SHM AA%'] = float(shm * 100)  # Add SHM percentage to the output
                # Return final output dictionary with mismatches, total mutations, and SHM%
                return output
            # Extract mismatches and length for the current region
            mismatches = float(alignment.get(' mismatches', 0))  # Default to 0 if not present
            length = float(alignment.get(' length', 0))  # Default to 0 if not present
            # Store mismatches by region type (e.g., FR1, CDR1, etc.)
            self.mismatches_summary[alignment_type] += mismatches
            # Accumulate total mismatches
            self.total_length += length
            # Accumulate total mismatches excluding CDR3
            self.total_mismatches_excl_cdr3 += mismatches
            return output

        # Check if the line matches the header to start processing
        matches = re.match(self.regex, line)
        if matches:
            # Extract the fields from the header
            self.fields = matches.group(1).split(',')
            self.fields.insert(0, 'type')  # Add 'type' at the start
            self.triggered = True  # Begin parsing
            return output

        return False


def translate_with_gaps(dna_sequence):
    """
    Translate a DNA sequence to a protein sequence while skipping gaps ('-').
    """
    # Initialize variables
    codon = []
    protein_seq = []

    # Iterate over the sequence, creating codons (groups of 3 nucleotides)
    for nucleotide in dna_sequence:
        if nucleotide != "-":  # Skip gaps
            codon.append(nucleotide)  # Add nucleotide to the codon list
            if len(codon) == 3:  # If we have a full codon (3 nucleotides)
                # Translate the codon into protein (as a string)
                protein_seq.append(str(Seq("".join(codon)).translate()))
                codon = []  # Reset codon for the next triplet

    # Join and return the protein sequence (it should be a string)
    return "".join(protein_seq)


def extract_sequences(output, seq):
    extracted_sequences = {}
    try:
        nucleotide_sequence = output['nucleotide_sequence'][seq]
        typeseq = seq
        positions = output['Positions']
        if 'FR1 first' in output['nucleotide_sequence']:
            fr1_first = int(output['nucleotide_sequence']['FR1 first'])
            fr1_first_ref = int(output['FR1 first ref'])
            frame = int(output['Reading Frame'])
        first_from_position = int(positions[next(iter(positions))]['from'])
        if 'CDR3 (V gene only)' not in positions and 'CDR3 Start' in output and 'CDR3 End' in output:
            positions["FR3"] = {'from': positions["FR3"]["from"], 'to': output['CDR3 Start']-1}
        
        # Normalizar la secuencia según las posiciones iniciales
        if fr1_first != first_from_position:
            nucleotide_sequence = nucleotide_sequence[((first_from_position - fr1_first)):]
            extracted_sequences['FR1 Nucleotide'] = nucleotide_sequence
    except KeyError:
        return {"Error": "Missing 'nucleotide_sequence' or 'Positions' key"}
    positions_ = [i for i, char in enumerate(nucleotide_sequence) if char == "-"]
    nucleotide_sequence = nucleotide_sequence.replace("-", "")
    # Normalizar posiciones con respecto a `first_from_position`
    normalized_positions = {}
    previous_to = None  # Mantener un seguimiento del último `to` ajustado

    if 'FR1' in positions:
        if frame == 2:
            positions['FR1']['from'] += 2
        elif frame == 3:
            positions['FR1']['from'] += 1
    for i, (region, pos) in enumerate(positions.items()):
        try:
            # Ajustamos las posiciones para que empiecen en base a `first_from_position`
            adjusted_from = pos['from'] - first_from_position
            adjusted_to = pos['to'] - first_from_position
            # Si hay un `previous_to`, ajustamos el inicio actual para evitar solapamiento
            if previous_to is not None:
                adjusted_from =  previous_to + 1
            # Ajuste si la longitud no es múltiplo de 3
            region_length = adjusted_to - adjusted_from + 1
            if region=="FR1" and region_length % 3 != 0:
                if region_length % 3 == 1:
                    adjusted_to -= 1  # Aumentamos 2 para que sea múltiplo de 3
                elif region_length % 3 == 2:
                    adjusted_to -= 2  # Aumentamos 1 para que sea múltiplo de 3
            if region!="FR1" and region_length % 3 != 0:
                if region_length % 3 == 1:
                    adjusted_to += 2  # Aumentamos 2 para que sea múltiplo de 3
                elif region_length % 3 == 2:
                    adjusted_to += 1  # Aumentamos 1 para que sea múltiplo de 3   
            # Actualizamos el `previous_to` para la siguiente región
            previous_to = adjusted_to
            # Guardamos la posición ajustada
            normalized_positions[region] = {'from': adjusted_from, 'to': adjusted_to}
        except KeyError:
            # Si falta alguna clave, la ignoramos
            continue
    if 'CDR3 AA' not in output:
        if 'CDR3 (V gene only)' not in normalized_positions:
            normalized_positions['CDR3'] = {'from': max(pos['to'] for pos in normalized_positions.values()) + 1, 'to': len(nucleotide_sequence)}
        else: 
            normalized_positions['CDR3'] = {'from': normalized_positions["CDR3 (V gene only)"]["from"], 'to': len(nucleotide_sequence)}
    # Extraer las secuencias ajustadas
    for region, pos in normalized_positions.items():
        try:
            adjusted_from = pos['from']
            adjusted_to = pos['to']
            # Extraer la secuencia nucleotídica para la región actual
            if region == 'CDR3 (V gene only)' and 'CDR3 Start' in output and 'CDR3 End' in output:
                start = normalized_positions["CDR3 (V gene only)"]["from"]
                end = normalized_positions["CDR3 (V gene only)"]["from"]+(len(output["CDR3 AA"])*3)-1
                extracted_sequences['CDR3 Nucleotide'] = nucleotide_sequence[start:end + 1]
                output["CDR3 End"]=end
                
            else:
                extracted_sequences[region + ' Nucleotide'] = nucleotide_sequence[adjusted_from:adjusted_to + 1]

        except KeyError as e:
            print(f"Error al procesar la región {region}: {e}")
    if 'CDR3 (V gene only)' not in normalized_positions and 'CDR3 Start' in output and 'CDR3 End' in output:
        output["CDR3 End"]= normalized_positions["FR3"]["to"]+len(output["CDR3 Nuc"])
  
    # Calcular FR4 si existe el final de CDR3
    if 'CDR3 End' in output:
        try:
            cdr3_end=output["CDR3 End"]+1
            if cdr3_end < len(nucleotide_sequence):
                extracted_sequences['FR4 Nucleotide'] = nucleotide_sequence[cdr3_end:]
        
        except KeyError:
            pass
    rna_protein_sequences = {}
    for region, seq in extracted_sequences.items():
        if region == "FR1 Nucleotide" and len(seq) % 3 != 0:
            seq = seq[1:]
            if len(seq) % 3 != 0:
                seq = seq[1:]
        try:
            rna_seq = seq.replace("T", "U")
            protein_seq = translate_with_gaps(seq)
            rna_protein_sequences[region.replace('Nucleotide', 'RNA')] = rna_seq
            rna_protein_sequences[region.replace('Nucleotide', 'Protein')] = protein_seq
        except Exception as e:
            continue
    extracted_sequences.update(rna_protein_sequences)
    extracted_sequences["numdif"]=len(output['nucleotide_sequence'][typeseq])-len(nucleotide_sequence[min(value['from'] for value in normalized_positions.values() if 'from' in value):])
    if nucleotide_sequence != output['nucleotide_sequence'][typeseq]:
        extracted_sequences['nucleotide_sequencenew'] = nucleotide_sequence

    fr4_protein_sequence = extracted_sequences.get('FR4 Protein', '')
    if (re.search(r'Ck|CK|CL|Cl', output["mAb name"]) and re.search(r'^[F]', fr4_protein_sequence)) or \
   (re.search(r'IgGi', output["mAb name"]) and re.search(r'^[W]', fr4_protein_sequence)):
        extracted_sequences['WGXG'] = 'YES'
    else:
        extracted_sequences['WGXG'] = 'NO'
    extracted_sequences["normalized_pos"]=normalized_positions
    
    return extracted_sequences




def compare_nucleotide_fragments(result, extracted_sequences):
    """Identifies nucleotides responsible for amino acid changes in specific fragments."""
    ns_mut_counts = {}  # Dictionary to store counts by fragment
    total_responsible_nucleotides = 0  # Total responsible nucleotides

    seq1 = result['nucleotide_sequence']['query_sequence']
    seq2 = result['nucleotide_sequence']['ref_sequence']
    endseq2 = seq2.rstrip('-')
    removed_hyphens_count = len(seq2)- len(endseq2)
    endseq1=seq1[:(len(seq1)-removed_hyphens_count)]
    seq1=endseq1
    seq2=endseq2

    if "-" not in seq1 and "-" not in seq2:
    # Translate the sequences of nucleotides to amino acids while skipping gaps
        protein1 = translate_with_gaps(seq1)
        protein2 = translate_with_gaps(seq2)

        responsible_nucleotide_positions = set()  # Unique nucleotide positions responsible

        # Compare amino acids
        for i, (aa1, aa2) in enumerate(zip(protein1, protein2)):
            if aa1 != aa2:  # If there is a change in the amino acid
                codon_start = i * 3  # Index of the first nucleotide of the codon
                codon1 = seq1[codon_start:codon_start + 3]
                codon2 = seq2[codon_start:codon_start + 3]
                if len(codon1) != 3 or len(codon2) != 3 or '-' in codon1 or '-' in codon2:
                    continue  # Skip invalid codons

                # Identify all responsible nucleotides
                for pos, (nuc1, nuc2) in enumerate(zip(codon1, codon2)):
                    if nuc1 != nuc2:
                        # Simulate substitution and verify if it affects the amino acid change
                        test_codon = list(codon1)
                        test_codon[pos] = nuc2
                        if Seq("".join(test_codon)).translate() != aa1:
                            responsible_nucleotide_positions.add(codon_start + pos)

        # Store the count of unique responsible nucleotides
        total_responsible_nucleotides += len(responsible_nucleotide_positions)
    else:
        for pos, (nuc1, nuc2) in enumerate(zip(seq1, seq2)):
            if nuc1 != nuc2:
                total_responsible_nucleotides += 1
    if total_responsible_nucleotides!=0:     
        shm_ns=total_responsible_nucleotides/len(seq1)if len(seq1) > 0 else 0
    else:
        shm_ns = 0  # Valor predeterminado si no existe o está vacío
    if 'Insertions/Deletions' in result:
        total_responsible_nucleotides=total_responsible_nucleotides-sum(result['Insertions/Deletions'].values())
    else:
        total_responsible_nucleotides=total_responsible_nucleotides
    # Return the final result
    result = {
        'NS MUT': total_responsible_nucleotides, 'SHM NS%': shm_ns
    }
    return result

# Function to parse the file using all parsers: SubRegionParser, AlignmentSummaryParser, VDJSummaryParser
def parse_file(file_path):
    query_parser = QueryParser()
    subregion_parser = SubRegionParser()
    alignment_parser = AlignmentSummaryParser()
    alignment_parser2 = AlignmentParser2()
    significant_parser = SignificantAlignmentParser()
    vdj_parser = VDJSummaryParser()
    nuc_parser = Nuc()

    output = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            query_parser.parse(line, output)
            vdj_parser.parse(line, output)
            subregion_parser.parse(line, output)
            previous_line_whitespace = line.isspace() or line == " "
            if "RefGene" in output:
                nuc_parser.set_refgene_regex(output["RefGene"])
            nuc_parser.parse(line)
            significant_parser.parse(line, output)
            alignment_parser2.parse(line, output, previous_line_whitespace)
            alignment_parser.parse(line, output)

    output['nucleotide_sequence'] = nuc_parser.finish()

    return output

        

def parse_protein_data(file_path):
    """Parse protein-specific alignment data from output1."""
    protein_parser = AlignmentSummaryParserProtein()
    output = {}
    # Intentar abrir el archivo y verificar cada línea
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            protein_parser.parse(line, output)
    return output

def parse_Cregion_ext5(output, input_folder, refgene):
    file_path = os.path.join(input_folder, output)
    Cregion5_parser = Cregion_ext5()
    nuc_parser = Nuc()
    output = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            nuc_parser.set_refgene_regex(refgene)
            Cregion5_parser.parse(line, output)
            nuc_parser.parse(line)
    output['nucleotide_sequenceext'] = nuc_parser.finish()
    return output


def write_fasta(output, fasta_folder):
    """Write output data to a FASTA file and return the original output dictionary."""
    os.makedirs(fasta_folder, exist_ok=True)
    query_name = output.get('mAb name', 'Query')
    fasta_path = os.path.join(fasta_folder, f"{query_name}.fasta")
    with open(fasta_path, 'w') as fasta_file:
        fasta_file.write(f">{query_name}\n")
        alignment_sequence = output.get('All seq', '')
        for i in range(0, len(alignment_sequence), 60):
            fasta_file.write(alignment_sequence[i:i+60] + '\n')
    return output, fasta_path  # Return the output dictionary for further use


def igblastp_execution(query_aa, seq_file):
    # Ruta base para las secuencias y otros archivos (ajusta según tu estructura)
    output_path = "C:/Users/csanchez1/Desktop/MAL67/Igblastp_output_noext"

    # Ruta base de las bases de datos germline
    germline_db_V = "C:/Users/csanchez1/Desktop/igblast-1.22.0/database3/imgt_Homo_sapiens_V_f_orf_p_prot"

    # Organismo y otros parámetros
    organism = "human"

    # Crear la carpeta de salida si no existe
    os.makedirs(output_path, exist_ok=True)

    # Bucle para procesar cada archivo de secuencia
    cmd = [
        "igblastp",
        "-germline_db_V", germline_db_V,
        "-query", query_aa,
        "-organism", organism,
        "-extend_align5end",
        "-extend_align3end",
        "-domain_system", "kabat",
        "-outfmt", "4"
    ]
    output_file = seq_file + ".txt"
    output_file_path = os.path.join(output_path, output_file)
    with open(output_file_path, 'w') as out_file:
        subprocess.run(cmd, stdout=out_file)
    return output_file_path

def compare_first_30_nucleotides(combined_output):
    """Encuentra la subsecuencia común más larga entre los primeros 30 nucleótidos de:
    'query_sequence' y 'query_sequenceext', así como entre 'ref_sequence' y 'query_sequenceext',
    y devuelve las secuencias antes de la región común en ambos casos."""

    # Obtener los primeros 30 nucleótidos de cada secuencia
    query_sequence = combined_output['nucleotide_sequence']['query_sequence'][:40]
    ref_sequence = combined_output['nucleotide_sequence']['ref_sequence'][:40]
    if "numdif" in combined_output:
        query_sequence=query_sequence[combined_output['numdif']:]
        ref_sequence=ref_sequence[combined_output['numdif']:]
    query_sequenceext = combined_output['nucleotide_sequenceext']['query_sequence'][:40]
    ref_sequenceext = combined_output['nucleotide_sequenceext']['ref_sequence'][:40]
    # Función para encontrar la secuencia sobrante antes de la región común
    def find_excess_before_common(main_seq, extended_seq):
        longest_common_subseq = ""
        common_start_index = -1
        for i in range(len(extended_seq)):
            for j in range(1, len(extended_seq) - i + 1):
                subseq = extended_seq[i:i+j]
                if subseq in main_seq:
                    if len(subseq) > len(longest_common_subseq):
                        longest_common_subseq = subseq
                        common_start_index = i

        if longest_common_subseq:
            # Extraer lo que sobra antes de la coincidencia en extended_seq
            return extended_seq[:common_start_index]
        return ""

    # Aplicar la función tanto para 'query_sequence' como para 'ref_sequence'
    before_common_in_queryext_for_query = find_excess_before_common(query_sequence, query_sequenceext)
    before_common_in_queryext_for_ref = find_excess_before_common(ref_sequence, ref_sequenceext)

    while len(before_common_in_queryext_for_query) % 3 != 0:
        before_common_in_queryext_for_query = before_common_in_queryext_for_query[1:]
        before_common_in_queryext_for_ref = before_common_in_queryext_for_ref[1:]
    min_length = min(len(before_common_in_queryext_for_query), len(before_common_in_queryext_for_ref))
    differences_count = 0
    for i in range(min_length):
        if before_common_in_queryext_for_query[i] != before_common_in_queryext_for_ref[i]:
            differences_count += 1
    # Almacenar los resultados en el diccionario 'combined_output'
    combined_output["before_common_in_queryext_for_query"] = before_common_in_queryext_for_query
    combined_output["before_common_in_queryext_for_ref"] = before_common_in_queryext_for_ref
    combined_output["Ext5 Mut"] = differences_count
    # Retornar el diccionario actualizado
    return combined_output



def parse_both_files(output_file1, fasta_folder):
    output1 = parse_file(output_file1)
    try:
        mAb_name = output1["mAb name"]
        output1, fasta_path = write_fasta(output1, fasta_folder)
        igblastp_output_path = igblastp_execution(fasta_path, mAb_name)
        output2 = parse_protein_data(igblastp_output_path)
        file_name = os.path.basename(output_file1)
        outputCregion5 = parse_Cregion_ext5(file_name, igblastCregion5, output1["RefGene"])
        combined_output = {**output1, **output2, **outputCregion5}
        combined_output["FR1 first ref"]=output1['nucleotide_sequence']["FR1 first ref"]
        combined_output["Reading Frame"] = (int(combined_output["FR1 first ref"]) - 1) % 3 + 1

        
    except KeyError:
        print("Error: 'mAb name' is missing.")
    except KeyError:
        # If "mAb name" is missing, skip igblastp execution and set output2 as an empty dictionary
        print(f"Warning: 'mAb name' is missing in {output_file1}. Skipping igblastp execution for this file.")
        combined_output = output1  # Use only output1 if "mAb name" is missing

    return combined_output



def process_multiple_files(input_folder, output_excel_path):
    combined_results = []
    # Loop over all .txt files in the input folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".txt"):
            file_path = os.path.join(input_folder, file_name)
            print(f"Processing {file_path}...")
            result = parse_both_files(file_path, fasta_folder)  # Replace with the actual protein file if needed
            extracted_sequences = extract_sequences(result, "query_sequence")
            extracted_sequences2 = extract_sequences(result, "ref_sequence")
            
            NS_mut=compare_nucleotide_fragments(result, extracted_sequences)
            result.update(extracted_sequences)
            result.update(NS_mut)
            if "CDR3 AA" not in result:
                result["CDR3 AA"] = result.get("CDR3 Protein", None)
            if result["CDR3 AA"]:
                result["Length"] = len(result["CDR3 AA"])
                result['(+)' + 'CDR3'] = sum(1 for aa in result["CDR3 AA"] if aa in {'R', 'K', 'H'})
                result['(-)' + 'CDR3'] = sum(1 for aa in result["CDR3 AA"] if aa in {'D', 'E'})
                result['TOTAL (+/-) CDR3'] = result['(+)' + 'CDR3'] - result['(-)' + 'CDR3']
                total_hydropathy = sum(dictionary_aa.get(aa, 0) for aa in result["CDR3 AA"])
                result['GRAVY CDR3'] = total_hydropathy
            else:
                result["Length"] = 0
                result['(+)' + 'CDR3'] = 0
                result['(-)' + 'CDR3'] = 0
                result['TOTAL (+/-) CDR3'] = 0
                result['GRAVY CDR3'] = 0
            if "nucleotide_sequencenew" in result:
                result['nucleotide_sequence']['query_sequence'] = result['nucleotide_sequencenew']
            nuc30 = compare_first_30_nucleotides(result)
            if "Mut" in result:
                result['TOTAL']= result['Mut']+result['Ext5 Mut']
            if 'CDR3 Start' not in result:
                result['Cut'] = "Yes"
            else:
                result['Cut'] = "No"

            combined_results.append(nuc30)
    return combined_results  # Ensure results are returned


import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

def calculate_identity(seq1, seq2):
    if not seq1 or not seq2 or len(seq1) != len(seq2):
        return 0  # Return 0 if sequences are empty or of different lengths
    
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    
    # Cargar la matriz de sustitución BLOSUM62
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')

    alignment = aligner.align(seq1, seq2)
    
    # Calcular el porcentaje de identidad basado en el mejor alineamiento
    match_count = sum(a == b for a, b in zip(alignment[0].target, alignment[0].query) if a != '-' and b != '-')
    identity_percentage = (match_count / len(seq1)) * 100
    
    return identity_percentage



def group_and_sort_by_identity(df, cdr3_column='CDR3 AA', vkl_column='V', jkl_column='J', threshold=80):
    """
    Agrupa las secuencias de acuerdo a la longitud, identidad y los valores de Vkl y Jkl.
    Las filas donde 'In frame' == 'NO' y 'Cut' == 'Yes' se les asigna el grupo 0 
    y no participan en el proceso de agrupación.
    """
    # Asegurarse de que la columna CDR3 AA no contiene valores nulos
    df = df[df[cdr3_column].notna()].copy()
    df[cdr3_column] = df[cdr3_column].astype(str)  # Convertir la columna CDR3 AA a string
    
    # Asignar grupo 0 a las filas donde 'In frame' == 'NO' y 'Cut' == 'Yes'
    df['Group'] = None
    condition = (df['In frame'] == "NO") | (df['Cut'] == "Yes")|(df['stop codon'] == "Yes")
    df.loc[condition, 'Group'] = 0
    
    # Filtrar filas que no cumplen la condición para el proceso de agrupación
    df_to_group = df[~condition].copy()
    
    group_counter = 1
    
    # Iterar sobre las filas y asignar grupos
    for i in range(len(df_to_group)):
        if pd.isnull(df_to_group.loc[df_to_group.index[i], 'Group']):  # Si no tiene grupo asignado
            df_to_group.loc[df_to_group.index[i], 'Group'] = group_counter
            seq = df_to_group.loc[df_to_group.index[i], cdr3_column]
            vkl = df_to_group.loc[df_to_group.index[i], vkl_column]
            jkl = df_to_group.loc[df_to_group.index[i], jkl_column]
            
            # Comparar con las secuencias siguientes
            for j in range(i + 1, len(df_to_group)):
                other_seq = df_to_group.loc[df_to_group.index[j], cdr3_column]
                other_vkl = df_to_group.loc[df_to_group.index[j], vkl_column]
                other_jkl = df_to_group.loc[df_to_group.index[j], jkl_column]
                
                # Verificar que las longitudes sean iguales y que Vkl y Jkl también coincidan
                if (pd.isnull(df_to_group.loc[df_to_group.index[j], 'Group']) and 
                    len(seq) == len(other_seq) and vkl == other_vkl and jkl == other_jkl):
                    if calculate_identity(seq, other_seq) >= threshold:
                        df_to_group.loc[df_to_group.index[j], 'Group'] = group_counter
            
            group_counter += 1  # Aumentar el contador de grupos
    
    # Combinar las filas agrupadas y las filas con grupo 0
    df.update(df_to_group)
    
    # Ordenar por el número de grupo para agrupar secuencias similares
    df = df.sort_values(by='Group').reset_index(drop=True)
    return df



def get_group_formats(workbook, groups_to_color):
    """
    Genera formatos de color solo para los grupos con más de una fila.
    """
    # Lista de 40 colores para asignar a los grupos
    color_codes = [
    '#FFE4E1', '#FFDAB9', '#FFFACD', '#E6E6FA', '#D3E9E4', '#D1E7DD', '#FFEFD5', '#FFF0E1',
    '#F0E68C', '#FADADD', '#E1EFFF', '#CCE5FF', '#D5F5E3', '#F8F9D7', '#FFF9E6', '#FDEBD0',
    '#E6F7FF', '#E8F0FE', '#EBF5FB', '#DFF3E3', '#E3DAC9', '#F3E9E4', '#D8D8D8', '#FFF5EB',
    '#F2E7E7', '#ECECEC', '#FAF3DD', '#FBF1E6', '#F5E6E8', '#FFF9F9', '#F0FFF4', '#FBFBFB',
    '#FFF4E6', '#FEF5E7', '#F5F5DC', '#FDF2E9', '#FCEDEB', '#F5EEF8', '#F9F3EE', '#FAF0E6'
]

    
    # Crear un diccionario de formatos de color para grupos con múltiples filas
    formats = {}
    for i, group in enumerate(groups_to_color):
        color = color_codes[i % len(color_codes)]  # Repetir colores si es necesario
        formats[group] = workbook.add_format({'bg_color': color})
    
    return formats

def filter_sequences_to_end(df):
    # Filtrar filas que deben ir al final según las condiciones en stop codon, In frame, y WGXG
    condition =  (df['In frame'] == 'NO')|(df['Cut'] == 'Yes')|(df['stop codon'] == 'Yes')  # | (df['WGXG'] == 'NO')
    filtered_df = df[condition].copy()  # Filas que cumplen la condición
    remaining_df = df[~condition].copy()  # Filas que no cumplen la condición
    return pd.concat([remaining_df, filtered_df], ignore_index=True)

# Modificar la función para aplicar borde rojo a las filas seleccionadas
def apply_color_to_letters(df, worksheet, workbook, combined):
    blue_format = workbook.add_format({'color': 'blue', 'font_name': 'Courier'})
    red_format = workbook.add_format({'color': 'red', 'font_name': 'Courier'})
    green_format = workbook.add_format({'color': 'green', 'font_name': 'Courier'})
    default_format = workbook.add_format({'color': 'black', 'font_name': 'Courier'})

    # Formato de borde rojo para las filas al final
    red_border_format = workbook.add_format({'border': 1, 'border_color': 'red'})
    if combined =="no_combined":
        # Contar los grupos y filtrar los que tienen más de una fila
        group_counts = df['Group'].value_counts()
        groups_to_color = group_counts[group_counts > 1].index.tolist()

        # Obtener formatos de fondo solo para los grupos con más de una fila
        group_formats = get_group_formats(workbook, groups_to_color)
        cdr3_counts = df['CDR3 AA'].value_counts()
        
        # Iterar filas y colorear solo los grupos seleccionados
        for row in range(1, len(df) + 1):
            group_number = df.iloc[row - 1]['Group']
            
            # Aplicar color solo si el grupo tiene más de una fila
            if group_number in group_formats:
                worksheet.set_row(row, None, group_formats[group_number])  # Colorea la fila completa

            # Colorear secuencias de CDR3 AA solo si tienen múltiples ocurrencias
            cdr3_aa_value = str(df.iloc[row - 1, df.columns.get_loc('CDR3 AA')])
            print(cdr3_aa_value)
            if cdr3_aa_value:
                formatted_text = []
                for letter in cdr3_aa_value:
                    # Aplicar formato de color según la letra
                    if letter in ['E', 'D']:
                        format = blue_format
                    elif letter in ['K', 'H', 'R']:
                        print("found")
                        format = red_format
                    else:
                        format = default_format
                    # Añadir el formato y la letra al arreglo para write_rich_string
                    formatted_text.extend([format, letter])
                # Escribir la cadena coloreada en la celda correspondiente
                worksheet.write_rich_string(row, df.columns.get_loc('CDR3 AA'), *formatted_text)
            else:
                # En caso de valor vacío, escribir sin formato especial
                worksheet.write(row, df.columns.get_loc('CDR3 AA'), cdr3_aa_value, default_format)

            # Colorear 'In frame' y 'stop codon' Colorear WGXG

            in_frame_value = df.iloc[row - 1, df.columns.get_loc('In frame')]
            stop_codon_value = df.iloc[row - 1, df.columns.get_loc('stop codon')]
            worksheet.write(row, df.columns.get_loc('In frame'), in_frame_value, green_format if in_frame_value == 'YES' else red_format)
            worksheet.write(row, df.columns.get_loc('stop codon'), stop_codon_value, red_format if stop_codon_value == 'Yes' else green_format)
            wgxg_value = df.iloc[row - 1, df.columns.get_loc('WGXG')]
            worksheet.write(row, df.columns.get_loc('WGXG'), wgxg_value, green_format if wgxg_value == 'YES' else red_format)

            # Aplicar borde rojo si la fila cumple alguna condición para estar al final
            if (df.iloc[row - 1]['In frame'] == 'NO' or  df.iloc[row - 1]['Cut'] == 'Yes' or df.iloc[row - 1]['stop codon'] == 'Yes'):      # or  df.iloc[row - 1]['stop codon'] == 'Yes') or df.iloc[row - 1]['WGXG'] == 'NO')
                worksheet.set_row(row, None, red_border_format)
          


    
    if combined=="combined":
        group_counts = df['Group'].value_counts()
        groups_to_color = group_counts[group_counts > 1].index.tolist()

        # Obtain background formats only for groups with more than one row
        group_formats = get_group_formats(workbook, groups_to_color)

        # Apply color to rows based on the group
        for row in range(1, len(df) + 1):  # Start from 1 since rows in Excel begin at 1
            row_index = row - 1  # Adjust for DataFrame's 0-based indexing
            group_number = df.iloc[row_index]['Group']  # Fetch the 'Group' value for the current row

            # Apply the group format if the group has more than one row
            if group_number in group_formats:
                worksheet.set_row(row, None, group_formats[group_number])  # Color the entire row

            # Colorear secuencias de CDR3 AA solo si tienen múltiples ocurrencias
            cdr3_aa_value_light = str(df.iloc[row - 1, df.columns.get_loc('CDR3 AA_light')])
            if cdr3_aa_value_light:
                formatted_text = []
                for letter in cdr3_aa_value_light:
                    # Aplicar formato de color según la letra
                    if letter in ['E', 'D']:
                        format = blue_format
                    elif letter in ['K', 'H', 'R']:
                        format = red_format
                    else:
                        format = default_format
                    # Añadir el formato y la letra al arreglo para write_rich_string
                    formatted_text.extend([format, letter])
                # Escribir la cadena coloreada en la celda correspondiente
                worksheet.write_rich_string(row, df.columns.get_loc('CDR3 AA_light'), *formatted_text)
                
            # Colorear secuencias de CDR3 AA solo si tienen múltiples ocurrencias
            cdr3_aa_value_heavy = str(df.iloc[row - 1, df.columns.get_loc('CDR3 AA_heavy')])
            if cdr3_aa_value_heavy:
                formatted_text = []
                for letter in cdr3_aa_value_heavy:
                    # Aplicar formato de color según la letra
                    if letter in ['E', 'D']:
                        format = blue_format
                    elif letter in ['K', 'H', 'R']:
                        format = red_format
                    else:
                        format = default_format
                    # Añadir el formato y la letra al arreglo para write_rich_string
                    formatted_text.extend([format, letter])
                # Escribir la cadena coloreada en la celda correspondiente
                worksheet.write_rich_string(row, df.columns.get_loc('CDR3 AA_heavy'), *formatted_text)
            





import sys

#with open('configpaths.json', 'r') as file:
#    config = json.load(file)

#input_folder = config['input_folder']
#fasta_folder = config['fasta_folder']
#igblastCregion5 = config['igblastCregion5']
#output_excel_path = config['output_excel_path']
#input_folder = "C:/Users/csanchez1/Desktop/MAL67/d5/Igblastn_output_noextend"
#fasta_folder = "C:/Users/csanchez1/Desktop/MAL67/d5/Sequences_fasta_aa"
#igblastCregion5 = "C:/Users/csanchez1/Desktop/MAL67/d5/Igblastn_output"
#output_excel_path = "C:/Users/csanchez1/Desktop/MAL67/d5/combined_output_d5.xlsx"
#output_excel_path = "C:/Users/csanchez1/Desktop/igblast-1.22.0/combined_output_d52.xlsx"
if len(sys.argv) != 5:
    print(f"Se esperaban 4 argumentos, pero se recibieron {len(sys.argv) - 1}.")
    print("Los argumentos esperados son: input_folder, fasta_folder, igblastCregion5, output_excel_path.")
    sys.exit(1)

input_folder = sys.argv[1]
fasta_folder = sys.argv[2]
igblastCregion5 = sys.argv[3]
output_excel_path = sys.argv[4]

#NS MUT when there are gaps, they are Mut normal because all aa are sometimes changed
if __name__ == "__main__":
    ordered_fields = [
        'mAb name', 'VK/VL/VH', 'V', 'V+allele','DH', 'D+allele','J', 'J+allele','E-value', 'Subtype', 'Strand', 'Reading Frame', 
        'CDR3 AA', 'Length', 'Mut','Ext5 Mut', 'TOTAL', 'SHM%', 'AA MUT', 'SHM AA%', 'NS MUT', 'SHM NS%', '(+)CDR3', 
        '(-)CDR3', 'TOTAL (+/-) CDR3', 'GRAVY CDR3', '(+)All seq', '(-)All seq', 
        'TOTAL (+/-) All seq', 'GRAVY All seq', 'In frame', 'stop codon', 'WGXG', 
        'FR1 Protein', 'CDR1 Protein', 'FR2 Protein', 'CDR2 Protein', 'FR3 Protein', 
        'FR4 Protein', "Cut",'All seq',  'AA MUT_region','NCL MUT', 'Insertions/Deletions'
    ]
    #input_folder = sys.argv[1]
    #fasta_folder = sys.argv[2]
    #igblastCregion5 = sys.argv[3]
    #output_excel_path = sys.argv[4]
    #input_folder = "C:/Users/csanchez1/Desktop/igblast-1.22.0/Igblastn_output_noextend"
    #fasta_folder = "C:/Users/csanchez1/Desktop/igblast-1.22.0/Sequences_fasta_aa"
    #igblastCregion5 = "C:/Users/csanchez1/Desktop/igblast-1.22.0/Igblastn_output"
   # input_folder = "C:/Users/csanchez1/Desktop/MAL67/Igblastn_output_noextend"
    #fasta_folder = "C:/Users/csanchez1/Desktop/MAL67/Sequences_fasta_aa"
    #igblastCregion5 = "C:/Users/csanchez1/Desktop/MAL67/Igblastn_output"
    #output_excel_path = "C:/Users/csanchez1/Desktop/MAL67/combined_output_d5.xlsx"
    #output_excel_path = "C:/Users/csanchez1/Desktop/igblast-1.22.0/combined_output_d5.xlsx"

    # Procesar archivos y crear DataFrame
    combined_results = process_multiple_files(input_folder, output_excel_path)
    df = pd.DataFrame(combined_results)    
    df = df.reindex(columns=ordered_fields, fill_value='')  # Llenar columnas faltantes
    # Replace NaN and inf values with empty strings in the DataFrame
    df = df.replace([float('inf'), float('-inf'), float('nan')], '')
    df = group_and_sort_by_identity(df, cdr3_column='CDR3 AA', threshold=80)
    df = filter_sequences_to_end(df)

    # Filtrar por 'mAb name' para cadenas ligeras
    light_chain_df = df[df['mAb name'].str.contains('ClXho|CLXho|Ck494|CK494', na=False)]
    light_chain_df = light_chain_df.drop(columns=['DH'], errors='ignore')  # Eliminar 'DH' si está presente
    light_chain_df = light_chain_df.drop(columns=['D+allele'], errors='ignore')  # Eliminar 'DH' si está presente

    light_chain_columns = {
        'V': 'Vkl','V+allele': 'Vkl+allele',
        'J': 'Jkl', 'J+allele': 'Jkl+allele'
        # Agregar otros cambios de nombres aquí si es necesario
    }
    light_chain_df.rename(columns=light_chain_columns, inplace=True)

    # Filtrar por 'mAb name' para IgGint
    iggint_df = df[df['mAb name'].str.contains('IgGint', na=False)]
    # Definir nombres de columnas personalizados para Heavy Chains
    heavy_chain_columns = {
        'V': 'VH', 'V+allele': 'VH+allele',
        'DH': 'DH', 'D+allele': 'DH+allele',
        'J': 'JH','J+allele': 'JH+allele',
        # Agregar otros cambios de nombres aquí si es necesario
    }
    iggint_df.rename(columns=heavy_chain_columns, inplace=True)

    # Crear un archivo Excel con múltiples hojas
    with pd.ExcelWriter(output_excel_path, engine='xlsxwriter') as writer:
        # Escribir la hoja para cadenas ligeras
        light_chain_df.to_excel(writer, index=False, sheet_name='Light Chains')
        light_chain_worksheet = writer.sheets['Light Chains']
        


        apply_color_to_letters(light_chain_df, light_chain_worksheet, writer.book, "no_combined")

        # Escribir la hoja para IgGint
        iggint_df.to_excel(writer, index=False, sheet_name='Heavy Chains')
        iggint_worksheet = writer.sheets['Heavy Chains']
        apply_color_to_letters(iggint_df, iggint_worksheet, writer.book, "no_combined")
        
        light_chain_df = light_chain_df.add_suffix('_light')
        iggint_df = iggint_df.add_suffix('_heavy')

        # Mantener `mAb name` como columna común para merge
        light_chain_df.rename(columns={'mAb name_light': 'mAb name'}, inplace=True)
        iggint_df.rename(columns={'mAb name_heavy': 'mAb name'}, inplace=True)

        # Crear clave común con el texto antes de '+'
        light_chain_df['merge_key'] = light_chain_df['mAb name'].str.split('+').str[0]
        iggint_df['merge_key'] = iggint_df['mAb name'].str.split('+').str[0]

        # Combinar en una sola hoja
        light_chain_df_filtered = light_chain_df[
            ~((light_chain_df['In frame_light'] == 'NO') |(light_chain_df['Cut_light'] == 'Yes')|(light_chain_df['stop codon_light'] == 'Yes') #(light_chain_df['stop codon_light'] == 'Yes')|(light_chain_df['WGXG_light'] == 'NO')   
            )
        ]

        iggint_df_filtered = iggint_df[
            ~((iggint_df['In frame_heavy'] == 'NO') | (iggint_df['Cut_heavy'] == 'Yes')|(iggint_df['stop codon_heavy'] == 'Yes')  #|  (iggint_df['stop codon_heavy'] == 'Yes')  |  (iggint_df['WGXG_heavy'] == 'NO')
            )
        ]

        # Ahora realiza el merge sobre los DataFrames filtrados
        merged_df = pd.merge(iggint_df_filtered, light_chain_df_filtered, on="merge_key", how="outer").drop(columns=["merge_key"])

        # Llenar valores NaN con cadena vacía
        merged_df = merged_df.fillna('')
        
    # Ordenar por `Group`, primero por Light, luego por Heavy
        
        # Escribir la hoja combinada
        merged_df_columns = {
        'mAb name_x': 'mAb name_heavy', 'mAb name_y': 'mAb name_light',
        # Agregar otros cambios de nombres aquí si es necesario
    }
        merged_df.rename(columns=merged_df_columns, inplace=True)

        merged_df['Group_heavy'] = merged_df['Group_heavy'].replace(['', ' '], np.nan)
        merged_df['Group_light'] = merged_df['Group_light'].replace(['', ' '], np.nan)

        merged_df['Group'] = merged_df.apply(
        lambda row: row['Group_heavy'] if pd.notna(row['Group_light']) and 
                    merged_df['Group_heavy'].value_counts().get(row['Group_heavy'], 0) > 1 
                    else (row['Group_light'] if pd.notna(row['Group_light']) else row['Group_heavy']),
        axis=1
    )
        group_counts = merged_df['Group'].value_counts()
        merged_df['Group_count'] = merged_df['Group'].map(group_counts)
        merged_df = merged_df.sort_values(by=['Group_count', 'Group'], ascending=[False, True])
        merged_df.drop(columns=['Group_count'], inplace=True)        
        merged_df.to_excel(writer, index=False, sheet_name='Combined Data')
        
        combined_worksheet = writer.sheets['Combined Data']
        apply_color_to_letters(merged_df, combined_worksheet, writer.book, "combined")
        
        pairs_df = merged_df[
            (merged_df["mAb name_heavy"].notna()) & 
            (merged_df["mAb name_light"].notna()) & 
            (merged_df["mAb name_heavy"] != "") & 
            (merged_df["mAb name_light"] != "")
        ]
        pairs_df.to_excel(writer, index=False, sheet_name='Pairs')
        pairs_df_worksheet = writer.sheets['Pairs']
        apply_color_to_letters(pairs_df, pairs_df_worksheet, writer.book, "combined")

        patterns_to_keep = [
            "mAb name", "VK/VL/VH", "VH", "DH", "JH","Vkl", "Jkl", "Subtype",
            "CDR3 AA", "Length", "Mut", "Ext5 Mut", "TOTAL", "SHM%", "(+)CDR3", "(-)CDR3"
        ]

        # Filtrar columnas que contienen alguno de los patrones
        filtered_columns = [
            col for col in merged_df.columns if any(pattern in col for pattern in patterns_to_keep)
        ]
        columns_to_exclude = ["TOTAL (+/-) CDR3_heavy", "TOTAL (+/-) All seq_heavy","TOTAL (+/-) CDR3_light", "TOTAL (+/-) All seq_light"]
        filtered_columns = [col for col in filtered_columns if col not in columns_to_exclude]

        # Crear un nuevo dataframe con las columnas filtradas
        filtered_df = merged_df[filtered_columns + ['Group']]
        filtered_df.to_excel(writer, index=False, sheet_name='Antibody Summary Combined')
        filtered_df_worksheet = writer.sheets['Antibody Summary Combined']
        apply_color_to_letters(filtered_df, filtered_df_worksheet, writer.book, "combined")

        
        patterns_to_keep = [
            "mAb name", "VK/VL/VH", "VH", "DH", "JH","Vkl", "Jkl", "Subtype",
            "CDR3 AA", "Length", "Mut", "Ext5 Mut", "TOTAL", "SHM%", "(+)CDR3", "(-)CDR3"
        ]

        # Filtrar columnas que contienen alguno de los patrones
        filtered_columns = [
            col for col in merged_df.columns if any(pattern in col for pattern in patterns_to_keep)
        ]
        columns_to_exclude = ["TOTAL (+/-) CDR3_heavy", "TOTAL (+/-) All seq_heavy","TOTAL (+/-) CDR3_light", "TOTAL (+/-) All seq_light"]
        filtered_columns = [col for col in filtered_columns if col not in columns_to_exclude]

        # Crear un nuevo dataframe con las columnas filtradas
        pairs_df = pairs_df[filtered_columns + ['Group']]
        pairs_df.to_excel(writer, index=False, sheet_name='Pairs Summary')
        pairs_df_worksheet = writer.sheets['Pairs Summary']
        apply_color_to_letters(pairs_df, pairs_df_worksheet, writer.book, "combined")

    print(f"Archivo guardado en {output_excel_path} con hojas para cadenas ligeras e IgGint y colores aplicados.")
