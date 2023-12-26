import argparse
 
class DNAAnalyzer:
    def __init__(self, file_path):
        self.file_path = file_path
        self.sequences = self.read_fasta_file()
 
    def read_fasta_file(self):
        sequences = {}
        with open(self.file_path, 'r') as file:
            sequence_id = ''
            sequence = ''
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if sequence_id:
                        sequences[sequence_id] = sequence
                        sequence = ''
                    sequence_id = line[1:]
                else:
                    sequence += line
            if sequence_id and sequence:
                sequences[sequence_id] = sequence
        return sequences
 
    def find_longest_orf(self, sequence):
        start_codon = 'ATG'
        stop_codons = ['TAG', 'TGA', 'TAA']
        longest_orf = ''
        longest_orf_index = (0, 0)
 
        for i in range(len(sequence)):
            if sequence[i:i + 3] == start_codon:
                for j in range(i + 3, len(sequence), 3):
                    codon = sequence[j:j + 3]
                    if codon in stop_codons:
                        current_orf = sequence[i:j + 3]
                        current_orf_index = (i, j + 3)
                        if len(current_orf) > len(longest_orf):
                            longest_orf = current_orf
                            longest_orf_index = current_orf_index
                        break
        return longest_orf, longest_orf_index
 
    def print_result_sequences(self, print_length=False):
        result_sequences = self.find_longest_orf_in_sequences()
        for seq_id, (longest_orf, index) in result_sequences.items():
            length_info = f", Length: {len(longest_orf)}" if print_length else ""
            print(f"Sequence ID: {seq_id}, ORF: {longest_orf}{length_info}, Index: {index}")
 
    def print_result_overall(self, print_length=False):
        result_overall, overall_index = self.find_overall_longest_orf()
        length_info = f", Length: {len(result_overall)}" if print_length else ""
        print(f"Overall longest ORF: {result_overall}{length_info}, Index: {overall_index}")
 
    def find_longest_orf_in_sequences(self):
        longest_orfs = {}
        for seq_id, sequence in self.sequences.items():
            longest_orf, index = self.find_longest_orf(sequence)
            longest_orfs[seq_id] = (longest_orf, index)
        return longest_orfs
 
    def find_overall_longest_orf(self):
        longest_orfs = self.find_longest_orf_in_sequences()
        overall_longest_orf, overall_index = max(longest_orfs.values(), key=lambda x: len(x[0]), default=('', (0, 0)))
        return overall_longest_orf, overall_index
 
 
def main():
    parser = argparse.ArgumentParser(description='Find ORFs in DNA sequences from a FASTA file.')
    parser.add_argument('file_path', help='Path to the FASTA file')
    parser.add_argument('--print_length', action='store_true', help='Print the length of the found ORFs')
    args = parser.parse_args()
 
    dna_analyzer = DNAAnalyzer(args.file_path)
    dna_analyzer.print_result_sequences(print_length=args.print_length)
    dna_analyzer.print_result_overall(print_length=args.print_length)
 
 
if __name__ == "__main__":
    main()