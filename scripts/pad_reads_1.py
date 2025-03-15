from Bio import SeqIO
from Bio.Seq import Seq
import sys

def pad_read(read, length):
    if len(read) < length:
        return read + 'N' * (length - len(read))
    else:
        return read[:length]

def pad_fasta(input_file, output_file, length):
    with open(output_file, 'w') as outfile:
        for record in SeqIO.parse(input_file, 'fasta'):
            padded_seq = pad_read(str(record.seq), length)
            record.seq = Seq(padded_seq)  # Convert string back to Seq object
            SeqIO.write(record, outfile, 'fasta')

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    length = int(sys.argv[3])
    pad_fasta(input_file, output_file, length)
