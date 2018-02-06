import argparse
from itertools import islice

parser = argparse.ArgumentParser(description='Generate sequence pairs')
parser.add_argument('-i', '--input_fp', help='input fasta file', required=True)
parser.add_argument('-o', '--output_fp', help='output fasta file', required=True)
args = parser.parse_args()

if __name__ == "__main__":
    headers = []
    reads = []
    with open(args.input_fp) as f:
        while True:
            next_n = list(islice(f, 2))
            if not next_n:
                break
            headers.append(next_n[0].strip()[1:])
            reads.append(next_n[1].strip())

    with open(args.output_fp, 'w') as fo:
        for i in range(len(headers)):
            for j in range(len(headers)):
                if i < j:
                    fo.write('>{}-{}\n{}\n'.format(
                        headers[i], headers[j], reads[i]))
                    fo.write('>{}-{}\n{}\n'.format(
                        headers[i], headers[j], reads[j]))


