import argparse
import random
from itertools import islice

parser = argparse.ArgumentParser(description='Sample sub-datasets')
parser.add_argument('-i', '--input_fp', help='input fasta file', required=True)
parser.add_argument('-o', '--output_fp', help='output fasta file', required=True)
parser.add_argument('-n', '--number', help='sample size', required=True, type=int)
parser.add_argument('-s', '--seed', help='seed', default=0, type=int)
args = parser.parse_args()

if __name__ == '__main__':
    seq_id_list = []
    read_list = []
    with open(args.input_fp) as f:
        cnt = 0
        while True:
            next_n = list(islice(f, 2))
            if not next_n:
                break
            seq_id = next_n[0].strip()
            read = next_n[1].strip()
            seq_id_list.append(seq_id)
            read_list.append(read)
            cnt += 1
    random.seed(args.seed)
    select_index = sorted(random.sample(range(cnt), args.number))
    with open(args.output_fp, 'w') as fo:
        for i in select_index:
            fo.write('{}\n{}\n'.format(seq_id_list[i], read_list[i]))

