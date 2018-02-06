import argparse
from itertools import islice

parser = argparse.ArgumentParser(description='Calculate alignment distance')
parser.add_argument('-i', '--input_fp', help='input file', required=True)
parser.add_argument('-o', '--output_fp', help='output file', required=True)
args = parser.parse_args()

def distance(str_x, str_y):
    assert(len(str_x) == len(str_y))
    match_cnt = 0
    miss_cnt = 0
    for char_x, char_y in zip(str_x, str_y):
        if char_x == char_y:
            match_cnt += 1
        else:
            miss_cnt += 1
    return float(miss_cnt) / (miss_cnt + match_cnt)

def get_pdist(input_fp, output_fp):
    with open(input_fp) as f, open(output_fp, 'w') as fo:
        while True:
            next_n = list(islice(f, 4))
            if not next_n:
                break
            fo.write('{}\t{:.4f}\n'.format(
                next_n[0].strip(),
                distance(next_n[1].strip(), next_n[2].strip())))

if __name__ == "__main__":
    get_pdist(args.input_fp, args.output_fp)
