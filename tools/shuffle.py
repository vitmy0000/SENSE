import argparse
import random
import os
from itertools import islice

parser = argparse.ArgumentParser(description='Shuffle training paris')
parser.add_argument('-p', '--pair_fp', help='pair file', required=True)
parser.add_argument('-d', '--dist_fp', help='dist file', required=True)
parser.add_argument('-s', '--seed', help='random seed', required=True, type=int)
args = parser.parse_args()


def shuffle_dist(fp, seed):
    dists = []
    with open(fp) as f:
        for line in f:
            dists.append(line)
    print(len(dists))
    random.seed(seed)
    random.shuffle(dists)
    split_fp = os.path.splitext(fp)
    with open(split_fp[0] + '_shuffle' + split_fp[1], 'w') as fo:
        for x in dists:
            fo.write(x)

def shuffle_pair(fp, seed):
    pairs = []
    with open(fp) as f:
        while True:
            next_n = list(islice(f, 4))
            if not next_n:
                break
            pairs.append(''.join(next_n))
    print(len(pairs))
    random.seed(seed)
    random.shuffle(pairs)
    split_fp = os.path.splitext(fp)
    with open(split_fp[0] + '_shuffle' + split_fp[1], 'w') as fo:
        for x in pairs:
            fo.write(x)
    
if __name__ == "__main__":
    shuffle_dist(args.dist_fp, args.seed)
    shuffle_pair(args.pair_fp, args.seed)
