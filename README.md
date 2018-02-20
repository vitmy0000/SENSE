# SENSE
Siamese neural network for sequence embedding

## Code Organization

This repo constains three major components:
* siamese.ipynb 
* select_training_data
* DNA_Align
* tools
* demo

`siamese.ipynb` is a notebook that contains the model definition and implementation in Pytorch.
`select_training_data` is the C++ implementation of the active landmark selection algorithm for preparing training data for SENSE.
`DNA_Align` constains the binary for evaluating the embedding results.
`tools` constains some useful python utilities.
`demo` constains data for demonstration.


## Requirements

* Clang
* Cmake
* Boost
* Pytorch
* CUDA

`CUDA` is needed only if you need GPU acceleration, but we highly recommend using it. The installation of these tools can be found on their official websites.

## Compile

To compile /DNA_Align/:
```bash
cd DNA_Align
mkdir build && cd build
cmake .. && make
```
The executable binary should be under `DNA_Align/build/src/`.


To compile /select_training_data/:
```bash
cd select_training_data
mkdir build && cd build
cmake .. && make
```
The executable binary should be under `select_training_data/build/src/`.

## Demo

The demo dataset contains 10,000 sequences sampled from RT988 dataset.
In this demo, we sampled 500 out of 10,000 sequences and compute their pairwise distances for evaluation.
This process may take a long time due to sequence alignment.
To prepare the evaluation data:

```bash
python tools/sample.py -i demo/seqs.fa -o demo/eval.fa -s 0 -n 500
python tools/pair.py -i demo/eval.fa -o demo/eval_pair.fa
./DNA_Align/build/src/nw demo/eval_pair.fa demo/eval_aligned.fa
python tools/dist.py -i demo/eval_aligned.fa -o demo/eval_dist.txt
```

In this demo, we prepare 20 * 500 training sequence paris and shuffle them.
To select training data:
```bash
./select_training_data/build/src/select_training_data -f demo/seqs.fa -s demo/seqs_ids.txt -p demo/pair.fa -d demo/dist.txt -a 1 -t 20 -n 500
python tools/shuffle.py -p demo/pair.fa -d demo/dist.txt -s 0
```

Here is the help for the options:
```cpp
options.add_options()
  ("f,fasta_file", "input fasta file", cxxopts::value<std::string>())
  ("s,seq_ids_file", "output seq ids file", cxxopts::value<std::string>())
  ("p,pairs_file", "output pairs file", cxxopts::value<std::string>())
  ("d,dists_file", "output dists file", cxxopts::value<std::string>())
  ("a,abundance_threshold", "abundance threshold", cxxopts::value<std::size_t>())
  ("t,target_num_landmarks", "target number of landmarks", cxxopts::value<std::size_t>())
  ("n,num_random_sample", "number of random sample", cxxopts::value<std::size_t>())
```

Run the jupyter notebook for defining, training and evaluating the model.
