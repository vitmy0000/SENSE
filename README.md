# SENSE
Siamese neural network for sequence embedding


## Code Organization

This repo constains three major component:
* siamese.ipynb 
* select_training_data
* DNA_Align
* tools
* demo

`siamese.ipynb` is a notebook which contains the model definition and implementation in Pytorch.
`select_training_data` is the C++ implementation of active landmark selection algorithm for preparing the training data of SENSE.
`DNA_Align` constains the binary for evaluating the embedding results.
`tools` constains some useful python utilities.
`demo` constains data for demostration.

## Requirements

* Clang
* Cmake
* Pytorch
* CUDA

`CUDA` is need only if you need GPU acceleration, but we highly recommend using it. The installation these tools can be found on their official websites.

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

To prepare the evaluation data (this may take a long time due to sequence alignment):

```bash
python tools/sample.py -i demo/seqs.fa -o demo/eval.fa -s 0 -n 500
python tools/pair.py -i demo/eval.fa -o demo/eval_pair.fa
./DNA_Align/build/src/nw demo/eval_pair.fa demo/eval_aligned.fa
python tools/pair.py -i demo/eval_aligned.fa -o demo/eval_dist.txt
```

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
