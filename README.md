# Longest Common Substring in Longest Common Subsequence's Solution Service: A Novel Hyper-Heuristic

Authors: Alireza Abdi, Masih Hajsaeedi, Mohsen Hooshmand

This code is written with python 3.8 for the paper with the above title. This version returns the LCS's length, not the LCS. 

The dataset contains six benchmarks, ACO-Rat, ACO-Random, ACO-Virus, ES, BB, and SARS-CoV-2.
These datasets are introduced in their original papers, and we provide them here.

To run the code in windows, do the following instruction:

1. extract datasets and code in the same folder
2. open cmd in the corresponding folder
3. type the following command in cmd: python UB-HH.py "dataset name" "beam width"

An example: python UB-HH.py ACO-Rat 200 
  
Note: When writing the name of the dataset, make sure that it exactly matches the above names.

Pleas, feel free to contact us if you have any questions about the paper or implementation.

The GMPSUM method (our competitor) code that we implemented is also available here, and its running instruction is like the UB-HH. Except for the lambda parameter that the user should specify.

An example: python GMPSUM.py ES 200 1.0

Email: alirezaabdi@iasbs.ac.ir

### citation

```bash
@article{abdi2023longest,
  title={Longest common substring in Longest Common Subsequence’s solution service: A novel hyper-heuristic},
  author={Abdi, Alireza and Hajsaeedi, Masih and Hooshmand, Mohsen},
  journal={Computational Biology and Chemistry},
  pages={107882},
  year={2023},
  publisher={Elsevier}
}
```
