# ECODE-ICLR-2023
Code for [EDCDE - Extended Discovery of Closed-Form Differential Equations](https://openreview.net/forum?id=EVz_vcZQvvg&referrer=%5BAuthor%20Console%5D(%2Fgroup%3Fid%3DICLR.cc%2F2023%2FTinyPapers%2FAuthors%23your-submissions)) and an extension of the algorithm.



## Installation

Clone this repository and all submodules (e.g. using `git clone --recursive`).
Python 3.6+ is recommended. Install dependencies as per [`requirements.txt`](./requirements.txt).

## Replicating Experiments

To run the extension algorithm simply open the [`final_extension.ipynb`](./final_extension.ipynb) notebook and run all cells. One can also just run the ['run_sensitivity_vi.py'](./run_sensitivity_vi.py) script to run the sensitivity analysis for the VI algorithm. The results will be saved in the [`results`](./results) folder.

## Citing

If you use this code, please cite the original paper:

```
@inproceedings{NEURIPS2021,
  author = {Qian, Zhaozhi and Kacprzyk, Krzysztof and van der Schaar, Mihaela},
  booktitle = {International Conference on Learning Representations},
  title = {D-CODE: Discovering Closed-form ODEs from Observed Trajectories},
  url = {https://openreview.net/pdf?id=wENMvIsxNN},
  volume = {10},
  year = {2022}
}
```
