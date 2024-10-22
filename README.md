# Cond2ST: General Frameworks for Conditional Two-Sample Testing

This repository contains the R codes for reproducing the experimental results presented in the paper "General Frameworks for Conditional Two-sample Testing" [link]. 
The paper introduces general methodologies for conditional two-sample testing, aiming to determine whether two populations have the same distributions after accounting for confounding variables. We introduce two approaches: (1) converting conditional independence tests to conditional two-sample tests and (2) density ratio-based testing. If you are planning to do tests for your own projects, we recommend to use R package named "Cond2ST". 

## Getting Started 

### Installation

To get started with the repository, clone it to your local machine:

```sh
git clone https://github.com/suman-cha/Cond2ST.git
cd Cond2ST
```

### Requirements

Make sure you have R installed (version 4.0 or higher) along with the necessary packages.

### Run Experiments 
The experiments in the paper can be reproduced using the R codes provided in the experiments/ directory. Each code corresponds to an experiment discussed in the paper.

1. **Prepare the Data**: The superconductivity dataset used in the real data experiments are included in the `real_examples/data/` folder. Or you can download data yourself from [superconductivity](https://archive.ics.uci.edu/dataset/464/superconductivty+data). Make sure the data is in the correct format before running the scripts. Other simulation data are generated in each code. 

2. **Run Scripts**: Navigate to the `experiments/` folder and run the desired experiment script. For example:
   ```sh
   Rscript experiments/scenarios/scenario1u.R
   ```
   The results will be saved in the `results/` folder.



## Citation

If you find this code useful for your research, please cite the paper:

TBD

## Contact

For questions or issues regarding the code, please do not hesitate to contact us!

Seongchan Lee* : statchan1106@yonsei.ac.kr

Suman Cha* : oldrain123@yonsei.ac.kr 

Ilmun Kim : ilmun@yonsei.ac.kr 

## License

TBD

## Acknowledgments

We acknowledge support from the Basic Science Research Program through the National Research Foundation of Korea (NRF) funded by the Ministry of Education (2022R1A4A1033384), and the Korea government (MSIT) RS-2023-00211073.
