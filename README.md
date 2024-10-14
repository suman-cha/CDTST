# CDTST: General Frameworks for Conditional Two-Sample Testing

This repository contains the code for reproducing the experimental results presented in the paper "General Frameworks for Conditional Two-sample Testing" [link]. 
The paper introduces general emthodologies for conditional two-sample testing, aiming to determine whether two populations have the same distributions after accounting for confounding variables. 

## Getting Started 

### Installation

To get started with the repository, clone it to your local machine:

```sh
git clone https://github.com/suman-cha/CDTST.git
cd CDTST
```

### Running Experiments

The experiments in the paper can be reproduced using the scripts provided in the experiments/ directory. Each script corresponds to an experiment discussed in the paper.

1. **Prepare the Data**: The superconductivity dataset used in the real data experiments are included in the `real_examples/data/` folder. Make sure the data is in the correct format before running the scripts.

2. **Run Scripts**: Navigate to the `experiments/` folder and run the desired experiment script. For example:
   ```sh
   Rscript experiments/scenarios/scenario1u.R
   ```
   The results will be saved in the `results/` folder for further analysis.


## Citation

If you find this code useful for your research, please cite the paper:

TBD

## Contact

For questions or issues regarding the code, please contact the authors:

TBD 

## License

TBD

## Acknowledgments

TBD
