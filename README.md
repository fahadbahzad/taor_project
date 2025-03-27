# Topics in Applied Operational Research - Group 23

This repository contains implementations of two algorithms: 
- **Prize-Collecting Travelling Salesman Problem (PC-TSP)** with solver callbacks.
- **Recursive Algorithm (RA)** with an optional stochastic element.

Additionally, the repository includes example files to guide users in running the algorithms, dataset generation, and timing analysis.

## Repository Structure

### Modules
- **`data.py`**: Handles reading and formatting the datasets, ensuring compatibility with the specifications required by both algorithms.
- **`pctsp.py`**: Implements the PC-TSP approach utilizing solver callbacks for efficient computation.
- **`ra.py`**: Implements the Recursive Algorithm including an optional stochastic approach.

### Example Notebooks
- **`pctsp_example.ipynb`**: Demonstrates how to run the PC-TSP algorithm module.
- **`ra_example.ipynb`**: Demonstrates how to run the Recursive Algorithm module.
- **`timing_results.ipynb`**: Provides the code for reproducing timing results from solving the 25 datasets in the `data/` folder.
- **`difficult_instance_example.ipynb`**: Showcases an example of a difficult instance and compares the performance of PC-TSP and RA, highlighting how PC-TSP is faster.
- **`stochastic_results_generation.ipynb`**: Runs the stochastic Recursive Algorithm on Dataset40.xlsx and Dataset50.xlsx with various parameters for testing.
- **`stochastic_results_assessment.ipynb`**: From stochastic_results.csv to evaluate the performance of the stochastic Recursive Algorithm.

### Datasets
- **`data/`**: This folder contains the datasets generated for testing and evaluation purposes.

## Getting Started

### Prerequisites
- Python 3.x
- Required libraries (see `requirements.txt`)

### Running the Algorithms
1. **PC-TSP**: Open and execute `pctsp_example.ipynb` for a guided walkthrough on running the PC-TSP algorithm.
2. **RA**: Open and execute `ra_example.ipynb` for a guided walkthrough on running the Recursive Algorithm.

### Reproducing Results
- Run the `timing_results.ipynb` notebook to replicate the timing analysis for the 25 datasets.
- Open `difficult_instance_example.ipynb` to explore a complex instance and observe performance differences between the algorithms.

---

Enjoy exploring the PC-TSP and RA algorithms!