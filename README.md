# Information-and-the-Living-Tree-of-Life

This repository contains all program and data files associated with the manuscript:  
**_Information and the Living Tree of Life: A Theory of measurement grounded in biology_**.

The code implements a stochastic, information-theoretic framework for understanding evolutionary dynamics, multifractal geometry, and observer-relative time.

> The code for generating animations can be found in a separate repository:  
> https://github.com/KevinAndrewHudnall/the-living-tree-of-life

---

## Data Availability

The full dataset (~2.5 GB) used in this study is archived on Zenodo:  
**https://doi.org/10.5281/zenodo.15749593**

This archive includes all data required to reproduce the manuscript's simulations, figures, and statistical analyses.

---

## Main Execution Pipeline

### `Main_Script.m`

**Warning: This script is computationally intensive and designed for high-memory machines with multiple CPU cores.**  
Reducing the `ITERATIONS` parameter lowers computational cost.

This is the primary driver for:

- Constructing a multifractal branching tree (`BuildMultifractalTreeFn`)
- Computing correlation dimensions of all paths (`CalculateFractalDimsFn`)
- Identifying valid pairwise comparisons based on convergence
- Computing pairwise quantities using `parfeval`, including:
  - Joint entropy \( H \)
  - Mutual information \( I \)
  - Information distance \( d \)
  - Fractal dimension at the most recent common ancestor (MRCA)
  - Real and complex solutions to the dilation equation

Results are saved chunk-by-chunk and later aggregated.

---

## Core Functions

### `BuildMultifractalTreeFn.m`

Implements a random iterated function system (RIFS) to generate a multifractal evolutionary tree.

- Outputs:
  - `S` – scale matrix (nodes × lineages)
  - `P` – offspring count matrix
  - `H` – entropy from intergenerational scale transitions
- Recursively spawns subtrees using a random branching and scaling process
- Basis for Figures 5, 7 and Videos 1–2

---

### `CalculateFractalDimsFn.m`

Estimates the fractal (correlation) dimension of each lineage using the ε-neighborhood method.

- Computes log–log regression of neighbor count vs. ε
- Outputs:
  - `D_F` – dimension matrix
  - `R_Squared` – fit quality
  - Summary statistics across generations
- Supports filtering by convergence and scaling diagnostics (used in Figure 6)

---

### `MfDfaFn.m`

Performs Multifractal Detrended Fluctuation Analysis (MF-DFA) on lineage paths.

- Inputs: `S`, `q_Values`, `Box_Sizes`
- Outputs:
  - `Generalized_Hurst_values`: matrix of H(q) exponents per path and moment q
- Used in Figure E1 to characterize multifractal scaling of time along biological lineages

---

### `GetCrosspathQuantities.m`

Processes one chunk of pairwise comparisons at a time to compute crosspath quantities.

- Inputs: `S`, `H`, `Conv_Tol`, and chunk index `k`
- Computes for each pair:
  - MRCA location and value
  - Joint entropy, mutual information, information distance
  - Fractal dimension at MRCA
  - Real/complex solution to the dilation equation
- Filters unresolved or indeterminant comparisons
- Used to generate Figures 8, 9, and 10

---

### `MakeRandomTree.m`

Generates a Galton–Watson random tree using a randomly sampled offspring distribution.

- Inputs:
  - `Max_Offspring`: maximum number of offspring per node
  - `Max_Gens`: maximum number of generations
- Output:
  - `Tree`: number of terminal nodes (leaves) in the resulting tree
- Used by `BuildMultifractalTreeFn` to instantiate the stochastic backbone of each multifractal subtree

Note: This function returns only the number of leaves. For tree structure visualization, see `MakeRandomTreeForVisual.m` in:  
[https://github.com/KevinAndrewHudnall/the-living-tree-of-life/tree/main/Functions](https://github.com/KevinAndrewHudnall/the-living-tree-of-life/tree/main/Functions)

---

## Statistical Validation

### `Statistical_Confirmation.m`

**Warning: This script is computationally expensive and should be run on a workstation.**

Validates that the proportions of coherent, divergent, convergent, unresolved, and indeterminant outcomes are statistically representative:

- Repeats the entire simulation and analysis pipeline multiple times
- For each run:
  - Constructs a multifractal tree
  - Samples a subset of paths (`n`)
  - Computes \( \binom{n}{2} \) pairwise comparisons
  - Calculates H, I, d, MRCA, and dilation solution
  - Tallies result types (coherent, divergent, real, imaginary, etc.)
- Outputs mean and standard deviation across runs
- Basis for Figures 10–12

---

## Visualization and Analysis

### `DataPlots.m`

Generates all manuscript figures from precomputed data:

- Nested structure and entropy dynamics (Figures 5, 7)
- Pathwise fractal dimension scaling (Figures 6, D1)
- Information quantities and dilation equation solutions (Figures 8–10)
- Observer-relative time rate distributions (Figure 12)
- MF-DFA analysis (Figure E1)

Assumes relevant variables are loaded into the MATLAB workspace.

---

## Summary Table

| Script / Function             | Purpose                                                    |
|------------------------------|-------------------------------------------------------------|
| `Main_Script.m`              | Runs the full simulation and analysis pipeline              |
| `BuildMultifractalTreeFn.m` | Generates multifractal evolutionary trees                   |
| `CalculateFractalDimsFn.m`  | Computes fractal dimension of paths                         |
| `MfDfaFn.m`                 | Computes Hurst exponents via MF-DFA                         |
| `GetCrosspathQuantities.m`  | Computes pairwise information & dilation results            |
| `MakeRandomTree.m`          | Samples a Galton–Watson branching tree                      |
| `Statistical_Confirmation.m`| Validates representativeness via simulations                |
| `DataPlots.m`               | Generates all manuscript figures                            |
