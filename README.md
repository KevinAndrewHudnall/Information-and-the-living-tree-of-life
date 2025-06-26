# Information-and-the-living-tree-of-life
Contains all the program and data files for the paper "Information and the living tree of life: A theory of time and measurement grounded in biology"
## Data Availability
The full dataset (≈2.5 GB) used in this study is archived on Zenodo:
10.5281/zenodo.15749593
This includes all data required to reproduce the figures, simulations, and analyses described in the manuscript *Information and the Living Tree of Life*.
### `BuildMultifractalTreeFn.m`
Implements a **random iterated function system (RIFS)** to generate a multifractal branching tree.  
- Returns three matrices:
  - `S` – scale values across generations and lineages  
  - `P` – number of progeny per node, preserving ancestral structure  
  - `H` – entropy estimated as the log-change in scale from generation to generation  
- At each iteration, the function draws a random scale and spawns recursive subtrees using a random branching process (up to `MaxOffspring` and `MaxGens`), resulting in a deeply nested multifractal structure.  
- The matrices explicitly reconstruct the full history of the tree, making them suitable for entropy analysis, observer-relative measurements, and figure generation.  
- Used as the core generative engine for the tree shown in **Videos 1, 2** and **Figures 5, 7**, of the manuscript.
  
### `CalculateFractalDimsFn.m`

Estimates the **fractal (correlation) dimension** of each lineage path in a multifractal tree.  
- Uses the epsilon-neighborhood method, counting point pairs closer than ε across generations.  
- Applies log–log linear regression to compute the correlation dimension \( D_F \) and goodness-of-fit (R²) for each path.  
- Outputs:
  - `D_F`: matrix of correlation dimensions across all paths and generations  
  - `R_Squared`: R² values indicating fit quality  
  - `Avg_R2_Vector`, `Std_Devs`: generation-wise statistics on fit quality  
  - A diagnostic analysis of a randomly selected path, including its dimension and correlation curve  
- Supports evaluation of scaling behavior and reliability across the multifractal structure  
- Used in quantifying scale-dependent complexity and verifying that pathwise growth exhibits fractal scaling

### `GenerateFiguresFn.m`

Generates all figures presented in the manuscript *Information and the Living Tree of Life* using precomputed data.  
- Produces plots of nested form and entropy dynamics (Figure 5 and 7)  
- Displays fractal dimension scaling across paths and pairwise comparisons (Figure 6 and D1)  
- Visualizes pairwise joint entropy, mutual information, and information distance (Figures 8–9)  
- Plots real and complex solutions to the dilation equation (Figure 10)  
- Depicts observer-relative rates of time elapse (Figure 12)  
- Renders MF-DFA results across q moments (Figure E1)  
- Assumes data variables (e.g., `Randomly_Reduced_S`, `Crosspath_H_Reduced`, `Dil_Eq_Real_Sols`) are available in the workspace  

This script serves as the visual output module for the entire analysis pipeline, ensuring full reproducibility of the figures used in publication.

### `GetCrosspathQuantities.m`

Computes all pairwise information-theoretic quantities and dilation equation solutions between leaves of a multifractal tree — one chunk at a time.  
- Used to parallelize and scale the computation of crosspath relationships between leaves.  
- For each chunk (indexed by `k`), returns a structured output (`Data_Chunk`) with:
  - Locations and values of most recent common ancestors (MRCA)
  - Pairwise entropy, mutual information, and information distance (H, I, d)
  - Crosspath fractal dimension estimates at MRCAs
  - Solutions to the dilation equation and their classification as real or complex
  - Counts of coherent (positive), divergent (negative), imaginary, and indeterminant comparisons  
- Filters out unresolved comparisons using a user-provided `Conv_Tol` matrix  
- Core function behind Figures 8, 9, and 10 in the manuscript, enabling scalable evaluation of pairwise information structure and observer-relative time dynamics

