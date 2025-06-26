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
- Used as the core generative engine for the tree shown in **Videos 1,2** and **Figures 5, 7**, of the manuscript.
