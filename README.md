# GECCO_2023

This directory contains the Tabu Search and Variable Function Search algorithms used in the experimentation of "New knowledge about the Elementary Landscape Decomposition for solving the Quadratic Assignment Problem". In order to execute the algorithms, please follow the following steps:

- Enter the corresponding directory (TabuSearch or VariableFunctionSearch).

- Compile the code using the "make" command.

- Execute the algorithms using the following commands:

    - Tabu Search: ./TS [Instance] [Random seed] [Tabu list size] [Number of solution evaluations] [Output header (yes=1)]

    - Variable Function Search: ./VFS [Instance] [Random seed] [Tabu list size] [Number of solution evaluations] [Output header (yes=1)]

- The algorithms output information about the solutions visited at each step of the search process (CSV format).
    
As an example, the results of the TS and VFS executions that were considered for the experimental analysis of Section 6 can be found in the Data directory.
