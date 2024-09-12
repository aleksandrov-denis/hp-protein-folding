# Objective
Solves the HP model protein structure folding problem using a custom Monte Carlo approach.

# Responsibilities
Responsible for implementing perturbations of given protein structures, developing a Monte-Carlo-based minimization algorithm for protein folding, and reporting improvements relative to the brute force method.

# Output
Outputs 10 iterations of the Monte-Carlo approach to solvoing the HP model protein folding problem. For each iteration, lists the generation at which the optimal score was found, and the structure and score associated with said generation.

# To Run
With jfreechart installed and appropriate classpath ```java -cp PATH Protein SEQUENCE``` where PATH is the classpath for the jfreechart.jar and SEQUENCE is something like HPHPHHPHPPHP, for example. This option shows a chart of the fold score for each Monte-Carlo generation.

Without jfreechart installed, in case you don't want a chart shown, comment out the imports relating to jfreechart and the graph method. Run ```java Protein SEQUENCE```.
