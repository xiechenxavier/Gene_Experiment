use io;

/* Reads instance data. */
function input() {
    local usage = "Usage: localsolver tsp.lsp "
        + "inFileName=inputFile [lsTimeLimit=timeLimit]";

    if (inFileName == nil) throw usage;
    local inFile = io.openRead(inFileName);

    // The input files follow the TSPLib "explicit" format.
    while (true) {
        local str = inFile.readln();
        if (str.startsWith("DIMENSION:")) {
            local dim = str.trim().split(":")[1];
            nbGenes = dim.toInt();
        } else if (str.startsWith("EDGE_WEIGHT_SECTION")) {
            break;
        }
    }

    // Distance from i to j
    estAvant[0..nbGenes - 1][0..nbGenes - 1] = inFile.readDouble();
}

/* Declares the optimization model. */
function model() {
    // A list variable: genes[i] is the index of the ith gene in the final solution
    genes <- list(nbGenes); 

    // All genes must be visited
    constraint count(genes) == nbGenes;

    // Maximize the total distance
    obj <- sum(0..nbGenes-2, i => sum(i+1..(nbGenes-1),j=>estAvant[genes[i]][genes[j]]));

    maximize obj;
}
/* Parameterizes the solver. */
function param() {
    if (lsTimeLimit == nil) lsTimeLimit = 100;  
}

