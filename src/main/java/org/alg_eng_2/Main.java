package org.alg_eng_2;

import static org.alg_eng_2.PARAMETERS.DEBUG_MODE;

public class Main {

    public static void printGraphInfo(DynGraph g) {
        System.err.printf("number of vertices: %d%n", g.numVertices());
        System.err.printf("number of edges: %d%n", g.numEdges());
    }

    public static void main(String[] args) {
        if (args.length > 1) {
            System.err.println("Pass at most one argument: a filepath or \"-\" for stdin.");
            return;
        }
        String filepath = args.length == 1 ? args[0] : "-";
        if (DEBUG_MODE) {
            System.err.printf("filepath: %s%n", filepath);
        }

        DynGraph g;
        try {
            g = new GraphBuilder(filepath).build();
        } catch (Exception ignored) {
            System.err.println("Something went wrong while trying to parse the graph.");
            return;
        }


        Solver.solveAndPrint(g);
    }
}
