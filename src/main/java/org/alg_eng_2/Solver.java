package org.alg_eng_2;

import it.unimi.dsi.fastutil.ints.IntIntImmutablePair;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;

import java.util.*;

import static org.alg_eng_2.Main.printGraphInfo;
import static org.alg_eng_2.PARAMETERS.*;

public class Solver {
    public static void solveAndPrint(DynGraph g) {

        IntOpenHashSet reduceVertices012 = g.reduceVertices012();

        if (g.numVertices() == 0) {
            reduceVertices012.forEach(System.out::println);
            return;
        }

        ArrayList<Integer> maxClique = new ArrayList<>();
        int maxDegeneracy = -1; // disable the Turing kernel
        DynGraph[] instances = TuringKernel.applyTo(g, maxDegeneracy);
        assert instances.length == g.numVertices() || instances.length == 1;
        // compute the maximum order of an instance
        ArrayList<Integer> orders = new ArrayList<>();
        orders.add(0); // to handle the null graph
        for (DynGraph h: instances) {
            orders.add(h.numVertices());
        }
        int maxOrder = Collections.max(orders);
        if (PRINT_IF_KERNEL_WAS_APPLIED) {
            boolean kernelApplied = instances.length == g.numVertices();
            System.err.printf("kernel applied: %b", kernelApplied);
            if (kernelApplied) {
                int degeneracy = maxOrder >= 1 ? maxOrder - 1 : 0;
                System.err.printf(", degeneracy: %d", degeneracy);
            } else {
                System.err.printf(", degeneracy: > %d", maxDegeneracy);
            }
        }
        for (DynGraph h: instances) {
            if (h.numVertices() <= maxClique.size()) {
                // no need to solve
                continue;
            }
            ArrayList<Integer> hMaxClique = solve(h);
            if (hMaxClique.size() > maxClique.size()) {
                maxClique = hMaxClique;
                if (maxClique.size() == maxOrder) {
                    // lower bound and upper bound coincide; we are done
                    break;
                }
            }
        }
        if (maxClique.size() < reduceVertices012.size()) {
            maxClique = new ArrayList<>(reduceVertices012);
        }
        maxClique.forEach(System.out::println);
    }

    public static IntOpenHashSet sortAndReduce(DynGraph g, boolean firstSort) {

        IntOpenHashSet reduceVertices012 = g.reduceVertices012();

        if (firstSort) {
            g.sortDegree();
        }

        long reduceStartTime = System.currentTimeMillis();

        int[] vertices = g.vertexSetCopy().toArray(new int[0]);
        Arrays.sort(vertices);
        IntOpenHashSet exclude = new IntOpenHashSet();
        IntOpenHashSet bestClique = new IntOpenHashSet();

        for (int i = vertices.length-1; i >= 0; i--) { // best way after sorting by degree
            if (System.currentTimeMillis() - reduceStartTime > MAX_REDUCE_TIME) {
                break;
            }
            int v = vertices[i];
            if (DEBUG_MODE) {
                System.err.printf("vertex: %d, minSize: %d\n",v, bestClique.size());
            }
            if (!g.hasVertex(v)) {
                continue;
            }
            if (exclude.contains(v)) {
                continue;
            }
            if (g.degree(v) <= bestClique.size()) {
                continue;
            }
            IntOpenHashSet clique;
            if (g.numVertices() > MIN_VERTICES_DYN) {
                clique = g.dynFastClique(v,bestClique.size());
            } else {
                clique = g.fastClique(v,bestClique.size());
            }
            if (clique.size() > bestClique.size()) {
                bestClique.clear();
                for (int c : clique) {
                    bestClique.add(g.to_original(c));
                }
                if (g.numVertices() > MIN_VERTICES_DYN) {
                    g.dynReduceVerticesAndEdges(bestClique.size());
                } else {
                    g.reduceVerticesAndEdges(bestClique.size());
                }
            }
            exclude.addAll(clique);
        }

        g.sortDegree();

        if (bestClique.size() > reduceVertices012.size()) {
            return bestClique;
        }

        return reduceVertices012;
    }

    public static ArrayList<Integer> solve(DynGraph g) {
        ArrayList<Integer> result = new ArrayList<>(); // contains a maximum clique when done

        if (DEBUG_MODE) {
            printGraphInfo(g);
        }

        IntOpenHashSet bestClique = sortAndReduce(g,true);

        if (DEBUG_MODE) {
            System.err.println(g.checkOriginalMappings());
            System.err.println(g.isDegreeSorted());
            System.err.println(bestClique);
            System.err.printf("Best Clique Size: %d\n",bestClique.size());
            System.err.println(g.isClique(bestClique));
            System.err.printf("Num Vertices: %d\n",g.numVertices());
            System.err.printf("Num Edges: %d\n",g.numEdges());
        }

        if (g.numVertices() == 0) {
            //System.err.printf("%d,%d", bestClique.size(),bestClique.size());
            return new ArrayList<>(bestClique);
        }

        StaticGraph staticGraph = new StaticGraph(g);

        int[] max_clique = staticGraph.maxClique();
        for (int i = 0; i < staticGraph.numVertices(); i++) {
            if (max_clique[i] == -1) {
                if (DEBUG_MODE) {
                    System.err.println(i);
                }
                break;
            }
            result.add(g.to_original(staticGraph.to_original(max_clique[i])));
        }

        if (bestClique.size() > result.size()) {
            return new ArrayList<>(bestClique);
        }

        return result;
    }
}
