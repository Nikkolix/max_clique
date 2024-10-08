package org.alg_eng_2;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;

import java.util.*;

public class TuringKernel {
    public static int[] degeneracyOrderingSlow(DynGraph g) {
        int[] ordering = new int[g.numVertices()];

        IntOpenHashSet remaining = g.vertexSetCopy();
        DynGraph clone = g.copy();
        int i = 0;
        while (i < g.numVertices()) {

            int minDeg = Integer.MAX_VALUE;
            int minDegV = -1;
            for (int u : remaining) {
                int degU = clone.degree(u);
                if (degU < minDeg) {
                    minDegV = u;
                    minDeg = degU;
                }
            }

            remaining.remove(minDegV);
            ordering[i] = minDegV;
            clone.removeVertex(minDegV);
            i++;
        }

        return ordering;
    }

    public static int[] degeneracyOrdering(DynGraph g, int maxDegeneracy) {
        if (maxDegeneracy < 0) {
            return null;
        }
        DynGraph h = g.copy();
        int n = h.numVertices();
        int[] ordering = new int[n];
        // construct mapping from degree to set of vertices with that degree
        Int2ObjectOpenHashMap<IntOpenHashSet> degreeToVertices = new Int2ObjectOpenHashMap<>();
        IntOpenHashSet degrees = new IntOpenHashSet();
        degreeToVertices.put(0, new IntOpenHashSet()); // to gracefully handle the null graph
        for (int v: h.vertexSetCopy()) {
            int degree = h.degree(v);
            if (!degreeToVertices.containsKey(degree)) {
                degreeToVertices.put(degree, new IntOpenHashSet());
                degrees.add(degree);
            }
            degreeToVertices.get(degree).add(v);
        }
        // turn mapping into an array list
        int maxDegree = Collections.max(degrees);
        ArrayList<IntOpenHashSet> degreeToVerticesFast = new ArrayList<>(maxDegree + 1);
        for (int d = 0; d <= maxDegree; d++) {
            IntOpenHashSet vertices = degreeToVertices.containsKey(d) ? degreeToVertices.get(d) : new IntOpenHashSet();
            degreeToVerticesFast.add(vertices);
        }
        // compute the degeneracy ordering
        for (int i = 0; i < n; i++) {
            // find a vertex of minimum degree
            int minDegreeVertex = -1;
            int minDegree = -1;
            for (int d = 0; d <= maxDegree; d++) {
                if (d > maxDegeneracy) {
                    // degeneracy is too high
                    return null;
                }
                IntOpenHashSet vertices = degreeToVerticesFast.get(d);
                if (!vertices.isEmpty()) {
                    minDegreeVertex = vertices.intIterator().nextInt();
                    minDegree = d;
                    break;
                }
            }
            assert minDegree != -1;
            // update neighbors
            for (int v: h.neighborhoodCopy(minDegreeVertex)) {
                int degree = h.degree(v);
                assert degree >= 1;
                degreeToVerticesFast.get(degree).remove(v);
                degreeToVerticesFast.get(degree - 1).add(v);
            }
            // remove the vertex that was found
            degreeToVerticesFast.get(minDegree).remove(minDegreeVertex);
            h.removeVertex(minDegreeVertex);
            // add the vertex to the ordering
            ordering[i] = minDegreeVertex;
        }
        return ordering;
    }

    public static DynGraph[] applyTo(DynGraph g, int maxDegeneracy) {
        int[] ordering = degeneracyOrdering(g, maxDegeneracy);
        if (ordering == null) {
            // degeneracy is too high; do not apply the Turing kernel
            DynGraph[] singleton = new DynGraph[1];
            singleton[0] = g;
            return singleton;
        }
        int n = g.numVertices();
        assert ordering.length == n;
        DynGraph[] result = new DynGraph[n];
        // construct mapping from any vertex to its index in the ordering
        Int2IntOpenHashMap vertexToIndex = new Int2IntOpenHashMap(n);
        for (int i = 0; i < n; i++) {
            vertexToIndex.put(ordering[i], i);
        }
        assert vertexToIndex.keySet().equals(g.vertexSetCopy());
        // fill result with induced subgraphs
        for (int i = 0; i < n; i++) {
            int v = ordering[i];
            IntOpenHashSet vertices = g.neighborhoodCopy(v);
            int iCopy = i; // necessary for lambda expression
            vertices.removeIf(u -> vertexToIndex.get(u) < iCopy);
            vertices.add(v);
            result[i] = g.inducedSubgraph(vertices);
        }
        return result;
    }
}
