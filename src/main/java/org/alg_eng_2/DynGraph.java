package org.alg_eng_2;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;

public interface DynGraph {
    int numVertices();
    int numEdges();
    boolean hasVertex(int v);
    int degree(int v);
    boolean hasEdge(int u, int v);
    IntOpenHashSet vertexSetCopy();
    IntOpenHashSet neighborhoodCopy(int v);
    void removeVertex(int v);
    void removeEdge(int u, int v);
    DynGraph inducedSubgraph(IntOpenHashSet vertices);
    void sortDegree();
    int to_original(int v);
    int from_original(int ov);
    int maxVertex();
    int maxDegree();
    int maxDegreeVertex();
    int minDegree();
    int minDegreeVertex();
    IntOpenHashSet fastClique(int vertex, int minSize);
    IntOpenHashSet dynFastClique(int vertex, int minSize);
    IntOpenHashSet reduceVertices012();
    void reduceVerticesAndEdges(int lowerDegreeBound);
    void dynReduceVerticesAndEdges(int lowerDegreeBound);
    DynGraph copy();
    boolean checkOriginalMappings();
    boolean isDegreeSorted();
    boolean isClique(IntOpenHashSet c);
}
