package org.alg_eng_2;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIntImmutablePair;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;

import java.util.Arrays;

public class MassiveGraph implements DynGraph {
    private IntArray[] neighbors;
    private int num_v;
    private int num_e;
    private IntOpenHashSet vertices;

    public MassiveGraph(int num_v, int num_e, IntArray edgesLeft, IntArray edgesRight) {
        this.num_v = num_v;
        this.num_e = num_e;
        this.neighbors = new IntArray[this.num_v+1];
        this.vertices = new IntOpenHashSet();
        for (int i = 0; i < this.num_v; i++) {
            this.vertices.add(i+1);
        }

        for (int e = 0; e < num_e; e++) {
            int u = edgesLeft.data[e];
            int v = edgesRight.data[e];

            assert u != v : "looping edge";
            assert !this.hasEdge(u,v) : "edge already present";

            this.neighbors[u].add(v);
            this.neighbors[v].add(u);
        }

        for (int i = 1; i <= this.num_v; i++) {
            Arrays.sort(this.neighbors[i].data);
        }
    }

    @Override
    public int numVertices() {
        return this.num_v;
    }

    @Override
    public int numEdges() {
        return this.num_e;
    }

    @Override
    public boolean hasVertex(int v) {
        return this.neighbors[v] != null;
    }

    @Override
    public int degree(int v) {
        return this.neighbors[v].len;
    }

    @Override
    public boolean hasEdge(int u, int v) {
        return Arrays.binarySearch(this.neighbors[u].data,v) >= 0;
    }

    @Override
    public IntOpenHashSet vertexSetCopy() {
        return new IntOpenHashSet(this.vertices);
    }

    @Override
    public IntOpenHashSet neighborhoodCopy(int v) {
        return new IntOpenHashSet(this.neighbors[v].data);
    }

    @Override
    public void removeVertex(int v) {

    }

    @Override
    public void removeEdge(int u, int v) {

    }

    @Override
    public DynGraph inducedSubgraph(IntOpenHashSet vertices) {
        return null;
    }

    @Override
    public void sortDegree() {

    }

    @Override
    public int to_original(int v) {
        return 0;
    }

    @Override
    public int from_original(int ov) {
        return 0;
    }

    @Override
    public int maxVertex() {
        return 0;
    }

    @Override
    public int maxDegree() {
        return 0;
    }

    @Override
    public int maxDegreeVertex() {
        return 0;
    }

    @Override
    public int minDegree() {
        return 0;
    }

    @Override
    public int minDegreeVertex() {
        return 0;
    }

    @Override
    public IntOpenHashSet fastClique(int vertex, int minSize) {
        return null;
    }

    @Override
    public IntOpenHashSet dynFastClique(int vertex, int minSize) {
        return null;
    }

    @Override
    public IntOpenHashSet reduceVertices012() {
        return null;
    }

    @Override
    public void reduceVerticesAndEdges(int lowerDegreeBound) {

    }

    @Override
    public void dynReduceVerticesAndEdges(int lowerDegreeBound) {

    }

    @Override
    public DynGraph copy() {
        return null;
    }

    @Override
    public boolean checkOriginalMappings() {
        return false;
    }

    @Override
    public boolean isDegreeSorted() {
        return false;
    }

    @Override
    public boolean isClique(IntOpenHashSet c) {
        return false;
    }
}
