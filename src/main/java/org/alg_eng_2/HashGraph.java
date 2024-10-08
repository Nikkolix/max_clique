package org.alg_eng_2;

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;

import java.util.*;

import static org.alg_eng_2.PARAMETERS.MAX_V_FOR_EDGE_REDUCTION;

public class HashGraph implements DynGraph {
    private final Int2ObjectOpenHashMap<IntOpenHashSet> neighbors;
    private final Int2IntOpenHashMap to_original;
    private final Int2IntOpenHashMap from_original;
    private final IntOpenHashSet vertices;
    private int num_e;
    private int num_v;

    //PRIVATE

    private HashGraph(HashGraph g) {
        this.num_v = g.num_v;
        this.num_e = g.num_e;
        this.vertices = g.vertices.clone();
        this.to_original = g.to_original.clone();
        this.from_original = g.from_original.clone();
        this.neighbors = new Int2ObjectOpenHashMap<>();
        for (int v : g.vertices) {
            this.neighbors.put(v,g.neighbors.get(v).clone());
        }
    }

    private HashGraph() {
        this.neighbors = new Int2ObjectOpenHashMap<>();
        this.to_original = new Int2IntOpenHashMap();
        this.from_original = new Int2IntOpenHashMap();
        this.vertices = new IntOpenHashSet();
        this.num_e = 0;
        this.num_v = 0;
    }

    private void addVertex(int v) {
        this.addVertex(v,v);
    }

    private void addVertex(int v, int ov) {
        assert !this.hasVertex(v);

        this.neighbors.put(v,new IntOpenHashSet());
        this.vertices.add(v);
        this.to_original.put(v,ov);
        this.from_original.put(ov,v);

        this.num_v++;
    }

    private void addEdge(int u, int v) {
        assert this.hasVertex(u);
        assert this.hasVertex(v);
        assert u != v;
        assert !this.hasEdge(u,v);

        this.neighbors.get(u).add(v);
        this.neighbors.get(v).add(u);
        this.num_e++;
    }

    //PUBLIC

    public HashGraph(int num_v, int num_e, IntArray edgesLeft, IntArray edgesRight) {
        this.neighbors = new Int2ObjectOpenHashMap<>();
        this.to_original = new Int2IntOpenHashMap();
        this.from_original = new Int2IntOpenHashMap();
        this.vertices = new IntOpenHashSet();
        this.num_e = num_e;
        this.num_v = num_v;

        for (int v = 1 ; v <= num_v; v++) {
            this.neighbors.put(v,new IntOpenHashSet());
            this.to_original.put(v,v);
            this.from_original.put(v,v);
            this.vertices.add(v);
        }

        for (int e = 0; e < num_e; e++) {
            int u = edgesLeft.data[e];
            int v = edgesRight.data[e];

            assert this.hasVertex(u);
            assert this.hasVertex(v);
            assert u != v;

            this.neighbors.get(u).add(v);
            this.neighbors.get(v).add(u);
        }
    }

    public int numVertices() {
        return this.num_v;
    }

    public int numEdges() {
        return this.num_e;
    }

    public boolean hasVertex(int v) {
        return this.vertices.contains(v);
    }

    public void removeVertex(int v) {
        assert this.hasVertex(v);

        this.num_e -= this.degree(v);
        for (int u: this.neighbors.get(v)) {
            this.neighbors.get(u).remove(v);
        }
        this.neighbors.remove(v);
        int ov = this.to_original.get(v);
        this.to_original.remove(v);
        this.from_original.remove(ov);
        this.vertices.remove(v);
        this.num_v--;
    }

    public void removeEdge(int u, int v) {
        assert this.hasVertex(u);
        assert this.hasVertex(v);
        assert this.hasEdge(u,v);

        this.neighbors.get(u).remove(v);
        this.neighbors.get(v).remove(u);

        this.num_e--;
    }

    public int degree(int v) {
        assert this.hasVertex(v);

        return this.neighbors.get(v).size();
    }

    public boolean hasEdge(int u, int v) {
        assert this.hasVertex(u);
        assert this.hasVertex(v);

        return this.neighbors.get(u).contains(v);
    }

    public IntOpenHashSet vertexSetCopy() {
        return new IntOpenHashSet(this.vertices);
    }

    public IntOpenHashSet neighborhoodCopy(int v) {
        assert this.hasVertex(v) : v + " must be in the vertices " + this.vertexSetCopy() + " of the graph";
        assert this.neighbors != null;
        assert this.neighbors.get(v) != null;

        return new IntOpenHashSet(this.neighbors.get(v));
    }

    public HashGraph inducedSubgraph(IntOpenHashSet vertices) {
        assert this.vertexSetCopy().containsAll(vertices);

        HashGraph g = new HashGraph();
        for (int v: vertices) {
            int ov = this.to_original.get(v);
            g.addVertex(v, ov);
        }
        for (int u: vertices) {
            IntOpenHashSet intersection = this.neighborhoodCopy(u);
            intersection.retainAll(vertices);
            for (int v: intersection) {
                g.addEdge(u, v);
            }
        }
        return g;
    }


    public void sortDegree() {
        ObjectOpenHashSet<IntIntImmutablePair> DV = new ObjectOpenHashSet<>();

        for (int v : this.vertexSetCopy()) {
            DV.add(new IntIntImmutablePair(this.degree(v),v));
        }

        IntIntImmutablePair[] DVA = DV.toArray(new IntIntImmutablePair[0]);
        Arrays.sort(DVA,Comparator.comparing(IntIntImmutablePair::leftInt));

        int[] vertices = this.vertices.toArray(new int[0]);
        Arrays.sort(vertices);

        Int2IntOpenHashMap mapping = new Int2IntOpenHashMap();
        for (int i = 0; i < DVA.length; i++) {
            mapping.put(DVA[i].rightInt(),vertices[i]);
        }

        Int2IntOpenHashMap to_original = this.to_original.clone();

        Int2ObjectOpenHashMap<IntOpenHashSet> neighbors = this.neighbors.clone();
        for (int v : vertices) {
            IntOpenHashSet neighborhood = new IntOpenHashSet();
            for (int n : neighbors.get(v)) {
                neighborhood.add(mapping.get(n));
            }
            this.neighbors.put(mapping.get(v),neighborhood);


            int ov = to_original.get(v);
            this.to_original.put(mapping.get(v),ov);
            this.from_original.put(ov,mapping.get(v));
        }
    }

    public int to_original(int v) {
        assert this.to_original.containsKey(v);

        return this.to_original.get(v);
    }

    public int from_original(int ov) {
        assert this.from_original.containsKey(ov);

        return this.from_original.get(ov);
    }

    public int maxVertex() {
        int name = -1;
        for (int v : this.vertexSetCopy()) {
            name = Math.max(name,v);
        }
        return name;
    }

    public int maxDegree() {
        int maxDeg = -1;
        for (int v : this.vertexSetCopy()) {
            maxDeg = Math.max(maxDeg, this.degree(v));
        }
        return maxDeg;
    }

    public int maxDegreeVertex() {
        int maxDeg = -1;
        int maxDegVertex = -1;
        for (int v : this.vertexSetCopy()) {
            if (maxDeg < this.degree(v)) {
                maxDeg = this.degree(v);
                maxDegVertex = v;
            }
        }
        return maxDegVertex;
    }

    public int minDegree() {
        int minDeg = this.numVertices()+1;
        for (int v : this.vertexSetCopy()) {
            minDeg = Math.min(minDeg, this.degree(v));
        }
        return minDeg;
    }

    public int minDegreeVertex() {
        int minDeg = this.numVertices()+1;
        int minDegVertex = -1;
        for (int v : this.vertexSetCopy()) {
            if (minDeg > this.degree(v)) {
                minDeg = this.degree(v);
                minDegVertex = v;
            }
        }
        return minDegVertex;
    }

    public HashGraph copy() {
        return new HashGraph(this);
    }

    public IntOpenHashSet dynFastClique(int vertex, int minSize) {

        int size = this.maxVertex() + 1;

        BitSet cliqueCandidates = new BitSet(size);
        for (int n : this.neighbors.get(vertex)) {
            cliqueCandidates.set(n, true);
        }

        int[] degrees = new int[size];
        for (int v : this.vertices) {
            degrees[v] = this.neighbors.get(v).size();
        }

        int outSize = 1;
        BitSet clique = new BitSet(size);
        clique.set(vertex, true);

        int[] vertexNeighbors = new int[0];

        while (true) {
            int bestNeighbor = -1;
            int bestNeighborSize = -1;

            if (outSize + this.degree(vertex) < minSize) {
                break;
            }

            // loop over vertices starting from max deg (if sorted before)
            vertexNeighbors = this.neighbors.get(vertex).toArray(vertexNeighbors);
            Arrays.sort(vertexNeighbors);
            for (int i = degrees[vertex] - 1; i >= 0; i--) {
                int neighbor = vertexNeighbors[i];

                // neighbor is already inside out clique
                if (clique.get(neighbor)) {
                    continue;
                }

                // neighbor is no valid candidate
                if (!cliqueCandidates.get(neighbor)) {
                    continue;
                }

                // neighbour has not enough neighbors for a clique bigger than min size
                if (degrees[neighbor] + outSize < minSize) {
                    continue;
                }

                // neighbor has not enough neighbors to have more common neighbors than the best before
                if (degrees[neighbor] <= bestNeighborSize) {
                    continue;
                }

                BitSet tmp = (BitSet) cliqueCandidates.clone();

                BitSet neighbors = new BitSet(size);
                for (int n : this.neighbors.get(neighbor)) {
                    neighbors.set(n, true);
                }

                tmp.and(neighbors);
                int in_candidates = tmp.cardinality();

                // if more common neighbors, than update best neighbor and its size
                if (bestNeighborSize < in_candidates) {
                    bestNeighborSize = in_candidates;
                    bestNeighbor = neighbor;
                }
            }

            if (bestNeighbor == -1) {
                break;
            }

            BitSet neighbors = new BitSet(size);
            for (int n : this.neighbors.get(bestNeighbor)) {
                neighbors.set(n, true);
            }

            cliqueCandidates.set(bestNeighbor, false);
            cliqueCandidates.and(neighbors);

            clique.set(bestNeighbor, true);
            outSize++;
            vertex = bestNeighbor;
        }

        IntOpenHashSet out = new IntOpenHashSet();
        for (int i = 0; i < clique.length(); i++) {
            if (clique.get(i)) {
                out.add(i);
            }
        }

        return out;
    }

//    public HashSet<Integer> tabuClique(int maxIter) {
//        int size = this.maxVertexName();
//
//        BitSet cliqueCandidates = new BitSet(size);
//        for (int v : this.vertexSet()) {
//            cliqueCandidates.set(v, true);
//        }
//
//        BitSet clique = new BitSet(size);
//        for (int v : this.vertexSet()) {
//            clique.set(v,false);
//        }
//        int candidateSize = size;
//        int maxCliquesize = 0;
//
//        while (candidateSize > 0) {
//            int vertex = -1;
//            Random rand = new Random();
//            int select = rand.nextInt(candidateSize);
//            int i = 0;
//            while (i < size) {
//                if (!cliqueCandidates.get(i)) {
//                i += 1;}
//                if (select > 0 && cliqueCandidates.get(i)) {
//                    i += 1;
//                    select --;
//                }
//                if (select == 0 && cliqueCandidates.get(i)) {
//                    vertex = i;
//                    i = size;
//                }
//            }
//            clique.set(i, true);
//            maxCliquesize += 1;
//            int j = 0;
//            while (j < size) {
//                if (!hasEdge(j, i) && cliqueCandidates.get(j)){
//                    cliqueCandidates.set(j, false);
//                    candidateSize --;
//                }
//                j ++;
//            }
//        }
//
//        BitSet largestClique = clique;
//        int cliqueSize = maxCliquesize;
//        int[] tabuList = new int[size];
//        int i = 0;
//        while (i < size) {
//            tabuList[i] = 0;
//            i += 1;
//        }
//        int iter = 0;
//        int[] NSzero = new int[size];
//        int NSzerolen = 0;
//        int[] NSone = new int[size];
//        int NSonelen = 0;
//        int[] NStwo = new int[size];
//        int NStwolen = 0;
//        int[] NSmore = new int[size];
//        int NSmorelen = 0;
//        int[] mapDegree = new int[size];
//        i = 0;
//        int[] expandDegree = new int[size];
//        int[] diverseDegree = new int[size];
//        while (i < size){
//            if (!clique.get(i)){
//                int deg = degree(i);
//                for (int n : this.neighbors.get(i)){
//                    if (!clique.get(n)){
//                        deg --;
//                    }
//                }
//                mapDegree[i] = maxCliquesize - deg;
//                diverseDegree[i] = size - mapDegree[i] - degree(i) - 1;
//                if (mapDegree[i] > 2) {
//                    NSmore[NSmorelen] = i;
//                    NSmorelen += 1;
//                }
//                else if (mapDegree[i] == 2) {
//                    NStwo[NStwolen] = i;
//                    NStwolen += 1;
//                }
//                else if (mapDegree[i] == 1) {
//                    NSone[NSonelen] = i;
//                    NSonelen += 1;
//                }
//            }
//            i ++;
//        }
//        i = 0;
//        while (i < size) {
//            if (clique.get(i)){
//                int j = 0;
//                while (j < size) {
//                    if (!clique.get(j) && !hasEdge(i, j) && mapDegree[j] == 1){
//                        expandDegree[i] += 1;
//                    }
//                    j +=1;
//                }
//            }
//        }
//
//        while (iter < maxIter) {
//            if (NSzerolen > 0) {
//                int vertex = NSzero[NSzerolen -1];
//                clique.set(vertex, true);
//                NSzerolen -= 1;
//                cliqueSize += 1;
//                if (cliqueSize > maxCliquesize) {
//                    largestClique = clique;
//                    maxCliquesize = cliqueSize;
//                }
//                mapDegree[j_1] = 0;
//                expandDegree[vertex] = 0;
//                i = 0;
//                while (i < size) {
//                    if ( !clique.get(i) && !hasEdge(i, vertex)) {
//                        diverseDegree[i] -= 1;
//                        mapDegree[i] += 1;
//                        if (mapDegree[i] == 1) {
//                            expandDegree[vertex] += 1;
//                            NSone[NSonelen] = i;
//                            NSonelen ++;
//                        }
//                        if (mapDegree[i] == 2) {
//                            expandDegree[vertex] -= 1;
//                            NStwo[NStwolen] = i;
//                            NStwolen ++;
//                            int j = 0;
//                            while (j < NSonelen) {
//                                if (NSone[j] == i) {
//                                    NSone[j] = NSone[NSonelen - 1];
//                                    NSonelen --;
//                                    j = NSonelen;
//                                }
//                                j ++;
//                            }
//                        }
//                        if (mapDegree[i] == 3) {
//                            NSmore[NSmorelen] = i;
//                            NSmorelen ++;
//                            int j = 0;
//                            while (j < NStwolen) {
//                                if (NStwo[j] == i) {
//                                    NStwo[j] = NStwo[NStwolen - 1];
//                                    NStwolen --;
//                                    j = NStwolen;
//                                }
//                                j ++;
//                            }
//                        }
//                    }
//                    i+=1;
//                }
//            }
//
//            else if (NSonelen != 0) {
//                BitSet notdisqualified = null;
//                if (NSonelen > NStwolen + NSmorelen) {
//                    notdisqualified = new BitSet(NSonelen);
//                    i = 0;
//                    while (i < NSonelen) {
//                        int vert = NSone[i];
//                        int j = 0;
//                        while (!clique.get(j) || hasEdge(vert, j)) {
//                            j+= 1;
//                        }
//                        int k = 0;
//                        while (k < size) {
//                            if (!hasEdge(k, j) && k != vert && k != j){
//                                notdisqualified.set(i, true);
//                            }
//                        }
//                        i += 1;
//                    }
//                }
//                int vertex = -1;
//                int maxexp = -1;
//                int maxdiv = -1;
//                i = 0;
//                int i_1 = -1;
//                int j_1 = -1;
//                while (i < NSonelen) {
//                    if (notdisqualified.get(NSone[i]) && tabuList[NSone[i]] <= iter) {
//                        int j = 0;
//                        while (!clique.get(j) && hasEdge(j, NSone[i])) {
//                            j ++;
//                        }
//                        if (expandDegree[j] > maxexp || (expandDegree[j] >= maxexp && diverseDegree[NSone[i]] > maxdiv)) {
//                            vertex = NSone[i];
//                            i_1 = i;
//                            j_1 = j;
//                            maxexp = expandDegree[j];
//                            maxdiv = diverseDegree[NSone[i]];
//                        }
//                    }
//                    i ++;
//                }
//                if (vertex != -1) {
//                    clique.set(j_1, false);
//                    clique.set(vertex, true);
//                    mapDegree[j_1] = 0;
//                    expandDegree[vertex] = 0;
//                    diverseDegree[j_1] = 0;
//                    NSzero[NSzerolen] = j_1;
//                    NSzerolen += 1;
//                    NSone[i_1] = NSone[NSonelen - 1];
//                    NSonelen -= 1;
//                    if (NSonelen < NStwolen + NSmorelen) {
//                        Random rand = new Random(NSonelen);
//                        tabuList[j_1] = iter + 10 + rand();
//                    }
//                    else {
//                        tabuList[j_1] = iter + NSonelen;
//                    }
//                    i = 0;
//                    while (i < size) {
//                        if ( !clique.get(i) && !hasEdge(i, vertex)) {
//                            diverseDegree[i] -= 1;
//                            mapDegree[i] += 1;
//                            if (mapDegree[i] == 1) {
//                                expandDegree[vertex] += 1;
//                                NSone[NSonelen] = i;
//                                NSonelen ++;
//
//                            }
//                            if (mapDegree[i] == 2) {
//                                expandDegree[vertex] -= 1;
//                                NStwo[NStwolen] = i;
//                                NStwolen ++;
//                                int j = 0;
//                                while (j < NSonelen) {
//                                    if (NSone[j] == i) {
//                                        NSone[j] = NSone[NSonelen - 1];
//                                        NSonelen --;
//                                        j = NSonelen;
//                                    }
//                                    j ++;
//                                }
//                            }
//                            if (mapDegree[i] == 3) {
//                                NSmore[NSmorelen] = i;
//                                NSmorelen ++;
//                                int j = 0;
//                                while (j < NStwolen) {
//                                    if (NStwo[j] == i) {
//                                        NStwo[j] = NStwo[NStwolen - 1];
//                                        NStwolen --;
//                                        j = NStwolen;
//                                    }
//                                    j ++;
//                                }
//                            }
//                        }
//                        if (i != j_1 && !clique.get(i) && !hasEdge(i, j_1)) {
//                            diverseDegree[i] += 1;
//                            mapDegree[i] --;
//                            diverseDegree[j_1] ++;
//                            if (mapDegree[i] == 0) {
//                                NSzero[NSzerolen] = i;
//                                NSzerolen ++;
//                                int j = 0;
//                                while (j < NSonelen) {
//                                    if (NSone[j] == i) {
//                                        NSone[j] = NSone[NSonelen - 1];
//                                        NSonelen --;
//                                        j = NSonelen;
//                                    }
//                                    j++;
//                                }
//                            }
//                            if (mapDegree[i] == 1) {
//                                NSone[NSonelen] = i;
//                                NSonelen ++;
//                                int j = 0;
//                                while (!clique.get(j) || hasEdge(j, i)) {
//                                    j+=1;
//                                }
//                                expandDegree[j] += 1;
//                                j = 0;
//                                while (j < NStwolen) {
//                                    if (NStwo[j] == i) {
//                                        NStwo[j] = NStwo[NStwolen - 1];
//                                        NStwolen --;
//                                        j = NStwolen;
//                                    }
//                                    j++;
//                                }
//                            }
//                            if (mapDegree[i] == 2) {
//                                NStwo[NStwolen] = i;
//                                NStwolen ++;
//                                int j = 0;
//                                while (j < NSmorelen) {
//                                    if (NSmore[j] == i) {
//                                        NSmore[j] = NSmore[NSmorelen - 1];
//                                        NSmorelen --;
//                                        j = NSmorelen;
//                                    }
//                                    j++;
//                                }
//                            }
//                        }
//                        i ++;
//                    }
//                }
//            }
//            if (vertex = -1 || (NSonelen = 0 && NSzerolen = 0)) {
//                if (NSonelen > NStwolen + NSmorelen && NSmorelen != 0) {
//                    int moreswap = 1;
//                }
//                else {
//                    if (NStwolen != 0) {
//                        int moreswap = 0;
//                    }
//                    if (NStwolen != 0 && NSmorelen != 0) {
//                        Random rand = new Random(1);
//                        moreswap = rand();
//                    }
//                }
//                vertex = -1;
//                locmax = -1;
//                i = 0;
//                if (moreswap == 0) {
//                    while (i < NStwolen) {
//                        if (tabuList[NStwo[i]] <= iter && diverseDegree[NStwo[i]] > locmax) {
//                            locmax = diverseDegree[NStwo[i]];
//                            vertex = NStwo[i];
//                            int i_2 = i
//                        }
//                        i ++;
//                    }
//                    if (vertex != -1) {
//                        NStwo[i_2] = NStwo[NStwolen - 1];
//                        NStwolen -= 1;
//                    }
//                }
//                if (moreswap == 1) {
//                    while (i < NSmorelen) {
//                        if (tabuList[NSmore[i]] <= iter && diverseDegree[NSmore[i]] > locmax) {
//                            locmax = diverseDegree[NSmore[i]];
//                            vertex = NSmore[i];
//                            int i_2 = i
//                        }
//                        i ++;
//                    }
//                    if (vertex != -1) {
//                        NSmore[i_2] = NSmore[NSmorelen - 1];
//                        NSmorelen -= 1;
//                    }
//                }
//                if (vertex != -1) {
//                    clique.set(vertex, true)
//                    cliqueSize ++;
//                    int[] swapSet = new int[size];
//                    int swapSize = 0;
//                    i = 0;
//                    while (i < size) {
//                        if (clique.get(i) && !hasEdge(i, vertex)) {
//                            clique.set(i, false);
//                            cliqueSize --;
//                            swapSet[swapSize] = i;
//                            swapSize ++;
//                            if (NSonelen < NStwolen + NSmorelen) {
//                                Random rand = new Random(NSonelen);
//                                tabuList[i] = iter + 10 + rand();
//                            }
//                            else {
//                                tabuList[i] = iter + NSonelen;
//                            }
//                        }
//                        i ++;
//                    }
//                }
//                if (vertex != -1) {
//                    i = 0;
//                    while (i < swapsize) {
//                        mapDegree[swapSet[i]] = 0;
//                        diverseDegree[swapSet[i]] = 0;
//                        NSzero[NSzerolen] = swapSet[i];
//                        NSzerolen += 1;
//                        i += 1;
//                    }
//                    expandDegree[vertex] = 0;
//                    i = 0;
//                    while (i < size) {
//                        if ( !clique.get(i) && !hasEdge(i, vertex)) {
//                            diverseDegree[i] -= 1;
//                            mapDegree[i] += 1;
//                            if (mapDegree[i] == 1) {
//                                expandDegree[vertex] += 1;
//                                NSone[NSonelen] = i;
//                                NSonelen ++;
//
//                            }
//                            if (mapDegree[i] == 2) {
//                                expandDegree[vertex] -= 1;
//                                NStwo[NStwolen] = i;
//                                NStwolen ++;
//                                j = 0;
//                                while (j < NSonelen) {
//                                    if (NSone[j] == i) {
//                                        NSone[j] = NSone[NSonelen - 1];
//                                        NSonelen --;
//                                        j = NSonelen;
//                                    }
//                                    j ++;
//                                }
//                            }
//                            if (mapDegree[i] == 3) {
//                                NSmore[NSmorelen] = i;
//                                NSmorelen ++;
//                                j = 0;
//                                while (j < NSmtwolen) {
//                                    if (NStwo[j] == i) {
//                                        NStwo[j] = NStwo[NStwolen - 1];
//                                        NStwolen --;
//                                        j = NStwolen;
//                                    }
//                                    j ++;
//                                }
//                            }
//                        }
//                        j = 0;
//                        while (j < swapsize) {
//                            if (i != swapSet[j] && !clique.get(i) && !hasEdge(i, swapSet[j])) {
//                                diverseDegree[i] += 1;
//                                mapDegree[i] --;
//                                diverseDegree[swapSet[j]] ++;
//                                if (mapDegree[i] == 0) {
//                                    NSzero[NSzerolen] = i;
//                                    NSzerolen ++;
//                                    int k = 0;
//                                    while (k < NSonelen) {
//                                        if (NSone[k] == i) {
//                                            NSone[k] = NSone[NSonelen - 1];
//                                            NSonelen --;
//                                            k = NSonelen;
//                                        }
//                                        k++;
//                                    }
//                                }
//                            }
//                            j +=1
//                        }
//                        if (i != swapSet[j] && !clique.get(i) && !hasEdge(i, swapSet[j])) {
//                            if (mapDegree[i] == 1) {
//                                NSone[NSonelen] = i;
//                                NSonelen ++;
//                                k = 0;
//                                while (!clique.get(k) || hasEdge(k, i)) {
//                                    k+=1
//                                }
//                                expandDegree[k] += 1;
//                                k = 0;
//                                while (k < NStwolen) {
//                                    if (NStwo[k] == i) {
//                                        NStwo[k] = NStwo[NStwolen - 1];
//                                        NStwolen --;
//                                        k = NStwolen;
//                                    }
//                                    k++;
//                                }
//                            }
//                            if (mapDegree[i] == 2) {
//                                NStwo[NStwolen] = i;
//                                NStwolen ++;
//                                k = 0;
//                                while (k < NSmorelen) {
//                                    if (NSmore[k] == i) {
//                                        NSmore[k] = NSmore[NSmorelen - 1];
//                                        NSmorelen --;
//                                        k = NSmorelen;
//                                    }
//                                    j++;
//                                }
//                            }
//                        }
//                        i += 1;
//                    }
//                }
//
//
//            }
//            iter ++;
//        }
//        return(maxCliquesize);
//    }

    public IntOpenHashSet fastClique(int vertex, int minSize) {

        int size = this.maxVertex()+1;

        BitSet cliqueCandidates = new BitSet(size);
        for (int n : this.neighbors.get(vertex)) {
            cliqueCandidates.set(n,true);
        }

        int[] degrees = new int[size];
        BitSet[] neighbors = new BitSet[size];
        for (int v : this.vertexSetCopy()) {
            neighbors[v] = new BitSet(size);
            for (int n : this.neighbors.get(v)) {
                neighbors[v].set(n,true);
            }
            degrees[v] = this.neighbors.get(v).size();
        }

        int outSize = 1;
        BitSet clique = new BitSet(size);
        clique.set(vertex,true);

        while (true) {
            int bestNeighbor = -1;
            int bestNeighborSize = -1;

            if (outSize + this.degree(vertex) < minSize) {
                break;
            }

            // loop over vertices starting from max deg (if sorted before)
            ArrayList<Integer> vertexNeighbors = new ArrayList<>(this.neighbors.get(vertex));
            for (int i = vertexNeighbors.size()-1; i >= 0 ; i--) {
                int neighbor = vertexNeighbors.get(i);

                // neighbor is already inside out clique
                if (clique.get(neighbor)) {
                    continue;
                }

                // neighbor is no valid candidate
                if (!cliqueCandidates.get(neighbor)) {
                    continue;
                }

                // neighbour has not enough neighbors for a clique bigger than min size
                if (degrees[neighbor] + outSize < minSize) {
                    continue;
                }

                // neighbor has not enough neighbors to have more common neighbors than the best before
                if (degrees[neighbor] <= bestNeighborSize) {
                    continue;
                }

                BitSet tmp = (BitSet) cliqueCandidates.clone();
                tmp.and(neighbors[neighbor]);
                int in_candidates = tmp.cardinality();

                // if more common neighbors, than update best neighbor and its size
                if (bestNeighborSize < in_candidates) {
                    bestNeighborSize = in_candidates;
                    bestNeighbor = neighbor;
                }
            }

            if (bestNeighbor == -1) {
                break;
            }

            cliqueCandidates.set(bestNeighbor,false);
            cliqueCandidates.and(neighbors[bestNeighbor]);

            clique.set(bestNeighbor,true);
            outSize++;
            vertex = bestNeighbor;
        }

        IntOpenHashSet out = new IntOpenHashSet();
        for (int i = 0; i < clique.length(); i++) {
            if (clique.get(i)) {
                out.add(i);
            }
        }

        return out;
    }

    public IntOpenHashSet reduceVertices0() {
        IntOpenHashSet out = new IntOpenHashSet();
        IntOpenHashSet vertices = this.vertexSetCopy();
        for (int v : vertices) {
            if (this.degree(v) == 0) {
                this.removeVertex(v);
                if (out.isEmpty()) {
                    out.add(v);
                }
            }
        }
        return out;
    }

    public IntOpenHashSet reduceVertices1() {
        IntOpenHashSet out = new IntOpenHashSet();
        IntOpenHashSet vertices = this.vertexSetCopy();
        for (int v : vertices) {
            if (this.hasVertex(v)) {
                if (this.degree(v) == 1) {
                    if (out.isEmpty()) {
                        out.add(v);
                        out.add((int) this.neighbors.get(v).toArray()[0]);
                    }
                    this.removeVertex(v);
                }
            }
        }
        return out;
    }

    public IntOpenHashSet reduceVertices2() {
        IntOpenHashSet out = new IntOpenHashSet();
        IntOpenHashSet vertices = this.vertexSetCopy();
        for (int v : vertices) {
            if (this.hasVertex(v)) {
                if (this.degree(v) == 2) {
                    int n1 = (int) this.neighbors.get(v).toArray()[0];
                    int n2 = (int) this.neighbors.get(v).toArray()[1];
                    this.removeVertex(v);
                    if (out.isEmpty()) {
                        out.add(v);
                        out.add(n1);
                        if (this.hasEdge(n1,n2)) {
                            out.add(n2);
                        }
                    } else if (out.size() < 3 && this.hasEdge(n1,n2)) {
                        out.clear();
                        out.add(v);
                        out.add(n1);
                        out.add(n2);
                    }

                }
            }
        }
        return out;
    }

    public IntOpenHashSet reduceVertices012() {
        int size = -1;

        IntOpenHashSet out = new IntOpenHashSet();
        IntOpenHashSet clique = this.reduceVertices2();
        if (!clique.isEmpty()) {
            out = clique;
        }
        clique = this.reduceVertices1();
        if (clique.size() > out.size()) {
            out = clique;
        }

        while (size != this.numVertices()) {
            size = this.numVertices();
            clique = this.reduceVertices2();
            if (clique.size() > out.size()) {
                out = clique;
            }
            clique = this.reduceVertices1();
            if (clique.size() > out.size()) {
                out = clique;
            }
        }

        clique = this.reduceVertices0();
        if (clique.size() > out.size()) {
            out = clique;
        }
        return out;
    }

    public void dynReduceVerticesAndEdges(int lowerDegBound) {
        int was_vertices = Integer.MAX_VALUE;
        int was_edges = Integer.MAX_VALUE;

        int[] vertices = this.vertices.toArray(new int[0]);
        Arrays.sort(vertices);
        while (this.vertices.size() < was_vertices || was_edges < this.numEdges()) {
            was_vertices = this.numVertices();
            was_edges = this.numEdges();

            for (int v : vertices) {
                if (!this.hasVertex(v)) {
                    continue;
                }
                if (this.neighbors.get(v).size() < lowerDegBound) {
                    this.removeVertex(v);
                }
            }

            if (this.vertices.size() == was_vertices) {
                for (int i = 0; i < vertices.length && i < MAX_V_FOR_EDGE_REDUCTION; i++) {
                    int v = vertices[i];
                    if (this.hasVertex(v)) {
                        IntOpenHashSet neighborhood = this.neighbors.get(v).clone();
                        for (int n : neighborhood) {
                            if (this.hasVertex(n)) {
                                int cSize = 0;

                                IntOpenHashSet neighs = this.neighbors.get(n);
                                for (int u : neighs) {
                                    if (this.neighbors.get(v).contains(u)) {
                                        cSize++;
                                    }
                                }

                                if (cSize < lowerDegBound - 1) {
                                    this.neighbors.get(v).remove(n);
                                    this.neighbors.get(n).remove(v);
                                    if (this.neighbors.get(v).size() < lowerDegBound) {
                                        this.removeVertex(v);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public void reduceVerticesAndEdges(int lowerDegreeBound) {
        int size = this.maxVertex()+1;
        BitSet[] neighbors = new BitSet[size];
        for (int v : this.vertexSetCopy()) {
            neighbors[v] = new BitSet(size);
            for (int n : this.neighbors.get(v)) {
                neighbors[v].set(n,true);
            }
        }

        IntOpenHashSet vertices = this.vertexSetCopy();
        int was_vertices = vertices.size();
        int was_edges = this.numEdges();
        for (int v : vertices) {
            if (this.degree(v) < lowerDegreeBound) {
                for (int n : this.neighbors.get(v)) {
                    neighbors[n].set(v,false);
                }
                neighbors[v] = null;
                this.removeVertex(v);
            } else {
                IntOpenHashSet neighborhood = this.neighborhoodCopy(v);
                for (int n : neighborhood) {
                    if (this.hasVertex(n)) {
                        BitSet current = (BitSet) neighbors[v].clone();
                        current.and(neighbors[n]);
                        if (current.cardinality() < lowerDegreeBound-1) {
                            neighbors[v].set(n,false);
                            neighbors[n].set(v,false);
                            this.removeEdge(n,v);
                        }
                    }
                }
            }
        }

        if (this.vertexSetCopy().size() < was_vertices || was_edges < this.numEdges()) {
            reduceVerticesAndEdges(lowerDegreeBound);
        }
    }

    // OTHER

    public boolean isClique(IntOpenHashSet clique) {
        for (int v : clique) {
            for (int u : clique) {
                if (u == v) {
                    continue;
                }
                if (!this.hasEdge(u,v)) {
                    return false;
                }
            }
        }
        return true;
    }

    public boolean isDegreeSorted() {
        int[] vertices = this.vertexSetCopy().toArray(new int[0]);
        Arrays.sort(vertices);
        for (int i = 0; i < vertices.length-1; i++) {
            if (this.degree(vertices[i]) > this.degree(vertices[i+1])) {
                return false;
            }
        }
        return true;
    }

    public boolean checkOriginalMappings() {
        HashSet<Integer> seen = new HashSet<>();
        if (this.to_original.size() != this.from_original.size()) {
            return false;
        }
        for (int v : this.vertexSetCopy()) {
            if (seen.contains(this.to_original.get(v))) {
                return false;
            }
            int ov = this.to_original.get(v);
            seen.add(ov);
            if (this.from_original.get(ov) != v) {
                return false;
            }
        }

        return true;
    }

    public int[] degrees() {
        int[] degrees = new int[this.num_v];
        int[] vertices = this.vertices.toArray(new int[0]);
        Arrays.sort(vertices);
        for (int i = 0; i < vertices.length; i++) {
            degrees[i] = this.degree(vertices[i]);
        }
        return degrees;
    }

    public ArrayList<HashSet<Integer>> components() {
        if (this.numVertices() == 0) {
            return new ArrayList<>();
        }

        int start = this.maxDegreeVertex();
        int[] stack = new int[this.maxVertex()+1];
        int stack_len = 1;
        stack[0] = start;

        boolean[] seen = new boolean[this.maxVertex()+1];
        seen[start] = true;
        int seen_num = 1;

        ArrayList<HashSet<Integer>> out = new ArrayList<>();
        out.add(new HashSet<>());
        out.get(0).add(start);

        while (seen_num < this.numVertices()) {
            if (stack_len > 0) {
                stack_len--;
                int vertex = stack[stack_len];

                for (int n : this.neighbors.get(vertex)) {
                    if (!seen[n]) {
                        stack[stack_len] = n;
                        stack_len++;
                        out.get(out.size() - 1).add(n);
                        seen[n] = true;
                        seen_num++;
                    }
                }
            } else {
                int new_start = -1;
                for (int v : this.vertexSetCopy()) {
                    if (!seen[v]) {
                        new_start = v;
                        break;
                    }
                }
                out.add(new HashSet<>());
                out.get(out.size() - 1).add(new_start);
                seen[new_start] = true;
                seen_num++;
                stack[stack_len] = new_start;
                stack_len++;
            }
        }

        return out;
    }

    // MIN CUT

    private boolean bfs(
            int s,
            int t,
            Int2IntOpenHashMap parent,
            Int2ObjectOpenHashMap<IntOpenHashSet> used
    ) {
        BitSet visited = new BitSet(this.num_v);

        IntArrayFIFOQueue q = new IntArrayFIFOQueue();
        q.enqueue(s);
        visited.set(s);
        parent.put(s,-1);

        while (!q.isEmpty()) {
            int v = q.dequeueInt();
            for (int n : this.neighbors.get(v)) {
                if (!used.get(v).contains(n) && !visited.get(n)) {
                    q.enqueue(n);
                    visited.set(n);
                    parent.put(n,v);
                }
            }
        }

        return visited.get(t);
    }

    private void dfs(int s, BitSet visited, Int2ObjectOpenHashMap<IntOpenHashSet> used) {
        IntArrayList stack = new IntArrayList(this.num_v);
        stack.add(s);

        while (!stack.isEmpty()) {
            int v = stack.popInt();
            visited.set(v);
            for (int n : this.neighbors.get(v)) {
                if (!used.get(n).contains(v) && !visited.get(n)) {
                    stack.add(n);
                }
            }
        }
    }

    public ObjectOpenHashSet<IntIntImmutablePair> minCut(int s, int t) {
        Int2IntOpenHashMap parent = new Int2IntOpenHashMap();

        Int2ObjectOpenHashMap<IntOpenHashSet> used = new Int2ObjectOpenHashMap<>();
        for (int vertex : this.vertices) {
            used.put(vertex, new IntOpenHashSet());
        }

        while (bfs(s, t, parent, used)) {
            for (int v = t; v != s; v = parent.get(v)) {
                int u = parent.get(v);
                used.get(v).add(u);
                used.get(u).add(v);
            }
        }

        BitSet isVisited = new BitSet(this.num_v);
        dfs(s, isVisited, used);

        ObjectOpenHashSet<IntIntImmutablePair> out = new ObjectOpenHashSet<>();

        for (int v : this.vertices) {
            for (int n : this.neighbors.get(v)) {
                if (isVisited.get(v) && !isVisited.get(n)) {
                    out.add(new IntIntImmutablePair(v,n));
                }
            }
        }

        return out;
    }


    // MAX SHORTEST PATH

    private void dfs(int s, BitSet visited, IntOpenHashSet component) {
        IntArrayList stack = new IntArrayList(this.num_v);
        stack.add(s);
        visited.set(s);

        while (!stack.isEmpty()) {
            int v = stack.popInt();
            component.add(v);
            for (int n : this.neighbors.get(v)) {
                if (!visited.get(n)) {
                    stack.add(n);
                    visited.set(n);
                }
            }
        }
    }

    private void floydWarshall(int[][] distances)
    {
        int[] vertices = this.vertices.toArray(new int[0]);
        for (int u : vertices) {
            for (int v : vertices) {
                if (distances[v][u] == Integer.MAX_VALUE) {
                    continue;
                }
                for (int w : vertices) {
                    if (w > v || distances[u][w] == Integer.MAX_VALUE) {
                        continue;
                    }
                    distances[v][w] = Math.min(distances[v][w], distances[v][u] + distances[u][w]);
                    distances[w][v] = distances[v][w];
                }
            }
        }
    }

    private IntIntImmutablePair maxInThisComponent(IntOpenHashSet component, int[][] distance)
    {
        int maxDistance = Integer.MIN_VALUE;
        int u = -1;
        int v = -1;

        for (int cu : component) {
            for (int cv : component) {
                if (cu <= cv) {
                    continue;
                }
                int d = distance[cu][cv];
                if (d > maxDistance) {
                    maxDistance = d;
                    u = cu;
                    v = cv;
                }
            }
        }

        return new IntIntImmutablePair(u,v);
    }

    public ObjectArrayList<IntIntImmutablePair> maxShortestDistanceVertices()
    {
        BitSet visited = new BitSet(this.num_v);
        ObjectArrayList<IntOpenHashSet> components =  new ObjectArrayList<>();

        for (int v : this.vertices) {
            if (!visited.get(v)) {
                IntOpenHashSet component = new IntOpenHashSet();
                this.dfs(v,visited,component);
                components.add(component);
            }
        }

        int[][] distances = new int[this.maxVertex()+1][this.maxVertex()+1];
        for (int[] distance : distances) {
            Arrays.fill(distance, Integer.MAX_VALUE);
        }

        for (int v : this.vertices) {
            for (int n : this.neighbors.get(v)) {
                distances[v][n] = 1;
                distances[n][v] = 1;
            }
        }

        this.floydWarshall(distances);

        ObjectArrayList<IntIntImmutablePair> result = new ObjectArrayList<>(components.size());

        for (IntOpenHashSet component : components) {
            result.add(this.maxInThisComponent(component,distances));
        }
        return result;
    }
}
