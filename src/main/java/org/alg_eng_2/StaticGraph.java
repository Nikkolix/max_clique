package org.alg_eng_2;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

import static org.alg_eng_2.PARAMETERS.*;

public class StaticGraph {
    /**
     * total number of vertices in the graph
     */
    private final int num_v;

    /**
     * total number of edges in the graph
     */
    private final int num_e;

    /**
     * list of neighbors for every vertex <br>
     * every list of neighbors is sorted in ascending order,
     * contains no duplicates and only vertices of this graph [0,num_v-1] <br>
     * length is always num_v x num_v
     */
    private final IntArray[] neighbors;

    /**
     * adjacency matrix for this graph <br>
     * adj[i][j] 1, if edge between vertex i and j exist, 0 else
     * length is always num_v x num_v
     */
    private final byte[][] adj;
    private final boolean hasAdj;

    /**
     * the original vertex in the HashGraph used for construction for every vertex <br>
     * length is always num_v
     */
    private final int[] to_original;

    /**
     * contains the new vertex name for the original vertex name for every vertex <br>
     * length is always num_v
     */
    private final Int2IntOpenHashMap from_original;

    /**
     * the minimum degree of the graphs vertices
     */
    private final int minDeg;

    /**
     * the maximum degree of the graphs vertices
     */
    private final int maxDeg;

    public StaticGraph(DynGraph g) {
        this.num_v = g.numVertices();
        this.num_e = g.numEdges();
        this.neighbors = new IntArray[this.num_v];
        this.to_original = new int[this.num_v];
        this.from_original = new Int2IntOpenHashMap();

        int[] vertices = g.vertexSetCopy().toArray(new int[0]);
        Arrays.sort(vertices);
        int i = 0;
        for (int v : vertices) {
            this.to_original[i] = v;
            this.from_original.put(v,i);
            i++;
        }

        this.hasAdj = this.numVertices() <= MAX_V_ADJ_STATIC;

        if (this.hasAdj) {
            this.adj = new byte[this.num_v][this.num_v];
        } else {
            this.adj = null;
        }

        int minDeg = this.num_v;
        int maxDeg = 0;
        for (int v = 0; v < this.num_v; v++) {
            int[] neighbors = g.neighborhoodCopy(this.to_original[v]).toArray(new int[0]); //TODO two copies
            for (int j = 0; j < neighbors.length; j++) {
                neighbors[j] = this.from_original.get(neighbors[j]);
            }
            Arrays.sort(neighbors);
            this.neighbors[v] = new IntArray(neighbors);

            if (this.hasAdj) {
                for (int n : neighbors) {
                    this.adj[v][n] = 1;
                }
            }

            minDeg = Math.min(minDeg, this.neighbors[v].len);
            maxDeg = Math.max(maxDeg, this.neighbors[v].len);
        }
        this.minDeg = minDeg;
        this.maxDeg = maxDeg;
    }

    public int numVertices() {
        return this.num_v;
    }

    public int numEdges() {
        return this.num_e;
    }

    public boolean hasVertex(int v) {
        return v >= 0 && v < this.num_v;
    }

    /**
     * hasEdge
     * @param u vertex
     * @param v vertex
     * @return true if this graph has an edge that connects the vertices u and v
     */
    public boolean hasEdge(int u, int v) {
        return this.adj[u][v] == 1;
    }

    public int[] neighbors(int v) {
        assert v >= 0;
        assert v < this.num_v;

        return this.neighbors[v].data;
    }

    public int to_original(int v) {
        assert v >= 0;
        assert v < this.num_v;

        return this.to_original[v];
    }

    public int from_original(int v) {
        return this.from_original.get(v);
    }

    public int degree(int v) {
        assert v >= 0;
        assert v < this.num_v;

        return this.neighbors[v].len;
    }

    public int minDeg() {
        return this.minDeg;
    }

    public int maxDeg() {
        return this.maxDeg;
    }

    public int[] fastColor(int start) {
        int[] coloring = new int[this.num_v];
        int[] stack = new int[this.num_v];
        this.fastColor(start,coloring,stack);
        return coloring;
    }

    private int fastColor(int start, int[] coloring, int[] stack) {
        assert start >= 0;
        assert start < this.num_v;

        // coloring for each vertex [0,num_v-1] initialized with -1
        Arrays.fill(coloring,-1);

        // stack for vertex exploration with capacity of num_v and initial len 1 containing only the start vertex
        int stack_len = 1;
        stack[0] = start;

        // boolean buffer for the colors (indices) of the neighbors of a node (not available colors)
        boolean[] n_colors = new boolean[this.num_v];

        // number of colored vertices
        int done = 0;

        // number of colors used
        int color_num = -1;

        // until all vertices are colored
        while (done < this.num_v) {

            // stack not empty
            if (stack_len > 0) {

                // pop stack
                int v = stack[stack_len-1];
                stack_len--;

                // initialize the not available colors all false
                Arrays.fill(n_colors, false);

                // for each neighbour
                for (int neighbor_index = 0; neighbor_index < this.neighbors[v].len; neighbor_index++) {
                    int neighbor = this.neighbors[v].data[neighbor_index];

                    // if already colored, set color as not available
                    if (coloring[neighbor] > -1) {
                        n_colors[coloring[neighbor]] = true;

                    // if not colored and not already on stack (not color value -2) set color to -2, and push neighbor on stack
                    } else if (coloring[neighbor] == -1) {
                        coloring[neighbor] = -2;
                        stack[stack_len] = neighbor;
                        stack_len++;
                    }
                }


                // get lowest color available
                int color = 0;
                for (int i = 0; i < this.num_v; i++) {
                    if (!n_colors[i]) {
                        color = i;
                        break;
                    }
                }

                // set that lowest color and increment number of colored vertices
                coloring[v] = color;
                color_num = Math.max(color_num,color);
                done++;

            // stack is empty
            } else {
                // get first not colored vertex and push on stack
                int v = 0;
                while (coloring[v] > -1) {
                    v++;
                }
                stack[stack_len] = v;
                stack_len++;
            }
        }

        // all vertices have been colored with colors >= 0
        return color_num;
    }

    public int[] slowColor(int start) {
        int[] coloring = new int[this.num_v];
        int[] stack = new int[this.num_v];
        this.slowColor(start,coloring,stack);
        return coloring;
    }

    private int slowColor(int start, int[] stack, int[] coloring) {
        assert start >= 0;
        assert start < this.num_v;

        // coloring for each vertex [0,num_v-1] initialized with -1
        //int[] coloring = new int[this.num_v];
        Arrays.fill(coloring,-1);
        int chromatic = 0;

        // stack for vertex exploration with capacity of num_v and initial len 1 containing only the start vertex
        int stack_len = 1;
        //int[] stack = new int[this.num_v];
        stack[0] = start;

        // boolean buffer for the colors (indices) of the neighbors of a node (not available colors)
        boolean[] n_colors = new boolean[this.num_v];

        // number of colored vertices
        int done = 0;

        // Counting number of colours in neighborhood for each node in Stack
        int[] viscol = new int[this.num_v];
        Arrays.fill(viscol, 0);

        // empty array for updating viscol
        int[] checkviscol = new int[this.num_v];
        int checkvislen = 0;

        // until all vertices are colored
        while (done < this.num_v) {

            // stack not empty
            if (stack_len > 0) {

                //find node with highest viscol
                int v = stack[0];
                int maxvis = viscol[v];
                for (int stack_index = 1; stack_index < stack_len; stack_index++) {
                    if (viscol[stack[stack_index]] > maxvis){
                        maxvis = viscol[stack[stack_index]];
                        v = stack[stack_index];
                        stack[stack_index] = stack[stack_len - 1];
                        stack[stack_len - 1] = v;
                    }
                }


                // pop stack
                v = stack[stack_len-1];
                stack_len--;

                // initialize the not available colors all false
                Arrays.fill(n_colors, false);

                // for each neighbour
                for (int neighbor_index = 0; neighbor_index < this.neighbors[v].len; neighbor_index++) {
                    int neighbor = this.neighbors[v].data[neighbor_index];

                    // if already colored, set color as not available
                    if (coloring[neighbor] > -1) {
                        n_colors[coloring[neighbor]] = true;
                        // Add to temporary List to check for increase in viscol
                        checkviscol[checkvislen] = neighbor;
                        checkvislen++;

                        // if not colored and not already on stack (not color value -2) set color to -2, and push neighbor on stack
                    } else if (coloring[neighbor] == -1) {
                        coloring[neighbor] = -2;
                        stack[stack_len] = neighbor;
                        stack_len++;
                        // increase viscol by one
                        viscol[neighbor] = 1;
                    }
                }


                // get lowest color available
                int color = 0;
                for (int i = 0; i < this.num_v; i++) {
                    if (!n_colors[i]) {
                        color = i;
                        break;
                    }
                }

                // set that lowest color and increment number of colored vertices
                coloring[v] = color;
                if (color + 1 > chromatic) {
                    chromatic = color + 1;
                }
                done++;

                // Update viscol
                while (checkvislen > 0){
                    int node = checkviscol[checkvislen - 1];
                    checkvislen--;
                    for (int neighbor_index = 0; neighbor_index < this.neighbors[node].len; neighbor_index++) {
                        int neighbor = this.neighbors[node].data[neighbor_index];
                        if (coloring[neighbor] == color) {
                            viscol[node] --;
                            neighbor_index = this.neighbors[node].len;
                        }

                    }
                    viscol[node] += 1;   
                }


                // stack is empty
            } else {
                // get first not colored vertex and push on stack
                int v = 0;
                while (coloring[v] > -1) {
                    v++;
                }
                stack[stack_len] = v;
                stack_len++;
            }
        }

        // all vertices have been colored with colors >= 0
        return chromatic;
    }

    public int[] maxClique() {
        long startTimeColoring = System.currentTimeMillis();
        int[] coloring = new int[this.num_v];
        int[] stack = new int[this.num_v];

        int minColor;
        if (this.num_v < MAX_V_SLOW_COLOR) {
            minColor = this.slowColor(this.num_v-1,stack,coloring);
        } else {
            minColor = this.fastColor(this.num_v-1,stack,coloring);
        }

        int num_coloring = 1;
        while (System.currentTimeMillis() - startTimeColoring < MAX_COLORING_TIME) {
            if (num_coloring > MAX_NUM_COLORING) {
                break;
            }
            int vertex = ThreadLocalRandom.current().nextInt(0,  this.num_v);

            int color;
            if (this.num_v < MAX_V_SLOW_COLOR) {
                color = this.slowColor(vertex,stack,coloring);
            } else {
                color = this.fastColor(vertex,stack,coloring);
            }

            minColor = Math.min(minColor,color);
            num_coloring++;
        }

        Arrays.fill(coloring, 0,this.num_v, minColor);

        for (int i = this.num_v-1; i >= 0; i--) {
            if (minColor < this.num_v - i) {
                break;
            }
            coloring[i] = this.num_v - i;
        }

        if (this.hasAdj) {
            return this.maxCliqueADJ(coloring, minColor);
        }
        return this.maxCliqueNoADJ(coloring, minColor);
    }

    public int[] maxCliqueADJ(int[] coloring, int minColor) {
        int[][] s_in = new int[minColor][this.maxDeg+1]; //new bounds
        int s_in_len = 0;

        int[] max = new int[this.num_v];
        int max_len = 0;

        Arrays.fill(max,0,this.num_v,-1);

        int[] current_max = new int[this.num_v];

        for (int i = this.num_v - 1; i >= 0; i--) {
            s_in_len = 0;

            for (int j = i+1; j < this.num_v; j++) {
                s_in[0][s_in_len] = j;
                s_in_len += this.adj[i][j];
            }

            current_max[0] = i;

            int rt = this.maxCliqueADJ(
                    s_in,
                    s_in_len,
                    max,
                    max_len,
                    current_max,
                    1,
                    coloring
            );

            max_len = Math.max(max_len,rt);

            coloring[i] = Math.min(max_len,coloring[i]);
        }

        return max;
    }

    private int maxCliqueADJ(
            int[][] u,
            int u_len,
            int[] max,
            int max_len,
            int[] current_max,
            int current_max_len,
            int[] bound
    ) {
        if (u_len == 0) {
            if (current_max_len > max_len) {
                System.arraycopy(current_max, 0, max, 0, current_max_len);
            }
            return current_max_len;
        }
        int u_index = 0;
        while (u_index < u_len) {
            if (current_max_len + u_len <= max_len) {
                return current_max_len;
            }

            int i = u[current_max_len-1][u_index];
            u_index++;

            if (current_max_len + bound[i] <= max_len) {
                return current_max_len;
            }

            byte[] adj_i = this.adj[i];

            int new_u_len = 0;
            for (int j = u_index; j < u_len; j++) {
                int v = u[current_max_len-1][j];
                u[current_max_len][new_u_len] = v;
                new_u_len += adj_i[v];
            }

            current_max[current_max_len] = i;

            int rt = this.maxCliqueADJ(
                    u,
                    new_u_len,
                    max,
                    max_len,
                    current_max,
                    current_max_len + 1,
                    bound
            );

            if (rt > max_len) {
                return rt;
            }
        }

        return max_len;
    }

    // NO ADJ MATRIX

    public int[] maxCliqueNoADJ(int[] coloring, int minColor) {
        int[][] s_in = new int[minColor][this.maxDeg+1]; //new bounds
        int s_in_len = 0;

        int[] max = new int[this.num_v];
        int max_len = 0;

        Arrays.fill(max,0,this.num_v,-1);

        int[] current_max = new int[this.num_v];

        for (int i = this.num_v - 1; i >= 0; i--) {
            s_in_len = 0;

            for (int j = i+1; j < this.num_v; j++) {
                s_in[0][s_in_len] = j;
                s_in_len += (Arrays.binarySearch(this.neighbors[i].data, j) >= 0 ? 1 : 0);
            }

            current_max[0] = i;

            int rt = this.maxCliqueNoADJ(
                    s_in,
                    s_in_len,
                    max,
                    max_len,
                    current_max,
                    1,
                    coloring
            );

            max_len = Math.max(max_len,rt);

            coloring[i] = Math.min(max_len,coloring[i]);
        }

        return max;
    }

    private int maxCliqueNoADJ(
            int[][] u,
            int u_len,
            int[] max,
            int max_len,
            int[] current_max,
            int current_max_len,
            int[] bound
    ) {
        if (u_len == 0) {
            if (current_max_len > max_len) {
                System.arraycopy(current_max, 0, max, 0, current_max_len);
            }
            return current_max_len;
        }
        int u_index = 0;
        while (u_index < u_len) {
            if (current_max_len + u_len <= max_len) {
                return current_max_len;
            }

            int i = u[current_max_len-1][u_index];
            u_index++;

            if (current_max_len + bound[i] <= max_len) {
                return current_max_len;
            }

            int new_u_len = 0;
            for (int j = u_index; j < u_len; j++) {
                int v = u[current_max_len-1][j];
                u[current_max_len][new_u_len] = v;
                new_u_len += (Arrays.binarySearch(this.neighbors[i].data, v) >= 0 ? 1 : 0);
            }

            current_max[current_max_len] = i;

            int rt = this.maxCliqueNoADJ(
                    u,
                    new_u_len,
                    max,
                    max_len,
                    current_max,
                    current_max_len + 1,
                    bound
            );

            if (rt > max_len) {
                return rt;
            }
        }

        return max_len;
    }
}
