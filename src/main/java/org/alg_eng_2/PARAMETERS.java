package org.alg_eng_2;

public final class PARAMETERS {
    public static final boolean DEBUG_MODE = false; // set this to false before submission
    public static final boolean PRINT_IF_KERNEL_WAS_APPLIED = false; // can be used to analyze instances
    public static final int MIN_VERTICES_DYN = 100_000;
    public static final int MIN_V_MASSIVE = Integer.MAX_VALUE;
    public static final int MIN_E_MASSIVE = Integer.MAX_VALUE;
    public static final long MAX_REDUCE_TIME = 5 * 1000;
    public static final long MAX_COLORING_TIME = 3 * 1000;
    public static final int MAX_NUM_COLORING = 32;
    public static final int MAX_V_ADJ_STATIC = 100_000;
    public static final int MAX_V_FOR_EDGE_REDUCTION = 1_000;
    public static final int MAX_V_SLOW_COLOR = 100_000;
}
