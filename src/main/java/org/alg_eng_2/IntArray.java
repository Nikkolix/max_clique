package org.alg_eng_2;

public class IntArray {
    public final int[] data;
    public int len;

    IntArray(int size) {
        this.data = new int[size];
        this.len = 0;
    }

    IntArray(int[] data)  {
        this.data = data;
        this.len = data.length;
    }

    public void add(int x) {
        this.data[this.len] = x;
        this.len++;
    }
}
