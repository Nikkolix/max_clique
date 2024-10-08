package org.alg_eng_2;

import java.io.File;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Scanner;
import java.util.regex.Pattern;

import static org.alg_eng_2.PARAMETERS.MIN_E_MASSIVE;
import static org.alg_eng_2.PARAMETERS.MIN_V_MASSIVE;

public class GraphBuilder {
    private final DynGraph g;

    public GraphBuilder(String filepath) throws Exception {
        Scanner scanner;
        if (filepath.equals("-")) {
            // access stdin
            scanner = new Scanner(System.in);
        } else {
            // access file
            File file = new File(filepath);
            Path path = file.toPath();
            InputStream inputStream = Files.newInputStream(path);
            scanner = new Scanner(inputStream);
        }
        // process the header
        String header = scanner.nextLine().trim();
        if (Character.isDigit(header.charAt(1))) {
            header = "# " + header.substring(1);
        }
        Pattern space = Pattern.compile("\\s+");
        String[] headerFields = space.split(header);

        int v = Integer.parseInt(headerFields[1].trim());
        int e = Integer.parseInt(headerFields[2].trim());

        IntArray edgesLeft = new IntArray(e);
        IntArray edgesRight = new IntArray(e);

        // process the remaining lines
        Pattern edge = Pattern.compile("^[0-9]+\\s+[0-9]+$");
        while (scanner.hasNextLine()) {
            String line = scanner.nextLine().trim();
            if (edge.matcher(line).matches()) {
                String[] lineFields = space.split(line);
                int left = Integer.parseInt(lineFields[0].trim());
                int right = Integer.parseInt(lineFields[1].trim());

                assert left != right : "looping edge";

                edgesLeft.add(left);
                edgesRight.add(right);
            }
        }


        if (v > MIN_V_MASSIVE || e > MIN_E_MASSIVE) {
            this.g = new MassiveGraph(v,e,edgesLeft,edgesRight);
        } else {
            this.g = new HashGraph(v,e,edgesLeft,edgesRight);
        }

        scanner.close();
    }

    public DynGraph build() {
        return g;
    }
}
