/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package moea;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextArea;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.EvolutionaryAlgorithm;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Settings;
import org.moeaframework.core.Solution;
import org.moeaframework.core.spi.AlgorithmFactory;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.examples.ga.tsplib.TSP2OptHeuristic;
import org.moeaframework.examples.ga.tsplib.TSPExample;
import static org.moeaframework.examples.ga.tsplib.TSPExample.solve;
import static org.moeaframework.examples.ga.tsplib.TSPExample.toTour;
import org.moeaframework.examples.ga.tsplib.TSPInstance;
import org.moeaframework.examples.ga.tsplib.TSPPanel;
import org.moeaframework.examples.ga.tsplib.Tour;
import org.moeaframework.util.Vector;
import org.moeaframework.util.io.CommentedLineReader;

/**
 *
 * @author User
 */
public class opgaproblem implements Problem {

    private static final Color lightGray = new Color(128, 128, 128, 64);

    private TSPInstance instance;

    private TSP2OptHeuristic heuristic;

    /**
     * The number of sacks.
     */
    private int nsacks;

    /**
     * The number of items.
     */
    private int nitems;

    /**
     * Entry {@code profit[i][j]} is the profit from including item {@code j} in
     * sack {@code i}.
     */
    private int[][] profit;

    /**
     * Entry {@code weight[i][j]} is the weight incurred from including item
     * {@code j} in sack {@code i}.
     */
    private int[][] weight;

    /**
     * Entry {@code capacity[i]} is the weight capacity of sack {@code i}.
     */
    private int[] capacity;

    public opgaproblem(TSPInstance instance) {
        //super(1, 1);
        this.instance = instance;

        heuristic = new TSP2OptHeuristic(instance);
    }

    /**
     * Constructs a multiobjective 0/1 knapsack problem instance loaded from the
     * specified file.
     *
     * @param file the file containing the knapsack problem instance
     * @throws IOException if an I/O error occurred
     */
    public opgaproblem(File file) throws IOException {
        this(new FileReader(file));
    }

    /**
     * Constructs a multiobjective 0/1 knapsack problem instance loaded from the
     * specified input stream.
     *
     * @param inputStream the input stream containing the knapsack problem
     * instance
     * @throws IOException if an I/O error occurred
     */
    public opgaproblem(InputStream inputStream) throws IOException {
        this(new InputStreamReader(inputStream));
    }

    /**
     * Constructs a multiobjective 0/1 knapsack problem instance loaded from the
     * specified reader.
     *
     * @param reader the reader containing the knapsack problem instance
     * @throws IOException if an I/O error occurred
     */
    public opgaproblem(Reader reader) throws IOException {
        super();

        load(reader);
    }

    /**
     * Loads the knapsack problem instance from the specified reader.
     *
     * @param reader the file containing the knapsack problem instance
     * @throws IOException if an I/O error occurred
     */
    private void load(Reader reader) throws IOException {
        CommentedLineReader lineReader = null;
        String line = null;
        Matcher matcher = null;
        try {
            lineReader = new CommentedLineReader(reader);
            line = lineReader.readLine(); // the problem specification line

        } finally {
            if (lineReader != null) {
                lineReader.close();
            }
        }
    }

    @Override
    public void evaluate(Solution solution) {
        boolean[] d = EncodingUtils.getBinary(solution.getVariable(0));
        double[] f = new double[nsacks];
        double[] g = new double[nsacks];

        // calculate the profits and weights for the knapsacks
        for (int i = 0; i < nitems; i++) {
            if (d[i]) {
                for (int j = 0; j < nsacks; j++) {
                    f[j] += profit[j][i];
                    g[j] += weight[j][i];
                }
            }
        }

        // check if any weights exceed the capacities
        for (int j = 0; j < nsacks; j++) {
            if (g[j] <= capacity[j]) {
                g[j] = 0.0;
            } else {
                g[j] = g[j] - capacity[j];
            }
        }

        // negate the objectives since Knapsack is maximization
        solution.setObjectives(Vector.negate(f));
        solution.setConstraints(g);
    }

    @Override
    public String getName() {
        return "Knapsack";
    }

    @Override
    public int getNumberOfConstraints() {
        return nsacks;
    }

    @Override
    public int getNumberOfObjectives() {
        return nsacks;
    }

    @Override
    public int getNumberOfVariables() {
        return 1;
    }

    @Override
    public Solution newSolution() {
        Solution solution = new Solution(1, 1);
        solution.setVariable(0, EncodingUtils.newPermutation(
                instance.getDimension()));
        return solution;
//        Solution solution = new Solution(1, nsacks, nsacks);
//        solution.setVariable(0, EncodingUtils.newBinary(nitems));
//        return solution;
    }

    @Override
    public void close() {
        //do nothing
    }

    /**
     * Solves this TSPLIB instance while displaying a GUI showing the
     * optimization progress.
     *
     * @param instance the TSPLIB instance to solve
     */
    public static void solve(TSPInstance instance) {
        TSPPanel panel = new TSPPanel(instance);
        panel.setAutoRepaint(false);

        // create other components on the display
        StringBuilder progress = new StringBuilder();
        JTextArea progressText = new JTextArea();

        JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        splitPane.setTopComponent(panel);
        splitPane.setBottomComponent(new JScrollPane(progressText));
        splitPane.setDividerLocation(300);
        splitPane.setResizeWeight(1.0);

        // display the panel on a window
        JFrame frame = new JFrame(instance.getName());
        frame.getContentPane().setLayout(new BorderLayout());
        frame.getContentPane().add(splitPane, BorderLayout.CENTER);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.setSize(500, 400);
        frame.setLocationRelativeTo(null);
        frame.setIconImages(Settings.getIconImages());
        frame.setVisible(true);

        // create the optimization problem and evolutionary algorithm
        Problem problem = new TSPExample.TSPProblem(instance);

        Properties properties = new Properties();
        properties.setProperty("swap.rate", "0.7");
        properties.setProperty("insertion.rate", "0.9");
        properties.setProperty("pmx.rate", "0.4");

        Algorithm algorithm = AlgorithmFactory.getInstance().getAlgorithm(
                "GA", properties, problem);

        int iteration = 0;

        // now run the evolutionary algorithm
        while (frame.isVisible()) {
            algorithm.step();
            iteration++;

            // clear existing tours in display
            panel.clearTours();

            // display population with light gray lines
            if (algorithm instanceof EvolutionaryAlgorithm) {
                EvolutionaryAlgorithm ea = (EvolutionaryAlgorithm) algorithm;

                for (Solution solution : ea.getPopulation()) {
                    panel.displayTour(toTour(solution), lightGray);
                }
            }

            // display current optimal solutions with red line
            Tour best = toTour(algorithm.getResult().get(0));
            panel.displayTour(best, Color.RED, new BasicStroke(2.0f));
            progress.insert(0, "Iteration " + iteration + ": "
                    + best.distance(instance) + "\n");
            progressText.setText(progress.toString());

            // repaint the TSP display
            panel.repaint();
        }
    }

    /**
     * Runs the example TSP optimization problem.
     *
     * @param file the file containing the TSPLIB instance
     * @throws IOException if an I/O error occurred
     */
    public static void solve(File file) throws IOException {
        solve(new TSPInstance(file));
    }

    /**
     * Runs the example TSP optimization problem.
     *
     * @param reader the reader containing the TSPLIB instance
     * @throws IOException if an I/O error occurred
     */
    public static void solve(Reader reader) throws IOException {
        solve(new TSPInstance(reader));
    }

    /**
     * Runs the example TSP optimization problem.
     *
     * @param stream the stream containing the TSPLIB instance
     * @throws IOException if an I/O error occurred
     */
    public static void solve(InputStream stream) throws IOException {
        solve(new TSPInstance(new InputStreamReader(stream)));
    }

}
