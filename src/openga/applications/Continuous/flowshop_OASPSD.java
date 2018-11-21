package openga.applications.Continuous;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.concurrent.CountDownLatch;
import openga.chromosomes.*;
import openga.operator.selection.*;
import openga.operator.crossover.*;
import openga.operator.mutation.*;
import openga.ObjectiveFunctions.*;
import openga.MainProgram.*;
import openga.Fitness.*;
import openga.util.fileWrite1;
import org.moeaframework.Executor;
import org.moeaframework.analysis.plot.Plot;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Population;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.Permutation;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;

/**
 * <p>
 * Title: The OpenGA project which is to build general framework of Genetic
 * algorithm.</p>
 * <p>
 * Description: OAS Flowshop Problem.</p>
 * <p>
 * Copyright: Copyright (c) 2018</p>
 * <p>
 * Company: Cheng-Shiu University</p>
 *
 * @author
 * @version 1.1
 */
public class flowshop_OASPSD extends AbstractProblem {

    /**
     * *
     * Basic variables of GAs.
     */
    int numberOfObjs = 1;
    double b = 0.1;
    static int variabletotal = 0;
    static int Objectivestotal = 0;

    populationI Population;
    SelectI Selection;
    CrossoverI Crossover, Crossover2;
    MutationI Mutation, Mutation2;
    ObjFunctionPFSSOAWT_PSDI[] ObjectiveFunction;
    FitnessI Fitness;
    MainI GaMain;

    /**
     * Parameters of the GA
     */
    int generations, length, initPopSize, fixPopSize;
    double crossoverRate, mutationRate;
    boolean[] objectiveMinimization; //true is minimum problem.
    boolean encodeType;  //binary of realize code
    int seed = 12345;
    int counter = 0;

//  Instance Data
    int piTotal;
    int machineTotal;
    int[] fristProfit;    //  revenue of order
    int[] di;        //  due-date
    double[] wi;     //  tardiness penalty weight
    int[][] processingTime;

    //Results
    double bestObjectiveValues[];
    populationI solutions;

    public int DEFAULT_generations = 1000,
            DEFAULT_PopSize = 100,
            DEFAULT_initPopSize = 100;

    public double DEFAULT_crossoverRate = 0.9,
            DEFAULT_mutationRate = 0.2,
            elitism = 0.2;     //the percentage of elite chromosomes

    int instance;
    double originalPoint[];
    double coordinates[][];
    double distanceMatrix[][];
    String instanceName = "";
    boolean applyLocalSearch = true;
    CountDownLatch latch;

    /**
     * The method is to modify the default value.
     */
    public void setParameter(double crossoverRate, double mutationRate, int counter, double elitism, int generation,
            int length, String instanceName, int piTotal, int machineTotal, int[] fristProfit, int[] di, double[] wi, int[][] processingTime, double b) {
        this.DEFAULT_crossoverRate = crossoverRate;
        this.DEFAULT_mutationRate = mutationRate;
        this.counter = counter;
        this.elitism = elitism;
        this.DEFAULT_generations = generation;
        this.instanceName = instanceName;
        this.length = length;

        this.piTotal = piTotal;
        this.machineTotal = machineTotal;
        this.fristProfit = fristProfit;
        this.di = di;
        this.wi = wi;
        this.processingTime = processingTime;
        this.b = b;
    }

    public void initiateVars() throws IOException {
        GaMain = new singleThreadGA();//singleThreadGAwithMultipleCrossover singleThreadGA adaptiveGA
        Population = new population();
        Selection = new binaryTournament();
        Crossover = new CyclingCrossoverP(); //twoPointCrossover2()  CyclingCrossoverP multiParentsCrossover()
        Crossover2 = new PMX();
        Mutation = new swapMutation();//shiftMutation
        Mutation2 = new shiftMutation();//inverseMutation
        ObjectiveFunction = new ObjFunctionPFSSOAWT_PSDI[numberOfObjs];
        ObjectiveFunction[0] = new ObjectiveFunctionMakespanFlowShop_OASPSD();//the first objective
        Fitness = new singleObjectiveFitness();//singleObjectiveFitness singleObjectiveFitnessByNormalize
        objectiveMinimization = new boolean[numberOfObjs];
        objectiveMinimization[0] = false;
        encodeType = true;

        ObjectiveFunction[0].setOASData(this.piTotal, this.machineTotal, this.fristProfit, this.di, this.wi, this.processingTime, this.b);

        //set the data to the GA main program.
        GaMain.setData(Population, Selection, Crossover, Mutation, ObjectiveFunction,
                Fitness, DEFAULT_generations, DEFAULT_initPopSize, DEFAULT_PopSize,
                length, DEFAULT_crossoverRate, DEFAULT_mutationRate, objectiveMinimization,
                numberOfObjs, encodeType, elitism);
        GaMain.setSecondaryCrossoverOperator(Crossover2, false);
        GaMain.setSecondaryMutationOperator(Mutation2, true);
    }

    public String start() {
        openga.util.timeClock timeClock1 = new openga.util.timeClock();
        timeClock1.start();
        GaMain.startGA();
        timeClock1.end();
        String implementResult = instanceName + "\t" + DEFAULT_crossoverRate + "\t" + DEFAULT_mutationRate + "\t"
                + elitism + "\t" + GaMain.getArchieve().getSingleChromosome(0).getObjValue()[0]
                + "\t " + timeClock1.getExecutionTime() / 1000.0 + "\n";
        writeFile("flowshop_OASPSD", implementResult);
        System.out.print(implementResult);

        return implementResult;
    }

    /**
     * Write the data into text file.
     */
    void writeFile(String fileName, String _result) {
        fileWrite1 writeResult = new fileWrite1();
        writeResult.writeToFile(_result, fileName + ".txt");
        writeResult.run();
//    Thread thread1 = new Thread(writeResult);
//    thread1.run();
    }

    public flowshop_OASPSD() {

        super(variabletotal, Objectivestotal);
    }

    @Override
    public void evaluate(Solution solution) {
        //double[] x = EncodingUtils.getReal(solution);
        int[] x =  ((Permutation) solution.getVariable(0)).toArray();
        double[] f = new double[numberOfObjectives];

        int k = numberOfVariables - numberOfObjectives + 1;

        double g = 0.0;
        for (int i = numberOfVariables - k; i < numberOfVariables; i++) {
            g += Math.pow(x[i] - 0.5, 2.0);
        }

        for (int i = 0; i < numberOfObjectives; i++) {
            f[i] = 1.0 + g;

            for (int j = 0; j < numberOfObjectives - i - 1; j++) {
                f[i] *= Math.cos(0.5 * Math.PI * x[j]);
            }

            if (i != 0) {
                f[i] *= Math.sin(0.5 * Math.PI * x[numberOfObjectives - i - 1]);
            }
        }
        solution.setObjectives(f);
//        double x = ((Permutation) solution.getVariable(0)).get(0);
//        double y = ((Permutation) solution.getVariable(0)).get(0);
//        double f1 = Math.pow(x - 2.0, 2.0) + Math.pow(y - 1.0, 2.0) + 2.0;
//        double f2 = 9.0 * x - Math.pow(y - 1.0, 2.0);
//        double c1 = Math.pow(x, 2.0) + Math.pow(y, 2.0) - 225.0;
//        double c2 = x - 3.0 * y + 10.0;
//
//        solution.setObjective(0, f1);
//        solution.setObjective(1, f2);
    }

    @Override
    public Solution newSolution() {
        int[] ints = new int[6];
        ints[0] = 0;
        ints[1] = 1;
        ints[2] = 2;
        ints[3] = 3;
        ints[4] = 4;
        ints[5] = 5;

        Solution solution = new Solution(getNumberOfVariables(),
                getNumberOfObjectives());
        for (int i = 0; i < variabletotal; i++) {
            for (int j = 1; j <= Objectivestotal; j++) {
                Permutation test = new Permutation(ints);
                solution.setVariable(i, test);
            }
        }
        return solution;
    }

    public static void main(String[] args) throws IOException, ParseException {

        String data = "instances/PFSS-OAWT-Data/p/", fileName;
        File f = new File(data);
        String[] fn = f.list();
        Population manualSolutions = new Population();

        double crossoverRate[], mutationRate[];
        double b = 0.1;   //0.1
        crossoverRate = new double[]{0.5, 0.9};
        mutationRate = new double[]{0.1, 0.5};
        int counter = 0;
        double elitism[] = new double[]{0.2};
        int generations[] = new int[]{1000};
        int instanceReplications = 1;
        Objectivestotal = 2;
        variabletotal = 2;

        for (int filelist = 0; filelist < fn.length; filelist++) {
            fileName = fn[filelist];
            if (fileName.substring(fileName.indexOf("_") + 1, fileName.indexOf("_") + 2).equals("0")) {
                openga.applications.data.readPFSSOAWT_flowshop fs = new openga.applications.data.readPFSSOAWT_flowshop();
                fs.setData(data, fileName);
                fs.readfile();
                int Length = fs.getPiTotal();
                Problem problem = new flowshop_OASPSD(); // put your problem here
                // load your solutions from file
                String line = null;
                BufferedReader reader = null;

                for (int k = 0; k < instanceReplications; k++) {
                    for (int m = 0; m < crossoverRate.length; m++) {
                        for (int n = 0; n < mutationRate.length; n++) {
                            for (int o = 0; o < elitism.length; o++) {
                                Solution solution = problem.newSolution();
//                                flowshop_OASPSD flowshop1 = new flowshop_OASPSD();
//                                flowshop1.setParameter(crossoverRate[m], mutationRate[n], counter, elitism[o], generations[0], length, fs.getfileName(),
//                                        fs.getPiTotal(), fs.getMachineTotal(), fs.getprofit(), fs.getdi(), fs.getwi(), fs.getprocessingTime(), b);
//                                flowshop1.initiateVars();
//                                flowshop1.start();
//                                counter++;
                                int[] ints = new int[6];
                                ints[0] = 0;
                                ints[1] = 1;
                                ints[2] = 2;
                                ints[3] = 3;
                                ints[4] = 4;
                                ints[5] = 5;
//                                for (int i = 0; i < problem.getNumberOfVariables(); i++) {
//                                    EncodingUtils.setPermutation(solution.getVariable(i), ints);
//                                }
                                EncodingUtils.setPermutation(solution.getVariable(0), ints);
                                manualSolutions.add(solution);

                                break;
                            }
                        }
                    }//end for
                }
                NondominatedPopulation nsgaResults = new Executor()
                        .withProblem(problem) // put your problem here
                        .withAlgorithm("NSGAII")
                        //.withProperty("populationSize", 30)
                        .withMaxEvaluations(100)
                        .run();
                for (Solution solution : nsgaResults) {
                    System.out.format("%.4f      %.4f%n",
                            solution.getObjective(0),
                            solution.getObjective(1));
                }
//                new Plot()
//                        //.add("External", manualSolutions)
//                        .add("NSGA-II", nsgaResults)
//                        .show();
//                try {
//                    reader = new BufferedReader(new FileReader(file));
//
//                    while ((line = reader.readLine()) != null) {
//                        Solution solution = problem.newSolution();
//                        String[] tokens = line.split("\\s+"); // split by whitespace
//                        int[] ints = new int[100];
//                        int icout = 0;
//                        for (int i = 99; i >= 0; i--) {
//                            ints[icout] = i;
//                            icout++;
//                        }
//                        for (int i = 0; i < problem.getNumberOfVariables(); i++) {
//                            EncodingUtils.setPermutation(solution.getVariable(i), ints);
//                        }
//                        manualSolutions.add(solution);
//                    }
//                } finally {
//                    if (reader != null) {
//                        reader.close();
//                    }
//                }

//                new Plot()
//                        //.add("External", manualSolutions)
//                        .add("NSGA-II", nsgaResults)
//                        .show();
//                Initialization initialization = new InjectedInitialization(
//                        problem,
//                        100);
//
//                TournamentSelection selection = new TournamentSelection(2,
//                        new ChainedComparator(
//                                new ParetoDominanceComparator(),
//                                new CrowdingComparator()));
//
//                Variation variation = new GAVariation(
//                        new SBX(1.0, 25.0),
//                        new PM(1.0 / problem.getNumberOfVariables(), 30.0));
//
//                NSGAII algorithm = new NSGAII(
//                        problem,
//                        new NondominatedSortingPopulation(),
//                        null, // no archive
//                        selection,
//                        variation,
//                        initialization);
//
//                while (algorithm.getNumberOfEvaluations() < 10000) {
//                    algorithm.step();
//                }
//
//                Population intermediateResult = algorithm.getPopulation();
//               
//                algorithm = new NSGAII(
//                        problem,
//                        new NondominatedSortingPopulation(),
//                        null, // no archive
//                        selection,
//                        variation,
//                        initialization);
//
//                while (algorithm.getNumberOfEvaluations() < 10000) {
//                    algorithm.step();
//                }
//NondominatedPopulation finalResult = algorithm.getResult();
//                int length = fs.getPiTotal();
//                int innerLoopSize = instanceReplications * crossoverRate.length * mutationRate.length * elitism.length;
//                CountDownLatch latch = new CountDownLatch(innerLoopSize);
//
//                for (int k = 0; k < instanceReplications; k++) {
//                    for (int m = 0; m < crossoverRate.length; m++) {
//                        for (int n = 0; n < mutationRate.length; n++) {
//                            for (int o = 0; o < elitism.length; o++) {
//                                flowshop_OASPSD flowshop1 = new flowshop_OASPSD();
//                                flowshop1.setParameter(crossoverRate[m], mutationRate[n], counter, elitism[o], generations[0], length, fs.getfileName(),
//                                        fs.getPiTotal(), fs.getMachineTotal(), fs.getprofit(), fs.getdi(), fs.getwi(), fs.getprocessingTime(), b, latch);
//                                flowshop1.initiateVars();
//                                flowshop1.start();
//                                counter++;
//                            }
//                        }
//                    }//end for
//                }
//                try {
//                    //Wait the all works are done. Then we process next instance.
//                    latch.await();
//                } catch (InterruptedException E) {
//                    E.printStackTrace();
//                }
            }
            break;
        }
    }
}
//
