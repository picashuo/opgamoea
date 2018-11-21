package openga.applications.singleMachineProblem;
import openga.chromosomes.*;
import openga.operator.selection.*;
import openga.operator.crossover.*;
import openga.operator.localSearch.*;
import openga.operator.mutation.*;
import openga.ObjectiveFunctions.*;
import openga.MainProgram.*;
import openga.ObjectiveFunctions.*;
import openga.Fitness.*;
import openga.util.printClass;
import openga.util.fileWrite1;
import openga.applications.data.*;
import openga.applications.singleMachine;

/**
 * <p>Title: The OpenGA project which is to build general framework of Genetic algorithm.</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2008</p>
 * <p>Company: Yuan-Ze University</p>
 * @author Chen, Shih-Hsin
 * @version 1.0
 */

public class singleMachineGLS extends singleMachine{
  public singleMachineGLS() {
  }
  /**
   * Local Search Operator
   */
  boolean applyLocalSearch = true;
  localSearchI localSearch1;
  int maxNeighborhood = 5;  //A default value of the maximum neighbors to search.

  public void setLocalSearchData(boolean applyLocalSearch, int maxNeighborhood){
    this.applyLocalSearch = applyLocalSearch;
    this.maxNeighborhood = maxNeighborhood;
  }

  public void initiateVars(){
    GaMain     = new singleThreadGAwithLocalSearch();
    Population = new population();
    Selection  = new binaryTournament();
    Crossover  = new twoPointCrossover2();
    Mutation   = new swapMutation();
    ObjectiveFunction = new ObjectiveFunctionScheduleI[numberOfObjs];
    ObjectiveFunction[0] = new ObjectiveEarlinessTardinessPenalty();
    Fitness    = new singleObjectiveFitness();
    objectiveMinimization = new boolean[numberOfObjs];
    objectiveMinimization[0] = true;
    encodeType = true;

    //set schedule data to the objectives
    ObjectiveFunction[0].setScheduleData(dueDay, processingTime, numberOfMachines);

    DEFAULT_generations = totalSolnsToExamine/(DEFAULT_PopSize);
    //set the data to the GA main program.
    GaMain.setData(Population, Selection, Crossover, Mutation, ObjectiveFunction, Fitness, DEFAULT_generations, DEFAULT_PopSize, DEFAULT_PopSize,
                   numberOfJob, DEFAULT_crossoverRate, DEFAULT_mutationRate, objectiveMinimization, numberOfObjs, encodeType, elitism);
    GaMain.setSecondaryCrossoverOperator(Crossover2, false);
    GaMain.setSecondaryMutationOperator(Mutation2, false);

    localSearch1 = new localSearchBy2Opt();//localSearchBy2Opt localSearchBySwap
    GaMain.setLocalSearchOperator(localSearch1, applyLocalSearch, maxNeighborhood);
  }

  public void startMain(){
    openga.util.timeClock timeClock1 = new openga.util.timeClock();
    timeClock1.start();
    GaMain.startGA();
    timeClock1.end();
    //to output the implementation result.
    String implementResult = "";
    int bestInd = getBestSolnIndex(GaMain.getArchieve());
    implementResult = fileName+"\t" +maxNeighborhood+"\t"
        +GaMain.getArchieve().getSingleChromosome(bestInd).getObjValue()[0]
        +"\t"+timeClock1.getExecutionTime()/1000.0+"\n";
    writeFile("singleMachineGLS_SKS_20080428", implementResult);
    System.out.print(implementResult);
  }

  /**
   * For single objective problem
   * @param arch1
   * @return
   */
  public int getBestSolnIndex(populationI arch1){
    int index = 0;
    double bestobj = Double.MAX_VALUE;
    for(int k = 0 ; k < GaMain.getArchieve().getPopulationSize() ; k ++ ){
      if(bestobj > GaMain.getArchieve().getObjectiveValues(k)[0]){
        bestobj = GaMain.getArchieve().getObjectiveValues(k)[0];
        index = k;
      }
    }
    return index;
  }

  public static void main(String[] args) {
    System.out.println("singleMachineGLS_SKS_20080428");
    int jobSets[] = new int[]{20, 30, 40, 50, 60, 90};//20, 30, 40, 50, 60, 90, 100, 200//20, 40, 60, 80
    int counter = 0;
    int repeatExperiments = 30;

    int popSize[] = new int[]{100};//50, 100, 155, 210 [100]
    double crossoverRate[] = new double[]{0.9},//0.6, 0.9 {0.9}
           mutationRate [] = new double[]{1.0},//0.1, 0.5 {0.5}
           elitism = 0.2;

    //local search
    boolean applyLocalSearch = true;
    int maxNeighborhood[] = new int[]{3};  //A default value of the maximum neighbors to search.

    //EDA parameters.
    double[] lamda = new double[]{0.1}; //learning rate{0.1, 0.5, 0.9}
    int numberOfCrossoverTournament[] = new int[]{2};//{1, 2, 4}
    int numberOfMutationTournament[] = new int[]{2};//{1, 2, 4}
    int startingGenDividen[] = new int[]{10};//{2, 4}{4}

    for (int j = 0; j < jobSets.length; j++) { //jobSets.length
      for (int k = 0; k < 49; k++) {
        /*
        if ( (jobSets[j] <= 50 &&
              (k == 0 || k == 3 || k == 6 || k == 21 || k == 24 || k == 27 ||
               k == 42 || k == 45 || k == 48)) || (jobSets[j] > 50 && k < 9)) {
            */
          if((jobSets[j] <= 50 && (k != 0 && k != 3 && k != 6 && k != 21 && k != 24 && k != 27 && k != 42 && k != 45 && k != 48)) ||  (jobSets[j] > 50 && k < 9)){
          for (int q = 0; q < maxNeighborhood.length; q++) {
            openga.applications.data.singleMachine readSingleMachineData1 = new
                openga.applications.data.singleMachine();
            int numberOfJobs = jobSets[j];
            String fileName = readSingleMachineData1.getFileName(
                numberOfJobs, k);
            System.out.print(fileName + "\t");
            readSingleMachineData1.setData("sks/" + fileName + ".txt");
            readSingleMachineData1.getDataFromFile();
            int dueDate[] = readSingleMachineData1.getDueDate();
            int processingTime[] = readSingleMachineData1.getPtime();

            for (int i = 0; i < repeatExperiments; i++) {
              System.out.println("Combinations: " + counter);
              singleMachineGLS singleMachine1 = new singleMachineGLS();
              singleMachine1.setData(numberOfJobs, dueDate, processingTime,
                                     fileName);
              singleMachine1.setLocalSearchData(applyLocalSearch,
                                                maxNeighborhood[q]);
              singleMachine1.initiateVars();
              singleMachine1.startMain();
              counter++;
            }
          }
        }
      }
    }

  }

}