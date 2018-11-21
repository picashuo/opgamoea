package openga.applications.singleMachineProblem;

import openga.chromosomes.*;
import openga.operator.selection.*;
import openga.operator.crossover.*;
import openga.operator.mutation.*;
import openga.ObjectiveFunctions.*;
import openga.MainProgram.*;
import openga.Fitness.*;
import openga.operator.localSearch.*;
/**
 * <p>Title: Self-Guided Inver-Over Operator for the single machine scheduling problem with setup and dynamic arrival
 * time in a common due day environment.</p>
 * <p>Description: Self-Guided Inver-Over Operator is constructed under the eSGGA. So it considers the position and interaction
 * of the variables.</p>
 * <p>Copyright: Copyright (c) 2015</p>
 * <p>Company: Cheng Shiu University University</p>
 * @author Chen, Shih-Hsin
 * @version 1.0
 */
public class singleMachineSetupDynamicArrivalSGIO extends singleMachineEDA2 {

    public singleMachineSetupDynamicArrivalSGIO() {
        
    }
    
    int processingTime[];
    int setupTime[][];
    int dynamicArrivalTime[];
    
    ObjectiveFunctionMatrixPTimeScheduleI ObjectiveFunction[];
    
    /**
     * Local Search Operator
     */
    boolean applyLocalSearch = true;
    localSearchI localSearch1;
    int maxNeighborhood = 3;  //A default value of the maximum neighbors to search.

    public void setData(int numberOfJobs, int processingTime[], int setupTime[][], 
            int dynamicArrivalTime[], String fileName){
      this.numberOfJob = numberOfJobs;      
      this.processingTime = processingTime;
      this.setupTime = setupTime;
      this.dynamicArrivalTime = dynamicArrivalTime;
      this.fileName = fileName;
    }    

    public void initiateVars() {
        GaMain = new singleThreadGAwithEDA2();//singleThreadGA singleThreadGAwithSecondFront singleThreadGAwithMultipleCrossover adaptiveGA
        Population = new population();
        Selection = new binaryTournament();//binaryTournament
        Crossover = new inverOverOperatorEDA2();//twoPointCrossover2 oneByOneChromosomeCrossover twoPointCrossover2withAdpative twoPointCrossover2withAdpativeThreshold
        Mutation = new dummyMutationEDA2();//shiftMutation shiftMutationWithAdaptive shiftMutationWithAdaptiveThreshold
        ObjectiveFunction = new ObjectiveFunctionMatrixPTimeScheduleI[numberOfObjs];
        ObjectiveFunction[0] = new ObjectiveETPenaltyDynamicArrval();
        Fitness = new singleObjectiveFitness();
        objectiveMinimization = new boolean[numberOfObjs];
        objectiveMinimization[0] = true;
        encodeType = true;
        //clone1 = new solutionVectorCloneWithMutation();//swap mutation
        //GaMain.setCloneOperatpr(clone1, true);
        //set schedule data to the objectives
        ObjectiveFunction[0].setScheduleData(processingTime, numberOfMachines);
        ObjectiveFunction[0].setScheduleData(setupTime, numberOfMachines);//We pass the setup time here.
        ((dynamicArrivalTimeI)ObjectiveFunction[0]).setDynamicArrivalTime(dynamicArrivalTime);
        DEFAULT_PopSize = 100;
        //System.out.println(DEFAULT_PopSize);
        //System.exit(0);
        DEFAULT_generations = totalSolnsToExamine / (DEFAULT_PopSize);
        //set the data to the GA main program.
        GaMain.setData(Population, Selection, Crossover, Mutation, ObjectiveFunction, Fitness, DEFAULT_generations, DEFAULT_PopSize, DEFAULT_PopSize,
                numberOfJob, DEFAULT_crossoverRate, DEFAULT_mutationRate, objectiveMinimization, numberOfObjs, encodeType, elitism);
        GaMain.setSecondaryCrossoverOperator(Crossover2, false);
        GaMain.setSecondaryMutationOperator(Mutation2, false);
        GaMain.setEDAinfo(lamda, beta, numberOfCrossoverTournament, numberOfMutationTournament, startingGenDividen);  //startingGenDividen here is as interval of EDA

        //Set local search by IG operator
        localSearch1 = new localSearchByVNS();//localSearchBy2Opt localSearchBySwap
        GaMain.setLocalSearchOperator(localSearch1, applyLocalSearch, maxNeighborhood);
    }

    public void startMain() {
        openga.util.timeClock timeClock1 = new openga.util.timeClock();
        timeClock1.start();
        GaMain.startGA();
        timeClock1.end();
        //to output the implementation result.
        String implementResult = "";
        int bestInd = getBestSolnIndex(GaMain.getArchieve());
        
        String fileNameArray[] = fileName.split("/");
        fileName = fileNameArray[2] + "\t"+ fileNameArray[3].substring(0, fileNameArray[3].indexOf("."));
        implementResult = fileName + "\t" + this.DEFAULT_crossoverRate+ "\t" + lamda + "\t" + beta + "\t" 
                + numberOfCrossoverTournament + "\t" + numberOfMutationTournament 
                + "\t" + startingGenDividen + "\t"                 
                + GaMain.getArchieve().getSingleChromosome(bestInd).getObjValue()[0] 
                + "\t" + timeClock1.getExecutionTime() / 1000.0  
                +"\n";
        writeFile("SGIO_20170319_Full", implementResult);
        System.out.print(implementResult);
    }

    /**
     * For single objective problem
     * @param arch1
     * @return
     */
    public int getBestSolnIndex(populationI arch1) {
        int index = 0;
        double bestobj = Double.MAX_VALUE;
        for (int k = 0; k < GaMain.getArchieve().getPopulationSize(); k++) {
            if (bestobj > GaMain.getArchieve().getObjectiveValues(k)[0]) {
                bestobj = GaMain.getArchieve().getObjectiveValues(k)[0];
                index = k;
            }
        }
        return index;
    }

    public static void main(String[] args) {
        System.out.println("SGIO_20170319_Full");
        //openga.applications.data.singleMachine singleMachineData = new openga.applications.data.singleMachine();
        int jobSets[] = new int[]{10, 15, 20, 25, 50, 100, 150, 200};//100, 200//10, 15, 20, 25, 50, 100, 150, 200
        int instanceReplication = 15;
        String types[] = new String[]{"low", "med", "high"};//"low", "med", "high"
        int counter = 0;
        int repeatExperiments = 30;
        int totalSolnsToExamine = 125000;

        int popSize[] = new int[]{100};//50, 100, 155, 210 [100]
        double crossoverRate[] = new double[]{0.1},//0.9, 0.5, [0.1]
                mutationRate[] = new double[]{0.5},//0.1, 0.5 {0.5}
                elitism = 0.1;

        //EDA parameters.
        double lamdalearningrate[] = new double[]{0.9}; //0.1, 0.5, [0.9]
        double betalearningrate[] = new double[]{0.1};  //[0.1], 0.5, 0.9
        int numberOfCrossoverTournament[] = new int[]{4};//1, 2, [4]
        int numberOfMutationTournament[] = new int[]{1};//There is no Guided Mutation.
        int startingGenDividen[] = new int[]{2};//[2], 4, 10

        for (int j = 0; j < jobSets.length; j++) {//jobSets.length
            for (int k = 1; k <= instanceReplication; k++) { 
              for(int a = 0 ; a < types.length ; a ++){
                openga.applications.data.singleMachineSetupDynamicData readSingleMachineData1 = new openga.applications.data.singleMachineSetupDynamicData();
                int numberOfJobs = jobSets[j];
                String fileName = "instances/SingleMachineSetupDynamicArrival/"+types[a]+"/"+jobSets[j]+"_"+k+".etp";
                System.out.print(fileName + "\t");
                readSingleMachineData1.setData(fileName, jobSets[j]);
                readSingleMachineData1.getDataFromFile();                                        
                int processingTime[] = readSingleMachineData1.getProcessingTime();
                int setupTime[][] = readSingleMachineData1.getSetupTime();
                int dynamicArrivalTime[] = readSingleMachineData1.getDynamicArrivalTime();   
                
                for (int lx = 0; lx < lamdalearningrate.length; lx++) {
                    for (int bx = 0; bx < betalearningrate.length; bx++) {
                        for (int m = 0; m < numberOfCrossoverTournament.length; m++) {
                            for (int n = 0; n < numberOfMutationTournament.length; n++) {
                                for (int p = 0; p < startingGenDividen.length; p++) {
                                  for(int q = 0 ; q < crossoverRate.length ; q ++){
                                    for (int i = 0; i < repeatExperiments; i++) {
                                        System.out.println("Combinations: " + counter);
                                        singleMachineSetupDynamicArrivalSGIO singleMachine1 = new singleMachineSetupDynamicArrivalSGIO();
                                        singleMachine1.setParameters(popSize[0], crossoverRate[q], mutationRate[0], totalSolnsToExamine);
                                        singleMachine1.setData(numberOfJobs, processingTime, setupTime, dynamicArrivalTime, fileName);
                                        singleMachine1.setEDAinfo(lamdalearningrate[lx], betalearningrate[bx], numberOfCrossoverTournament[m], numberOfMutationTournament[n], startingGenDividen[p]);
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
              }
            }
        }
    }
}
