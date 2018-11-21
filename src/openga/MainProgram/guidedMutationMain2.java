package openga.MainProgram;
import openga.chromosomes.*;
import openga.operator.selection.*;
import openga.operator.crossover.*;
import openga.operator.mutation.*;
import openga.operator.clone.*;
import openga.ObjectiveFunctions.*;
import openga.Fitness.*;
import openga.operator.miningGene.*;

/**
 * <p>Title: The OpenGA project which is to build general framework of Genetic algorithm.</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2008</p>
 * <p>Company: Yuan-Ze University</p>
 * @author Chen, Shih-Hsin
 * @version 1.0
 */

public class guidedMutationMain2 extends singleThreadGAwithProbabilityMatrix {
    public guidedMutationMain2() {
    }
    PBIL PBIL1;
    int strategy = 1;
    double container[][];
    double lamda = 0.1; //learning rate
    double beta = 0.9;
    GuidedMuation GuidedMuation1;
    populationI artificialPopulation = new population();
    String generationResults="";
    double genoDiversityValue=0;
    
    public void setGuidedMutationInfo(double lamda, double beta){
        this.lamda = lamda;
        this.beta = beta;
    }
    
    public void startGA(){
        Population = initialStage();   //test
        
        intialArttificalPopulation();
        //evaluate the objective values and calculate fitness values
        ProcessObjectiveAndFitness();
        //intialOffspringPopulation();
        archieve = findParetoFront(Population, 0);
        PBIL1 = new PBIL(Population, lamda);
        container = PBIL1.getContainer();
        
        for(int i = 0 ; i < generations ; i ++ ){
            //System.out.println("generations "+i);
            // output diversity
            
            diversityMeasure();
            //generationResults = "";
            //generationResults = i+"\t"+getBest()+"\t"+genoDiversityValue+"\n";
            //System.out.print(generationResults);
            //writeFile("EAGSGA_obj_Diversity_952", generationResults);
            
            currentGeneration = i;
            
/*    generation_result
      generationResults = i+"\t"+getBest()+"\n";
      System.out.print(generationResults);
      writeFile("EAG_Generation_XXX_125000" , generationResults);
 */
            if (i % 6 == 0){
                Population = selectionStage(Population);
                
                //collect gene information, it's for mutation matrix
                PBIL1.setData(Population);
                PBIL1.startStatistics();
                container = PBIL1.getContainer();
                ProbabilityMatrix();
                
            } else{
                Population = selectionStage(Population);
                if(genoDiversityValue>0.3){
                    Population = crossoverStage(Population);
                    //System.out.println("crossover1, " + genoDiversityValue );
                }else{
                    Population = crossover2Stage(Population);
                    //System.out.println("crossover2, " + genoDiversityValue );
                }
                Population = mutationStage(Population);
                
                //clone
                if(applyClone == true){
                    Population = cloneStage(Population);
                }
                
                //evaluate the objective values and calculate fitness values
                ProcessObjectiveAndFitness();
            }
            
            populationI tempFront = (population)findParetoFront(Population,0);
            archieve = updateParetoSet(archieve,tempFront);
            
            
            
            
/*
      if(i == 500){
        String generationResults = "";
        generationResults = i+"\t"+beta+"\t"+getBest()+"\n";
        System.out.print(generationResults);
        writeFile("EAG_0223_" + i, generationResults);
      }
      if(i == 750){
        String generationResults = "";
        generationResults = i+"\t"+beta+"\t"+getBest()+"\n";
        System.out.print(generationResults);
        writeFile("EAG_0223_" + i, generationResults);
      }
      if(i == 1000){
        String generationResults = "";
        generationResults = i+"\t"+beta+"\t"+getBest()+"\n";
        System.out.print(generationResults);
        writeFile("EAG_0223_" + i, generationResults);
      }
 */
            
        }
    }
    
    public void ProbabilityMatrix(){
        Population = selectionStage(Population);
        openga.operator.miningGene.GuidedMuation probMatrix1 = new openga.operator.miningGene.GuidedMuation();
        probMatrix1.setData(Population, artificialPopulation);//archieve, Population, artificialPopulation
        probMatrix1.setDuplicateRate(1);
        probMatrix1.setStrategy(strategy);
        probMatrix1.setEDAinfo(container, archieve, beta);
        probMatrix1.startStatistics();
        artificialPopulation = cloneStage(artificialPopulation);
        evalulatePop(artificialPopulation);
        Population = replacementStage(Population, artificialPopulation);
    }
    
    public void intialArttificalPopulation(){
        artificialPopulation.setGenotypeSizeAndLength(encodeType, fixPopSize, length, numberOfObjs);
        artificialPopulation.initNewPop();
    }
    
    public populationI crossover2Stage(populationI Population){
        Crossover2.setData(crossoverRate, Population);
        Crossover2.startCrossover();
        Population = Crossover2.getCrossoverResult();
        
        return Population;
    }
    
    public void diversityMeasure(){
        
        openga.adaptive.diversityMeasure.diversityMeasureI diversityMeasure1 = new openga.adaptive.diversityMeasure.HammingDistance();  //GenoTypePositionDifference() HammingDistance()
        diversityMeasure1.setData(Population);
        //diversityMeasure1.CalcDiversity();
        diversityMeasure1.startCalcDiversity();
        genoDiversityValue = diversityMeasure1.getPopulationDiversityValue();
        
        
    }
    
    
    
    
    
    
}