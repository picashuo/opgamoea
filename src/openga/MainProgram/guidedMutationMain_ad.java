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

public class guidedMutationMain_ad extends singleThreadGAwithProbabilityMatrix {
  public guidedMutationMain_ad() {
  }
  PBIL PBIL1;
  int strategy = 1;
  double container[][];
  double lamda; //learning rate
  double beta;
  GuidedMuation GuidedMuation1;
  populationI artificialPopulation = new population();

  boolean executeAdaptive;

  public void setGuidedMutationInfo(double lamda, double beta){
    this.lamda = lamda;
    this.beta = beta;
  }

  public void setApaptive(boolean executeAdaptive){
      this.executeAdaptive = executeAdaptive;
  }
  
  
  public void startGA(){
    Population = initialStage();
    intialArttificalPopulation();
    //evaluate the objective values and calculate fitness values
    ProcessObjectiveAndFitness();
    //intialOffspringPopulation();
    archieve = findParetoFront(Population, 0);
    PBIL1 = new PBIL(Population, lamda);
    container = PBIL1.getContainer();
    double temp = beta;
    for(int i = 0 ; i < generations ; i ++ ){
      
      currentGeneration = i;
      Population = selectionStage(Population);


      if(executeAdaptive){
          beta = temp * (double)(generations - currentGeneration)/(double)generations;
          System.out.println(beta);
      }
      
      
      //collect gene information, it's for mutation matrix
      PBIL1.setData(Population);
      PBIL1.startStatistics();
      container = PBIL1.getContainer();
      ProbabilityMatrix();

      populationI tempFront = (population)findParetoFront(Population, 0);
      archieve = updateParetoSet(archieve, tempFront);
      
      
   /*   if(i % 10 == 0){
        String generationResults = "";
        generationResults = "EAG_1250_"+ length+"\t"+i+"\t"+beta+"\t"+getBest()+"\n";
        System.out.print(generationResults);
        writeFile("EAG_beta_9_0_125000_" + length, generationResults);    
      }   */
      
      
    }
  }

  public void ProbabilityMatrix(){
      Population = selectionStage(Population);
      openga.operator.miningGene.GuidedMuation probMatrix1 = new openga.operator.miningGene.GuidedMuation();
      probMatrix1.setData(Population, artificialPopulation);//archieve, Population, artificialPopulation
      probMatrix1.setDuplicateRate(1);
      probMatrix1.setStrategy(strategy);
      probMatrix1.setEDAinfo(container, archieve, beta);
      //System.out.println(beta);
      probMatrix1.startStatistics();
      artificialPopulation = cloneStage(artificialPopulation);
      evalulatePop(artificialPopulation);
      Population = replacementStage(Population, artificialPopulation);
  }

  public void intialArttificalPopulation(){
    artificialPopulation.setGenotypeSizeAndLength(encodeType, fixPopSize, length, numberOfObjs);
    artificialPopulation.initNewPop();
  }

}