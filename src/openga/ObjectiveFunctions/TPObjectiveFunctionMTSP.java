package openga.ObjectiveFunctions;

import openga.chromosomes.*;

/**
 *
 * @author YU-TANG CHANG
 */
public class TPObjectiveFunctionMTSP extends ObjectiveFunctionMTSP {

  int numberOfSalesmen;

  public void calcObjective() {
    double obj;
    double objectives[];

    for (int i = 0; i < population.getPopulationSize(); i++) {
      objectives = population.getObjectiveValues(i);
      obj = evaluateAll(population.getSingleChromosome(i), numberOfSalesmen);
//      System.out.println(obj);
      objectives[indexOfObjective] = obj;
      population.setObjectiveValue(i, objectives);
//      for(int j=0; j<population.getObjectiveValues(i).length; j++){
//        System.out.println(population.getObjectiveValues(i)[j]);
//      }
    }
  }

  public double evaluateAll(chromosome _chromosome1, int numberOfSalesmen) {
    TPforMTSPDistanceCalculation TPforDistanceCalculation1 = new TPforMTSPDistanceCalculation();

    if (typeName.equals("TotalDistance")) {
      TPforDistanceCalculation1 = new TPforMTSPDistanceCalculation();
    } else {//MaxDistance
      TPforDistanceCalculation1 = new TPforMTSPMaxDistanceCalculation();
    }

    TPforDistanceCalculation1.setData(distanceToOriginal, distanceMatrix, _chromosome1, numberOfSalesmen);
    TPforDistanceCalculation1.calcObjective();
//    System.out.println(TPforDistanceCalculation1.getObjectiveValue());
    return TPforDistanceCalculation1.getObjectiveValue();
  }

  @Override
  public void setTSPData(double originalPoint[], double coordinates[][], int numberOfSalesmen) {
    this.numberOfSalesmen = numberOfSalesmen;
    this.originalPoint = originalPoint;
    this.coordinates = coordinates;
    openga.ObjectiveFunctions.util.generateMatrix_EuclideanDist generateMatrix_EuclideanDist1
            = new openga.ObjectiveFunctions.util.generateMatrix_EuclideanDist();
    generateMatrix_EuclideanDist1.setData(coordinates);
    generateMatrix_EuclideanDist1.setDistanceMatrixData();
    distanceMatrix = generateMatrix_EuclideanDist1.getMatrix();

    openga.ObjectiveFunctions.util.generateMatrix_EuclideanDist2 generateMatrix_EuclideanDist21
            = new openga.ObjectiveFunctions.util.generateMatrix_EuclideanDist2();
    generateMatrix_EuclideanDist21.setData(originalPoint, coordinates);
    generateMatrix_EuclideanDist21.setDistanceMatrixData();
    distanceToOriginal = generateMatrix_EuclideanDist21.getMatrix2();
  }
}
