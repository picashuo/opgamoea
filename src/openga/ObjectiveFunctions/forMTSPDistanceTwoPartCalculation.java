/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package openga.ObjectiveFunctions;
import openga.chromosomes.*;
/**
 *
 * @author UlyssesZhang
 */
public class forMTSPDistanceTwoPartCalculation extends forMTSPDistanceCalculation{
    
      int numberOfSalesmen;
  
    public void setData(double distanceToOriginal[], double distanceMatrix[][], chromosome chromosome1, int numberOfSalesmen){
      this.distanceToOriginal = distanceToOriginal;
      this.distanceMatrix = distanceMatrix;
      this.chromosome1 = chromosome1;
      this.numberOfSalesmen = numberOfSalesmen;
      length = distanceMatrix.length;
    }
  
    public void calcObjective(){
      //to get the distance of each salesmen
      int currentPosition = 0;//To record the position of the Part I chromosome
      
      for(int k = 0 ; k < numberOfSalesmen ; k ++){
        //from the original point to the first position.
        objVal += distanceToOriginal[chromosome1.genes[currentPosition]];
        
        int numberOfCities = length - numberOfSalesmen;
        int stopPosition = numberOfCities + currentPosition - 1;

        for(int i = currentPosition ; i <= stopPosition ; i ++ ){
          if(i < stopPosition){
            int index1 = chromosome1.genes[i];
            int index2 = chromosome1.genes[i+1];
            objVal += distanceMatrix[index1][index2];       
          }
          else{
            //The last point then go back to original point.
            objVal += distanceToOriginal[chromosome1.genes[currentPosition]];
          }
          currentPosition ++;
        }
      }
    }
}
