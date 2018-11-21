package openga.adaptive.diversityMeasure;
import openga.chromosomes.*;
/**
 * <p>Title: The OpenGA project which is to build general framework of Genetic algorithm.</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2005</p>
 * <p>Company: Yuan-Ze University</p>
 * @author Chen, Shih-Hsin
 * @version 1.0
 */

public interface diversityMeasureI {
  public void setData(populationI population1);
  public void startCalcDiversity();
   public double getPopulationDiversityValue();
}