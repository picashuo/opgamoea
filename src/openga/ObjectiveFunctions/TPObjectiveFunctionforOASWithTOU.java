package openga.ObjectiveFunctions;

import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import openga.chromosomes.*;

/**
 *
 * @author Kuo Yu-Cheng
 */
public class TPObjectiveFunctionforOASWithTOU extends TPObjectiveFunctionMTSP implements ObjectiveFunctionOASWithTOUI{
  
  ObjFunctionPFSSOAWTWithTOUTariffs OFPFSSOAWT = new ObjFunctionPFSSOAWTWithTOUTariffs();
  chromosome chromosome1;
  
  //calcCost
  private double on_peakPrice = 0.1327;
  private double mid_peakPrice = 0.0750;
  private double off_peakPrice = 0.0422;
//  private int power = 20;
  private double startTime;
  private double endTime;
  private double midPeakStart;
  private double midPeakEnd;
  private double onPeakEnd;
  private double midPeakEnd2;
  
  private double TOUcostStartTime = 360.0;
  private double TOUcostStartTimeTemp;
  private double TOUcostProcessingTime;
  private double TOUcostCompleteTime;
  private double[] TOUcost;
  private String[] TOUcostStartTimeStr;
  private String[] TOUcostCompleteTimeStr;
  private double totalCost;
  
  private int[] Sequence;

  //Objective Value
  double minimumCost;
  double Revenue;
  double maximumRevenue;
  boolean havePunish;

  //Instance Data
  double[] r;       //  release date.(arrival time)
  double[] p;       //  processing time
  double[] d;       //  due-date
  double[] d_bar;   //  deadline
  double[] e;       //  revenue
  double[] w;       //  weight
  double[] power;       //  power
  double[][] s;     //  setup times
  double[] C;       //completion Time

  public static void main(String[] args) throws IOException {
  
    TPObjectiveFunctionforOASWithTOU TPOAS = new TPObjectiveFunctionforOASWithTOU();
    
    openga.applications.data.OASInstancesWithTOU OAS = new openga.applications.data.OASInstancesWithTOU();
    OAS.setData("@../../instances/SingleMachineOAS/10orders/Tao5/R7/Dataslack_10orders_Tao5R7_10.txt", 10);
    OAS.getDataFromFile();
    
    TPOAS.r = OAS.getR();
    TPOAS.p = OAS.getP();
    TPOAS.d = OAS.getD();
    TPOAS.d_bar = OAS.getD_bar();
    TPOAS.e = OAS.getE();
    TPOAS.w = OAS.getW();
    TPOAS.power = OAS.getPower();
    TPOAS.s = OAS.getS();

    TPOAS.C = new double[TPOAS.r.length];

    chromosome _chromosome1 = new chromosome();
//    chromosome _chromosome2 = new chromosome();
    _chromosome1.setGenotypeAndLength(true, 12, 2);
//    _chromosome2.setGenotypeAndLength(true, 5, 2);
    _chromosome1.generateTwoPartPop(12, 2);
//    _chromosome2.generateTwoPartPop(5, 2);

//    int[] soln = new int[]{0, 2, 1, 3, 0};
//    _chromosome1.setSolution(soln);
//    TPOAS.evaluateAll(_chromosome1, 2);

//    for(int i = 0 ; i < _chromosome1.getLength() ; i++)
//    {
//      System.out.println(_chromosome1.genes[i]);
//    }

    System.out.println("\n" + TPOAS.evaluateAll(_chromosome1, 2));
//    System.out.println("ObjValue1" + _chromosome1.getObjValue()[0]);
//    System.out.println("ObjValue1" + TPOAS.maximumRevenue);
//    TPOAS.maximumRevenue = 0;
//    int[] soln2 = new int[]{0, 4, 3, 8, 7, 2, 9, 1, 6, 5, 7, 3};
//    _chromosome2.setSolution(soln2);
//    TPOAS.evaluateAll(_chromosome2, 2);
//    System.out.println("ObjValue2" + _chromosome2.getObjValue()[0]);
//    System.out.println("ObjValue2" + TPOAS.maximumRevenue);
  }

  @Override
  public void calcObjective() {
//    System.out.print("calcObjective");
    double obj;
    double objectives[];
//    for (int i = 0; i < 1; i++) {

    for (int i = 0; i < population.getPopulationSize(); i++) {
//      System.out.print("population(before): ");
//      for (int j = 0; j < population.getSingleChromosome(i).genes.length; j++) {
//        System.out.print(population.getSingleChromosome(i).genes[j] + " ");
//      }
//      System.out.println("End");
      objectives = population.getObjectiveValues(i);
      obj = evaluateAll(population.getSingleChromosome(i), numberOfSalesmen);
      objectives[indexOfObjective] = obj;
      population.setObjectiveValue(i, objectives);
      chromosome1.setObjValue(objectives);
      population.setSingleChromosome(i, chromosome1);
//      System.out.print("population(after):  ");
//      for (int j = 0; j < population.getSingleChromosome(i).genes.length; j++) {
//        System.out.print(population.getSingleChromosome(i).genes[j] + " ");
//      }
//      System.out.println("End");
    }
  }

  public List<Integer> chromosometoList(chromosome _chromosome1) {
    List<Integer> soln = new ArrayList<Integer>();
    for (int i = 0; i < _chromosome1.genes.length; i++) {
      soln.add(_chromosome1.genes[i]);
    }
    return soln;
  }

  @Override
  public double evaluateAll(chromosome _chromosome1, int numberOfSalesmen) {
    this.setData(_chromosome1, numberOfSalesmen);
    this.calcMaximumRevenue();
    
    return this.getMaximumRevenue() - evaluateAllwithTOUcost(_chromosome1.genes) ;
  }

  public void calcMaximumRevenue() {
    Revenue = 0;
    maximumRevenue = 0;
    int numberOfCities = length - numberOfSalesmen,
            currentPosition = 0, stopPosition = 0,
            index, lastindex = 0;
    double time = 0, dayGap;
    List<Integer> _chromosome1 = new ArrayList<>();
    List<Integer> accept = new ArrayList<>();
    List<Integer> reject = new ArrayList<>();
    List<Integer> salesmen = new ArrayList<>();
    _chromosome1.addAll(chromosometoList(chromosome1));
    for (int i = 0; i < numberOfSalesmen; i++) {
      stopPosition += _chromosome1.get(numberOfCities + i);
      for (int j = currentPosition; j < stopPosition; j++) {
        index = chromosome1.genes[j];
        if (time < r[index]) {
          time = r[index];
        }
        if (j == 0) {
          time += s[0][index];
        }
        if (j != 0) {
          time += s[lastindex+1][index];
        }
        time += p[index];
        C[index] = time;

//        dayGap = Math.max(0, Math.min((d_bar[index] - d[index]), (d_bar[index] - C[index])));
        if (C[index] <= d[index]) {
          Revenue = e[index];
        } else if (C[index] > d[index] && C[index] <= d_bar[index]) {
          Revenue = e[index] - (C[index] - d[index]) * w[index];
        } else {
          Revenue = 0;
        }

        maximumRevenue += Revenue;
//        if(maximumRevenue>90){
//          printResult(time,index,lastindex);
//          System.exit(0);
//        }
        
//              System.out.println("obj: " + (index));
//          System.out.println("release-date: " + r[index]);
//        if (j == 0) {
//          System.out.println("setup-times: " + s[0][index]);
//        }
//        if (j != 0) {
//          System.out.println("setup-times: " + s[lastindex+1][index]);
//        }
//        System.out.println("processing-time: " + p[index]);
//        System.out.println("due-date: " + d[index]);
//        System.out.println("deadline: " + d_bar[index]);
//        System.out.println("completion-time: " + time);
//        System.out.println("weight: " + w[index]);
//        System.out.println("time: " + time);
//        System.out.println("Revenue: " + Revenue);
//        System.out.println("maximumRevenue: " + maximumRevenue);
        

        if (Revenue == 0) {
          
          reject.add(_chromosome1.get(j));
        } else {
          accept.add(_chromosome1.get(j));
        }
        lastindex = index;
      }
      currentPosition += _chromosome1.get(numberOfCities + i);
    }
    salesmen.add(accept.size());
    salesmen.add(reject.size());
    _chromosome1.clear();
    _chromosome1.addAll(accept);
    _chromosome1.addAll(reject);
    _chromosome1.addAll(salesmen);
    chromosome1.setSolution(_chromosome1);
    
//    System.out.print("\nreject : ");
//    for(int i = 0 ; i < reject.size() ; i++)
//    {
//      System.out.print(reject.get(i) + ",");
//    }
//    System.out.print("\naccept : ");
//    for(int i = 0 ; i < accept.size() ; i++)
//    {
//      System.out.print(accept.get(i) + ",");
//    }
    
  }

  public void calcMinimumCost() {
  }
  
  public double evaluateAllwithTOUcost(int[] Sequence) {
    
    totalCost = 0.0;
    
    this.Sequence = Sequence;
    TOUcostStartTimeStr = new String[Sequence.length];
    TOUcostCompleteTimeStr = new String[Sequence.length];
    TOUcost = new double[Sequence.length];
    
    calculateTOUMachineCost();
    
    return totalCost;
  }
  
  public void setData( String startTime, String endTime) throws ParseException {
    SimpleDateFormat simple = new SimpleDateFormat();
    simple.applyPattern("HH:mm");
    this.startTime = (double) simple.parse(startTime).getTime() / (1000 * 60);
    this.endTime = (double) simple.parse(endTime).getTime() / (1000 * 60);
  }
  
  public void calculateTOUMachineCost()
  {
    try {
        int dayTime = 1440 ;
        int acceptedOrders = Sequence[Sequence.length-2] ; // acceptedOrders =  Sequence[8] ������
//        System.out.println("\n" + (Sequence.length-2));
        for (int i = 0; i < acceptedOrders; i++) {
              double tempCompleteTime,tempStartTime;
              tempCompleteTime = 0.0;
              
//              if(accept[i])
//              if(true)
//              {
                TOUcostProcessingTime = C[i] - p[Sequence[i]];
                TOUcostCompleteTime = C[i] + TOUcostStartTime;
                TOUcostStartTimeTemp = (TOUcostStartTime + TOUcostProcessingTime);
                
                if(TOUcostCompleteTime > dayTime && TOUcostStartTimeTemp < dayTime)//�p�G���� > 24 && �}�l < 24
                {
                    tempStartTime = TOUcostStartTimeTemp;
                    tempCompleteTime = TOUcostCompleteTime - dayTime + TOUcostStartTime;

                    TOUcostCompleteTime = dayTime;//�}�l _ ����
                    TOUcostStartTimeStr[i] = ((int)(TOUcostStartTimeTemp / 60) + ":" + ((int)(TOUcostStartTimeTemp % 60)));
                    TOUcostCompleteTimeStr[i] = (((int)(TOUcostCompleteTime / 60)) + ":" + ((int)(TOUcostCompleteTime % 60)));
                    setData(TOUcostStartTimeStr[i] , TOUcostCompleteTimeStr[i]);
                    TOUcost[i] = calculateTOUCost(i);

                    TOUcostStartTimeTemp = TOUcostStartTime;//���� _ ����
                    TOUcostCompleteTime = tempCompleteTime;
                    TOUcostStartTimeStr[i] = ((int)(TOUcostStartTimeTemp / 60) + ":" + ((int)(TOUcostStartTimeTemp % 60)));
                    TOUcostCompleteTimeStr[i] = (((int)(TOUcostCompleteTime / 60)) + ":" + ((int)(TOUcostCompleteTime % 60)));
                    TOUcost[i] += calculateTOUCost(i);
                    
                    TOUcostStartTimeTemp = tempStartTime;
                    TOUcostCompleteTime = tempCompleteTime;
                    TOUcostStartTimeStr[i] = ((int)(TOUcostStartTimeTemp / 60) + ":" + ((int)(TOUcostStartTimeTemp % 60)));
                    TOUcostCompleteTimeStr[i] = (((int)(TOUcostCompleteTime / 60)) + ":" + ((int)(TOUcostCompleteTime % 60)));

                  
                }else if(TOUcostCompleteTime > dayTime && TOUcostStartTimeTemp > dayTime)//�p�G���� > 24 && �}�l > 24
                {
                  while(TOUcostCompleteTime > dayTime && TOUcostStartTimeTemp > dayTime)//�h�Y���24�p�ɤ�
                  {
                    TOUcostCompleteTime = (TOUcostCompleteTime - (dayTime - TOUcostStartTime));
                    TOUcostStartTimeTemp = (TOUcostStartTimeTemp - (dayTime - TOUcostStartTime));
                  }
                  
                  if(TOUcostCompleteTime > dayTime && TOUcostStartTimeTemp < dayTime)//�p�G���� > 24 && �}�l < 24
                  {
                    tempStartTime = TOUcostStartTimeTemp;
                    tempCompleteTime = TOUcostCompleteTime - dayTime + TOUcostStartTime;
                  
                    TOUcostCompleteTime = dayTime;//�}�l _ ����
                    TOUcostStartTimeStr[i] = ((int)(TOUcostStartTimeTemp / 60) + ":" + ((int)(TOUcostStartTimeTemp % 60)));
                    TOUcostCompleteTimeStr[i] = (((int)(TOUcostCompleteTime / 60)) + ":" + ((int)(TOUcostCompleteTime % 60)));
                    setData(TOUcostStartTimeStr[i] , TOUcostCompleteTimeStr[i]);
                    TOUcost[i] = calculateTOUCost(i);

                    TOUcostStartTimeTemp = TOUcostStartTime;//���� _ ����
                    TOUcostCompleteTime = tempCompleteTime;
                    TOUcostStartTimeStr[i] = ((int)(TOUcostStartTimeTemp / 60) + ":" + ((int)(TOUcostStartTimeTemp % 60)));
                    TOUcostCompleteTimeStr[i] = (((int)(TOUcostCompleteTime / 60)) + ":" + ((int)(TOUcostCompleteTime % 60)));
                    TOUcost[i] += calculateTOUCost(i);
                    
                    TOUcostStartTimeTemp = tempStartTime;
                    TOUcostCompleteTime = tempCompleteTime;
                    TOUcostStartTimeStr[i] = ((int)(TOUcostStartTimeTemp / 60) + ":" + ((int)(TOUcostStartTimeTemp % 60)));
                    TOUcostCompleteTimeStr[i] = (((int)(TOUcostCompleteTime / 60)) + ":" + ((int)(TOUcostCompleteTime % 60)));
                    
                  }else
                  {
                    TOUcostStartTimeStr[i] = ((int)(TOUcostStartTimeTemp / 60) + ":" + ((int)(TOUcostStartTimeTemp % 60)));
                    TOUcostCompleteTimeStr[i] = (((int)(TOUcostCompleteTime / 60)) + ":" + ((int)(TOUcostCompleteTime % 60)));

                    setData( TOUcostStartTimeStr[i] , TOUcostCompleteTimeStr[i]);
                    TOUcost[i] = calculateTOUCost(i);
                    
                  }
                  
                }else
                {
                  TOUcostStartTimeStr[i] = ((int)(TOUcostStartTimeTemp / 60) + ":" + ((int)(TOUcostStartTimeTemp % 60)));
                  TOUcostCompleteTimeStr[i] = (((int)(TOUcostCompleteTime / 60)) + ":" + ((int)(TOUcostCompleteTime % 60)));

                  setData( TOUcostStartTimeStr[i] , TOUcostCompleteTimeStr[i]);
                  TOUcost[i] = calculateTOUCost(i);
                  
                }
                
                
                
                
//                System.out.println(" ���� : " + Sequence[i] + " ���x : " + j + " �}�l�ɶ� : " + TOUcostStartTimeStr[i] + " �����ɶ� : " + TOUcostCompleteTimeStr[i] + " �ӶO���� : " + TOUcost[i]);
//              }
//              else
//              {
//                TOUcost[i] = 0;
//                TOUcostStartTimeStr[i] = "0";
//                TOUcostCompleteTimeStr[i] = "0";
////                System.out.println(" ���� : " + Sequence[i] + " ���x : " + j + " �}�l�ɶ� : " + TOUcostStartTimeStr[i] + " �����ɶ� : " + TOUcostCompleteTimeStr[i] + " �ӶO���� : " + TOUcost[i]);
//              }
              
              if(TOUcost[i] >= 0)
              {
              }else
              {
                TOUcost[i] = 0;
              }
              totalCost += TOUcost[i];
//              System.out.println(TOUcost[i]);

        }
        
    } catch (ParseException ex) {
              Logger.getLogger(ObjFunctionPFSSOAWTWithTOUTariffs.class.getName()).log(Level.SEVERE, null, ex);
            }
  }
  
  public double calculateTOUCost(int sequenceCount) throws ParseException {
    SimpleDateFormat simple = new SimpleDateFormat();
    simple.applyPattern("HH:mm");
    double totalPrice = 0.0;
    double totalPeriod = endTime - startTime;
    midPeakStart = (double) simple.parse("07:00").getTime() / (1000 * 60);
    midPeakEnd = (double) simple.parse("15:00").getTime() / (1000 * 60);
    onPeakEnd = (double) simple.parse("20:00").getTime() / (1000 * 60);
    midPeakEnd2 = (double) simple.parse("22:00").getTime() / (1000 * 60);
    
    if (totalPeriod < 0) {
      System.out.println("endTime can not be greater than startTime!!");
      return -1;
    } else {
      if (startTime < midPeakStart)//Off_peak
      {
        if (endTime > midPeakStart) {
          totalPrice += ((((midPeakStart - startTime) / 60) * off_peakPrice) * (((midPeakStart - startTime)  / totalPeriod) * power[sequenceCount]));
          startTime = midPeakStart;
        } else {
          totalPrice += ((((endTime - startTime) / 60) * off_peakPrice) * (((endTime - startTime) / totalPeriod) * power[sequenceCount]));
        }
      }
      if (startTime >= midPeakStart && startTime < midPeakEnd)//Mid-peak 07:00 ~ 15:00
      {
        if (endTime > midPeakEnd) {
          totalPrice += ((((midPeakEnd - startTime) / 60) * mid_peakPrice) * (((midPeakEnd - startTime) / totalPeriod) * power[sequenceCount]));
          startTime = midPeakEnd;
        } else {
          totalPrice += ((((endTime - startTime) / 60) * mid_peakPrice) * (((endTime - startTime) / totalPeriod) * power[sequenceCount]));
        }
      }
      if (startTime >= midPeakEnd && startTime < onPeakEnd)//On-peak 15:00 ~ 20:00
      {
        if (endTime > onPeakEnd) {
          totalPrice += ((((onPeakEnd - startTime) / 60) * on_peakPrice) * (((onPeakEnd - startTime) / totalPeriod) * power[sequenceCount]));
          startTime = onPeakEnd;
        } else {
          totalPrice += ((((endTime - startTime) / 60) * on_peakPrice) * (((endTime - startTime) / totalPeriod) * power[sequenceCount]));
        }
      }
      if (startTime >= onPeakEnd && startTime < midPeakEnd2)//Mid-peak 20:00 ~ 22:00
      {
        if (endTime > midPeakEnd2) {
          totalPrice += ((((midPeakEnd2 - startTime) / 60) * mid_peakPrice) * (((midPeakEnd2 - startTime) / totalPeriod) * power[sequenceCount]));
          startTime = midPeakEnd2;
        } else {
          totalPrice += ((((endTime - startTime) / 60) * mid_peakPrice) * (((endTime - startTime) / totalPeriod) * power[sequenceCount]));
        }
      }
      if (startTime >= midPeakEnd2)//Off_peak
      {
        totalPrice += ((((endTime - startTime) / 60) * off_peakPrice) * (((endTime - startTime) / totalPeriod) * power[sequenceCount]));
      }
      return totalPrice;
    }
  }
  
  public void printResult(double time, int index, int lastindex) {
    for(int i = 0; i < chromosome1.genes.length; i++){
      System.out.println("obj: " + (index));
        if (time < r[index]) {
          System.out.println("release-date: " + r[index]);
        }
        if (i == 0) {
          System.out.println("setup-times: " + s[0][index]);
        }
        if (i != 0) {
          System.out.println("setup-times: " + s[index + 1][lastindex]);
        }
        System.out.println("processing-time: " + p[index]);
        System.out.println("due-date: " + d[index]);
        System.out.println("deadline: " + d_bar[index]);
        System.out.println("completion-time: " + time);
        System.out.println("weight: " + w[index]);
        System.out.println("time: " + time);
        System.out.println("Revenue: " + Revenue);
        System.out.println("maximumRevenue: " + maximumRevenue);
    }
    
  }

  public void setOASData(double[] r, double[] p, double[] d, double[] d_bar, double[] e, double[] w, double[] power, double[][] s, int numberOfSalesmen) {
    this.r = r;
    this.p = p;
    this.d = d;
    this.d_bar = d_bar;
    this.e = e;
    this.w = w;
    this.power = power;
    this.s = s;
    this.numberOfSalesmen = numberOfSalesmen;
    C = new double[p.length];
  }

  @Override
  public void setData(chromosome chromosome1, int numberOfSalesmen) {
    this.chromosome1 = chromosome1;
    this.numberOfSalesmen = numberOfSalesmen;
    length = chromosome1.getLength();
  }

  public double getMinimumCost() {
    return minimumCost;
  }

  public double getMaximumRevenue() {
    return maximumRevenue;
  }

  @Override
  public double[] getObjectiveValues(int index) {
    double objectives[];
    objectives = chromosome1.getObjValue();
    double obj = evaluateAll(chromosome1, numberOfSalesmen);
    objectives[0] = obj;
    chromosome1.setObjValue(objectives);
    return chromosome1.getObjValue();
  }

  public void setPowerData(double[] power) {
    throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
  }

  @Override
  public void setPowerData(double[] p, double[] power) {
    throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
  }
  
}
