/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package openga.ObjectiveFunctions;

import openga.applications.flowshopProblem.*;
import openga.ObjectiveFunctions.ObjectiveFunctionFlowShopScheduleI;

/**
 *
 * @author Kuo Yu-Cheng
 */
public interface ObjFunctionPFSSOAWTI extends ObjectiveFunctionFlowShopScheduleI {
  void setOASData(int piTotal , int machineTotal , int[] fristProfit , int[] di , double[] wi , int[][] processingTime);
  void setWriteData(String fileo100x10_0txt);
  void output();
}
