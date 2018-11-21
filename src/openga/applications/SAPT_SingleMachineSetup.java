package openga.applications;
import openga.chromosomes.chromosome;
import openga.ObjectiveFunctions.*;
//import openga.Selection.*;
import openga.MainProgram.*;
import openga.util.fileWrite1;

import openga.applications.data.singleMachineSetupData;
/**
 * <p>Title: </p>
 * <p>Description: The SAPT heuristic is to construct initial solutions for
 * single machine schedule problem with setup considerations.</p>
 * <p>Copyright: Copyright (c) 2006</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class SAPT_SingleMachineSetup{
    public SAPT_SingleMachineSetup() {
    }
    
    int numberOfJobs = 8;
    double processingTime[][];
    double backupProcessingTime[][];
    public int B[];
    int sequence[];
    int bestsequence[];
    double obj;//earliness + taridness
    //r is the middle job.
    int r = 4;
    double previousMinAP = 0;
    
    double miniAP = 0;
    int i_index=0;
    int j_index=0;
    int current=1;
    
    int formerIndex;
    double objFormer1 = 0;
    int laterIndex;
    double objLater1 = 0;
    int gpitimes = 0;
    
    
    public void setData(int numberOfJobs, double processingTime[][]){
        this.numberOfJobs = numberOfJobs;
        this.processingTime = processingTime;
        sequence = new int[numberOfJobs];
        B = new int[numberOfJobs];
        current=1;
        i_index = 0;
        j_index = 0;
        bestsequence = new int[numberOfJobs];
        gpitimes=0;
    }
    
    public int[] getSmallestIndex2(double matrix[][]){
        int index[] = new int[2];
        boolean srh=true;
        for(int i = i_index ; i < matrix.length ; i ++ ){
            for(int j = j_index ; j < matrix[0].length ; j ++ ){
                if(i != j && matrix[i][j] != -1 && matrix[i][j] == miniAP){
                    index[0] = i;
                    index[1] = j;
                    srh = false;
                    break;
                }
            }
            if(srh == false){
                break;
            }else{
                j_index = 0;
            }
        }
//                System.out.println(index[0]+","+index[1]+" " +miniAP );
//        System.exit(0);
        
        previousMinAP = miniAP;
        i_index = index[0];
        j_index = index[1]+1 ;
        return index;
        
    }
    
    public int[] getSmallestIndex(double matrix[][]){
        int index[] = new int[2];
        double minVal = Integer.MAX_VALUE;
        for(int i = 0 ; i < matrix.length ; i ++ ){
            for(int j = 0 ; j < matrix[0].length ; j ++ ){
                if(i != j && matrix[i][j] != -1 && matrix[i][j] < minVal && matrix[i][j] >= previousMinAP){
                    minVal = matrix[i][j];
                    index[0] = i;
                    index[1] = j;
                }
            }
        }
        previousMinAP = minVal;
        return index;
    }
    
    public int getsmallesTimes(double matrix[][]){
        int times=0;
        double minVal = Double.MAX_VALUE;
        for(int i = 0 ; i < matrix.length ; i ++ ){
            for(int j = 0 ; j < matrix[0].length ; j ++ ){
                if(i != j && matrix[i][j] < minVal && matrix[i][j] > miniAP){
                    minVal = matrix[i][j];
                }
//                System.out.print(matrix[i][j] + " ");
            }
//            System.out.println();
        }
        
        miniAP = minVal;
        
//        System.out.println(minVal);
//        System.out.println("small index "+i_index +","+j_index);
        for(int i = 0 ; i < matrix.length ; i ++ ){   //find the same times
            for(int j = 0 ; j < matrix[0].length ; j ++ ){
                if(i != j && (int)matrix[i][j] == minVal){
                    times++;
                }
            }
        }
        
//        System.out.println("small index times"+times);
        return times;
    }
    
    public void setindex(){
        i_index = 0;
        j_index = 0;
    }
    
    public void SAPT(){
        for(int i = 0 ; i < numberOfJobs ; i ++ ){
            sequence[i] = i;
        }
        for(int i = 0 ; i < numberOfJobs ; i ++ ){
            B[i]= 0 ;
        }
        int index[] = new int[2];
        findMiddlePosition(sequence);

//        if (current==1){
//            index[0] = i_index;
//            index[1]= j_index;
//            current++;
//        }else{
        index = getSmallestIndex2(processingTime);
        
//            index[0] = i_index;
//            index[1]= j_index;
        
//        }
//        System.out.println(index[0]+ " "+index[1]);
        //   int index[] = getSmallestIndex(processingTime);
        B[r+1] = index[1];
        B[r] = index[0];
        
        dropRow(index[0], processingTime);
        dropColumn(index[1], processingTime);
        processingTime[index[0]][index[1]] = -1;
        dropJobs(index[0]);
        dropJobs(index[1]);
        
//        printSequence(sequence);
//        printSequence(B);
        
        
        int sizeB = 1, sizeA = 1;
        while(sizeB <= r || sizeA < (numberOfJobs - r-1)){
//            index = getSmallestIndex(processingTime);
//            int jobIndex = index[0];
//            int jobIndex2 = index[1];
            
            int jobJ = B[r - sizeB + 1];
            int jobI = B[r + sizeA];
            int jobIndex = getRandomJobsi(jobJ);
            int jobIndex2 = getRandomJobsj(jobI);
            
            int penaltyB = r - sizeB;
            int penaltyA = numberOfJobs - (r + sizeA);
            //penaltyB = 1;
            //penaltyA = 1;
            
            if(r - sizeB < 0){
                penaltyB = Integer.MAX_VALUE;
            }
            
            if(r + sizeA + 1 == numberOfJobs){
                penaltyA = Integer.MAX_VALUE;
            }
            
            //start to test //int jobFix, int job1, int job2, int penalty
//            getFormerPosValue(jobI, jobIndex, jobIndex2, penaltyB);
//            getLaterPosValue(jobJ, jobIndex, jobIndex2, penaltyA);
            objFormer1 = processingTime[jobIndex][jobJ]*penaltyB;
            objLater1 = processingTime[jobI][jobIndex2]*penaltyA;
//            System.out.println(jobIndex+","+jobIndex2+"  "+objFormer1+","+objLater1);
            
            if(objFormer1 <= objLater1){//to put the job ahead at position r - sizeB
                B[r - sizeB] = jobIndex;
                dropRow(jobIndex, processingTime);
                dropColumn(jobJ, processingTime);
                sizeB += 1;
                dropJobs(jobIndex);
            } else{//to put the job at [r + rizeA + 1]
                B[r + sizeA + 1] = jobIndex2;
                dropRow(jobI, processingTime);
                dropColumn(jobIndex2, processingTime);
                sizeA += 1;
                dropJobs(jobIndex2);
                
            }
      /*
      if(processingTime[jobIndex][jobI]*(penaltyB) <= processingTime[jobJ][jobIndex]*(penaltyA)){
        B[r - sizeB] = jobIndex;
        dropRow(jobIndex, processingTime);
        dropColumn(jobI, processingTime);
        sizeB += 1;
      }
      else{
        B[r + sizeA + 1] = jobIndex;
        dropRow(jobJ, processingTime);
        dropColumn(jobIndex, processingTime);
        sizeA += 1;
      }
      dropJobs(jobIndex);
       */
            
            //System.out.println(sizeB+" "+sizeA);
//            printSequence(B);
        }//end while
//        System.out.println("result");
//        printSequence(B);
        calcObj(B);
//        System.out.println("obj "+obj);
//        System.exit(0);
    }
    
    public void SAPT3(){
        for(int i = 0 ; i < numberOfJobs ; i ++ ){
            sequence[i] = i;
        }
        int index[] = new int[2];
        findMiddlePosition(sequence);
        
//        if (current==1){
//            index[0] = i_index;
//            index[1]= j_index;
//            current++;
//        }else{
        index = getSmallestIndex2(processingTime);
//            index[0] = i_index;
//            index[1]= j_index;
        
//        }
//        System.out.println(index[0]+ " "+index[1]);
        //   int index[] = getSmallestIndex(processingTime);
        B[r+1] = index[1];
        B[r] = index[0];
        
        dropRow(index[0], processingTime);
        dropColumn(index[1], processingTime);
        dropJobs(index[0]);
        dropJobs(index[1]);
        
//        printSequence(sequence);
//        printSequence(B);
        
        
        int sizeB = 1, sizeA = 1;
        while(sizeB <= r || sizeA < (numberOfJobs - r-1)){
//            index = getSmallestIndex(processingTime);
//            int jobIndex = index[0];
//            int jobIndex2 = index[1];
            
            int jobJ = B[r - sizeB + 1];
            int jobI = B[r + sizeA];
            int jobIndex = getRandomJobsi(jobJ);
            int jobIndex2 = getRandomJobsj(jobI);
            
            int penaltyB = r - sizeB;
            int penaltyA = numberOfJobs - (r + sizeA);
            //penaltyB = 1;
            //penaltyA = 1;
            
            if(r - sizeB < 0){
                penaltyB = Integer.MAX_VALUE;
            }
            
            if(r + sizeA + 1 == numberOfJobs){
                penaltyA = Integer.MAX_VALUE;
            }
            
            //start to test //int jobFix, int job1, int job2, int penalty
//            getFormerPosValue(jobI, jobIndex, jobIndex2, penaltyB);
//            getLaterPosValue(jobJ, jobIndex, jobIndex2, penaltyA);
            objFormer1 = processingTime[jobIndex][jobJ]*penaltyB;
            objLater1 = processingTime[jobI][jobIndex2]*penaltyA;
//            System.out.println(jobIndex+","+jobIndex2+"  "+objFormer1+","+objLater1);
            
            
            if(objFormer1 <= objLater1){//to put the job ahead at position r - sizeB
                B[r - sizeB] = jobIndex;
                dropRow(jobIndex, processingTime);
                dropColumn(jobJ, processingTime);
                sizeB += 1;
                dropJobs(jobIndex);
            } else{//to put the job at [r + rizeA + 1]
                B[r + sizeA + 1] = jobIndex2;
                dropRow(jobI, processingTime);
                dropColumn(jobIndex2, processingTime);
                sizeA += 1;
                dropJobs(jobIndex2);
                
            }
      /*
      if(processingTime[jobIndex][jobI]*(penaltyB) <= processingTime[jobJ][jobIndex]*(penaltyA)){
        B[r - sizeB] = jobIndex;
        dropRow(jobIndex, processingTime);
        dropColumn(jobI, processingTime);
        sizeB += 1;
      }
      else{
        B[r + sizeA + 1] = jobIndex;
        dropRow(jobJ, processingTime);
        dropColumn(jobIndex, processingTime);
        sizeA += 1;
      }
      dropJobs(jobIndex);
       */
            
            //System.out.println(sizeB+" "+sizeA);
//            printSequence(B);
        }//end while
//        System.out.println("result");
//        printSequence(B);
        calcObj(B);
//        System.out.println("obj "+obj);
//        System.exit(0);
    }
    
    public int getRandomJobsi(int js){
        double tmpap = Double.MAX_VALUE;
        int itmp = 0;
        for (int jobIndex=0;jobIndex<numberOfJobs;jobIndex++){
            if(jobIndex!=js && processingTime[jobIndex][js]!=-1.0 && sequence[jobIndex]!=-1){
                if(tmpap > (processingTime[jobIndex][js])){//to put the job ahead at position r - sizeB
                    tmpap = processingTime[jobIndex][js];
                    itmp = jobIndex;
                }
            }
        }
        return itmp;
    }
    
    
    public int getRandomJobsj(int is){
        double tmpap = Double.MAX_VALUE;
        int jtmp = 0;
        for (int jobIndex=0;jobIndex<numberOfJobs;jobIndex++){
            if(jobIndex!=is && processingTime[is][jobIndex]!=-1.0 && sequence[jobIndex]!=-1){
                if(tmpap > (processingTime[is][jobIndex])){//to put the job ahead at position r - sizeB
                    tmpap = processingTime[is][jobIndex];
                    jtmp = jobIndex;
                }
            }
        }
        return jtmp;
    }
    
    
    public void SAPT2(){
        for(int i = 0 ; i < numberOfJobs ; i ++ ){
            sequence[i] = i;
        }
        int index[] = new int[2];
        findMiddlePosition(sequence);
        
//        if (current==1){
//            index[0] = i_index;
//            index[1]= j_index;
//            current++;
//        }else{
        index = getSmallestIndex2(processingTime);
//            index[0] = i_index;
//            index[1]= j_index;
        
//        }
//        System.out.println(index[0]+ " "+index[1]);
        //   int index[] = getSmallestIndex(processingTime);
        B[r+1] = index[1];
        B[r] = index[0];
        
        dropRow(index[0], processingTime);
        dropColumn(index[1], processingTime);
        dropJobs(index[0]);
        dropJobs(index[1]);
        
//        printSequence(sequence);
//        printSequence(B);
        
        
        int sizeB = 1, sizeA = 1;
        while(sizeB <= r || sizeA < (numberOfJobs - r-1)){
//            index = getSmallestIndex(processingTime);
//            int jobIndex = index[0];
//            int jobIndex2 = index[1];
            int jobIndex = getRandomJobs();
            int jobIndex2 = getRandomJobs();
            int jobI = B[r - sizeB + 1];
            int jobJ = B[r + sizeA];
            
            int penaltyB = r - sizeB;
            int penaltyA = numberOfJobs - (r + sizeA);
            //penaltyB = 1;
            //penaltyA = 1;
            
            if(r - sizeB < 0){
                penaltyB = Integer.MAX_VALUE;
            }
            
            if(r + sizeA + 1 == numberOfJobs){
                penaltyA = Integer.MAX_VALUE;
            }
            
            //start to test //int jobFix, int job1, int job2, int penalty
            getFormerPosValue(jobI, jobIndex, jobIndex2, penaltyB);
            getLaterPosValue(jobJ, jobIndex, jobIndex2, penaltyA);
            
            if(objFormer1 <= objLater1){//to put the job ahead at position r - sizeB
                B[r - sizeB] = formerIndex;
                dropRow(formerIndex, processingTime);
                dropColumn(jobI, processingTime);
                sizeB += 1;
                dropJobs(formerIndex);
            } else{//to put the job at [r + rizeA + 1]
                B[r + sizeA + 1] = laterIndex;
                dropRow(jobJ, processingTime);
                dropColumn(laterIndex, processingTime);
                sizeA += 1;
                dropJobs(laterIndex);
            }
      /*
      if(processingTime[jobIndex][jobI]*(penaltyB) <= processingTime[jobJ][jobIndex]*(penaltyA)){
        B[r - sizeB] = jobIndex;
        dropRow(jobIndex, processingTime);
        dropColumn(jobI, processingTime);
        sizeB += 1;
      }
      else{
        B[r + sizeA + 1] = jobIndex;
        dropRow(jobJ, processingTime);
        dropColumn(jobIndex, processingTime);
        sizeA += 1;
      }
      dropJobs(jobIndex);
       */
            
            //System.out.println(sizeB+" "+sizeA);
            //printSequence(B);
        }//end while
//        printSequence(B);
        calcObj(B);
        //System.out.println("obj "+obj);
    }
    
    public int[] GPI(int _sequence[]){
        double obj_pre = 0 , obj_after = 0;
        for (int i =0;i<numberOfJobs-1;i++){
            for(int j=i+1;j<numberOfJobs;j++){
                obj_pre = calcObj2(_sequence);
                _sequence = swapJobs2(_sequence,i,j);
                obj_after = calcObj2(_sequence);
                gpitimes++;
                if (obj_after >= obj_pre){
                    _sequence = swapJobs2(_sequence,i,j);
                }else{
                    i=0;
                    j=1;
                }
            }
        }
//        printSequence(_sequence);
        return _sequence;
    }
    
    public double calcObj2(int _sequence[]){
        double obj2 = 0;
        
        //to calculate the Total earliness
        for(int i = 0 ; i < r ; i ++ ){
            obj2 += (i+1)*processingTime[_sequence[i]][_sequence[i+1]];
        }
        
        //to calculate the Total tardiness
        for(int i = r ; i < numberOfJobs - 1 ; i ++ ){
            obj2 += (numberOfJobs - i - 1)*processingTime[_sequence[i]][_sequence[i+1]];
        }
        return obj2;
    }
    
    public int[] swapJobs2(int _sequence[], int pos1, int pos2){
        int temp = _sequence[pos1];
        _sequence[pos1] = _sequence[pos2];
        _sequence[pos2] = temp;
        return _sequence;
    }
    
    public double findMiddlePosition(int _seq[]){
        if(_seq.length % 2 == 0){
            r = _seq.length / 2;
        } else{
            r = (_seq.length+1) / 2;
        }
        
        if(r > 0){
            r -= 1;
        }
        return r;
    }
    
    
    public int getRandomJobs(){
        int jobIndex = (int)(Math.random()*numberOfJobs);
        while(sequence[jobIndex] == -1){
            jobIndex = (jobIndex + 1) % numberOfJobs;
        }
        return jobIndex;
    }
    
    public void getFormerPosValue(int jobFix, int job1, int job2, int penalty){
        if(processingTime[job1][jobFix] <= processingTime[job2][jobFix]){
            formerIndex = job1;
            objFormer1 = processingTime[job1][jobFix]*penalty;
        } else{
            formerIndex = job2;
            objFormer1 = processingTime[job2][jobFix]*penalty;
        }
    }
    
    
    public void getLaterPosValue(int jobFix, int job1, int job2, int penalty){
        if(processingTime[job1][jobFix] <= processingTime[job2][jobFix]){
            laterIndex = job1;
            objLater1 = processingTime[jobFix][job1]*penalty;
        } else{
            laterIndex = job2;
            objLater1 = processingTime[jobFix][job2]*penalty;
        }
    }
    
    
    public void dropJobs(int jobIndex){
        sequence[jobIndex] = -1;
    }
    
    public double[][] dropColumn(int index, double matrix[][]){
        for(int i = 0 ; i < matrix.length ; i ++ ){
            matrix[i][index] = -1;
        }
        return matrix;
    }
    
    public double[][] dropRow(int index, double matrix[][]){
        for(int i = 0 ; i < matrix.length ; i ++ ){
            matrix[index][i] = -1;
        }
        return matrix;
    }
    
    public void swapJobs(int pos1, int pos2){
        int temp = B[pos1];
        B[pos1] = B[pos2];
        B[pos2] = temp;
    }
    
    public void calcObj(int _sequence[]){
        obj = 0;
        restoreAP();
        
        //to calculate the Total earliness
        for(int i = 0 ; i < r ; i ++ ){
            obj += (i+1)*processingTime[_sequence[i]][_sequence[i+1]];
        }
        
        //to calculate the Total tardiness
        for(int i = r ; i < numberOfJobs - 1 ; i ++ ){
            obj += (numberOfJobs - i - 1)*processingTime[_sequence[i]][_sequence[i+1]];
        }
    }
    
    public int getgpitimes(){
        return gpitimes;
    }
    
    public void setgpitimes(){
        gpitimes = 1;
    }
    
    public double getObjValue(){
        return obj;
    }
    
    public int[] getSequence(){
        return B;
    }
    
    public int[] getBestSequence(){
        return bestsequence;
    }
    
    public void dumpSequence(int _sequence[]){
        for(int i = 0 ; i < numberOfJobs ; i ++ ){
            bestsequence[i] = _sequence[i];
        }
        
    }
    
    public void backupAP(){
        backupProcessingTime = new double[numberOfJobs][numberOfJobs];
        for(int n = 0 ; n < numberOfJobs ; n ++ ){
            for (int k = 0; k < numberOfJobs ; k++) { //15
                backupProcessingTime[n][k] = processingTime[n][k];
            }
        }
    }
    
    public void restoreAP(){
        for(int n = 0 ; n < numberOfJobs ; n ++ ){
            for (int k = 0; k < numberOfJobs ; k++) { //15
                processingTime[n][k] = backupProcessingTime[n][k];
            }
        }
        for(int i = 0 ; i < numberOfJobs ; i ++ ){
            sequence[i] = i;
        }
    }
    
    
    public void printSequence(int _B[]){
        for(int k = 0 ; k < _B.length ; k ++ ){//15
            System.out.print(_B[k]+" ");
        }
        System.out.print("\n");
    }
    
    public static void main(String[] args) {
        int numberOfJobs;
       double processingTime[][];
        int jobSets[] = new int[]{25};//10, 15, 20, 25
        String type[] = new String[]{"low"};//"low", "med", "high"
        
        for(int replications = 0 ; replications < 1 ; replications ++ ){
            for(int m = 0 ; m < jobSets.length ; m ++ ){//jobSets.length
                for(int n = 0 ; n < type.length ; n ++ ){
                    for(int k = 1 ; k <= 1 ; k ++ ){//15
                        singleMachineSetupData singleMachineData1 = new singleMachineSetupData();
                        SAPT_SingleMachineSetup singleMachine1 = new SAPT_SingleMachineSetup();
                        String fileName = "instances\\SingleMachineSetup\\"+type[n]+"\\"+jobSets[m]+"_"+k+".etp";
                        //fileName = "Data\\SMSetupTime8.txt";//for test
                        //System.out.println(fileName);
                        singleMachineData1.setData(fileName);
                        singleMachineData1.getDataFromFile();
                        numberOfJobs = singleMachineData1.getSize();
                        processingTime = singleMachineData1.getProcessingTime();
                        double bestobj = Double.MAX_VALUE;
                        double obj = Double.MAX_VALUE;
                        int currentSoluion[] = new int[numberOfJobs];
                        homework.util.timeClock timeClock1 = new homework.util.timeClock();
                        timeClock1.start();
                        double times = singleMachine1.getsmallesTimes(processingTime);
                        
                        for(int i = 0 ; i < times ; i ++ ){//i initial solutions
                            if(i == 0){
                                //singleMachine1.setData(numberOfJobs, processingTime);
                                singleMachine1.backupAP();
                            } else{
                                singleMachine1.restoreAP();
                            }
                            
                            singleMachine1.SAPT();
                            singleMachine1.calcObj(singleMachine1.getBestSequence());
                            System.out.println("SAPT_composite " + i);
                            singleMachine1.printSequence(singleMachine1.getBestSequence());
                            String result1 = "best " +"\t"+bestobj+"\n";
                            System.out.println(result1);
                            
                            singleMachine1.GPI(singleMachine1.getBestSequence());
                            singleMachine1.calcObj(singleMachine1.getBestSequence());
                            System.out.println("GPI_composite " + i);
                            singleMachine1.printSequence(singleMachine1.getBestSequence());
                            result1 = "best " +"\t"+bestobj+"\n";
                            System.out.println(result1);
                            System.out.println("obj in main "+singleMachine1.getObjValue());
                            
                            if(bestobj > singleMachine1.getObjValue()){
                                bestobj = singleMachine1.getObjValue();
                                singleMachine1.dumpSequence(singleMachine1.getSequence());
                            }
                        }
                        timeClock1.end();
                        System.out.println("first stategy");
                        singleMachine1.printSequence(singleMachine1.getBestSequence());
                        singleMachine1.calcObj(singleMachine1.getBestSequence());
                        String result1 = "best " + type[n]+"\t"+jobSets[m]+"\t"+k+"\t"+bestobj+"\n";
                        System.out.println(result1);
                        singleMachine1.GPI(singleMachine1.getBestSequence());
                        singleMachine1.calcObj(singleMachine1.getBestSequence());
                        System.out.println("second stategy");
                        singleMachine1.printSequence(singleMachine1.getBestSequence());
                        String result2 = "best " + type[n]+"\t"+jobSets[m]+"\t"+k+"\t"+bestobj+"\n";
                        System.out.println(result2);
                        //System.out.println(singleMachine1.getObjValue());
                        
                        //singleMachine1.writeFile("oneMachineSetup1018", result);
                        
            /*
            System.out.print(replications+" "+obj+" [");
            for(int j = 0 ; j < numberOfJobs ; j ++ ){
              System.out.print((currentSoluion[j]+1)+ " ");
            }
            System.out.print("]\n");
             */
                    }
                }
            }
        }
    }//end main
    
    
    
}