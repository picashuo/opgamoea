/* Copyright 2009-2016 David Hadka
 *
 * This file is part of the MOEA Framework.
 *
 * The MOEA Framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * The MOEA Framework is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the MOEA Framework.  If not, see <http://www.gnu.org/licenses/>.
 */
import org.moeaframework.Executor;
import org.moeaframework.algorithm.NSGAII;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.NondominatedSortingPopulation;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.CrowdingComparator;
import org.moeaframework.core.comparator.ParetoDominanceComparator;
import org.moeaframework.core.operator.GAVariation;
import org.moeaframework.core.operator.InjectedInitialization;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.core.operator.TournamentSelection;
import org.moeaframework.core.operator.real.PM;
import org.moeaframework.core.operator.real.SBX;

/**
 * Demonstrates using an Executor to solve the UF1 test problem with NSGA-II,
 * one of the most widely-used multiobjective evolutionary algorithms.
 */
public class Example1 {

	public static void main(String[] args) {
             Initialization initialization = new RandomInitialization(
                problem,
                100);

        TournamentSelection selection = new TournamentSelection(2,
                new ChainedComparator(
                        new ParetoDominanceComparator(),
                        new CrowdingComparator()));

        Variation variation = new GAVariation(
                new SBX(1.0, 25.0),
                new PM(1.0 / problem.getNumberOfVariables(), 30.0));

        NSGAII algorithm = new NSGAII(
                problem,
                new NondominatedSortingPopulation(),
                null, // no archive
                selection,
                variation,
                initialization);

        while (algorithm.getNumberOfEvaluations() < 10000) {
            algorithm.step();
        }

        Population intermediateResult = algorithm.getPopulation();
        initialization = new InjectedInitialization(
                problem,
                100,
                intermediateResult);

        algorithm = new NSGAII(
                problem,
                new NondominatedSortingPopulation(),
                null, // no archive
                selection,
                variation,
                initialization);

        while (algorithm.getNumberOfEvaluations() < 10000) {
            algorithm.step();
        }

        NondominatedPopulation finalResult = algorithm.getResult();

//		//configure and run this experiment
//                new Executor().withAlgorithm("");
//		NondominatedPopulation result = new Executor()
//				.withProblem("UF11")
//				.withAlgorithm("NSGAII")
//				.withMaxEvaluations(10000)
//				.run();
//		
//		//display the results
//		System.out.format("Objective1  Objective2%n");
//		
//		for (Solution solution : result) {
//			System.out.format("%.4f      %.4f%n",
//					solution.getObjective(0),
//					solution.getObjective(1));
//		}
	}

}
