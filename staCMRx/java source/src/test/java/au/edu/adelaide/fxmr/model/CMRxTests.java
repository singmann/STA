package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.Random;

import org.junit.Test;

import au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.jet.math.Functions;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class CMRxTests {

	private static final double[][] MODEL_A = { { 1, 0 }, { 0, 1 }, { 1, 1 } };
	private static final double[][] MODEL_B = { { 1, 0 }, { 0, 1 }, { -1, -1 } };

	private static final double[][] MEANS = {
			{ 0.681434511967325, 0.628569652137633, 0.858527246374456, -0.606809499137584, 0.232089352293278,
					-0.296680985874007, 0.170528182305449, 0.834387327659620, 0.514400458221443, -0.239108306049287 },
			{ -0.491435642056938, -0.512950062550021, -0.300032468030383, -0.497832284047938, -0.053422302194541,
					0.661657255792582, 0.099447216582279, -0.428321962359253, 0.507458188556991, 0.135643281450442 },
			{ 0.189998869910387, 0.115619589587611, 0.558494778344073, -1.104641783185522, 0.178667050098737,
					0.364976269918575, 0.269975398887728, 0.406065365300367, 1.021858646778433, -0.103465024598844 } };

	double[][] MODEL_C = { { 1, 1 }, { 0, 1 }, { 1, 1 }, { 0, 1 } };
	double[][] MEANS_C = { { 99.2, 99.8, 99.4, 99.6, 92.7, 99.2, 97.3, 91.4 },
			{ 58.7, 62.7, 52, 49.8, 74.1, 64.6, 64.7, 65.3 }, { 94.2, 99.3, 95.5, 92.9, 61.7, 49.8, 86.4, 55.6 },
			{ 65.6, 48.6, 71.9, 62.4, 88.9, 99, 78, 89.8 } };
	DoubleMatrix2D[] WEIGHTS_C = {
			new DenseDoubleMatrix2D(new double[][] { { 250.000000000000, 0, 0, 0, 0, 0, 0, 0 },
					{ 0, 250.000000000000, 0, 0, 0, 0, 0, 0 }, { 0, 0, 160, 0, 0, 0, 0, 0 },
					{ 0, 0, 0, 160, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0.0990074503106359, 0, 0, 0 },
					{ 0, 0, 0, 0, 0, 250.000000000000, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0.657462195923734, 0 },
					{ 0, 0, 0, 0, 0, 0, 0, 0.0881659282770173 } }),
			new DenseDoubleMatrix2D(new double[][] { { 0.0167965600644988, 0, 0, 0, 0, 0, 0, 0 },
					{ 0, 0.0171461149046890, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0.0239118608807934, 0, 0, 0, 0, 0 },
					{ 0, 0, 0, 0.0160641925132831, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0.0221453287197232, 0, 0, 0 },
					{ 0, 0, 0, 0, 0, 0.0180309320639557, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0.0179546107440391, 0 },
					{ 0, 0, 0, 0, 0, 0, 0, 0.0198411714227608 } }),
			new DenseDoubleMatrix2D(new double[][] { { 0.425124880433627, 0, 0, 0, 0, 0, 0, 0 },
					{ 0, 62.5000000000000, 0, 0, 0, 0, 0, 0 }, { 0, 0, 2.26757369614512, 0, 0, 0, 0, 0 },
					{ 0, 0, 0, 0.324648973297622, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0.0424407686023194, 0, 0, 0 },
					{ 0, 0, 0, 0, 0, 0.0277008310249307, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0.182615047479912, 0 },
					{ 0, 0, 0, 0, 0, 0, 0, 0.0418931515170558 } }),
			new DenseDoubleMatrix2D(new double[][] { { 0.0469131169074873, 0, 0, 0, 0, 0, 0, 0 },
					{ 0, 0.0274115292892190, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0.0532793435984869, 0, 0, 0, 0, 0 },
					{ 0, 0, 0, 0.0660982219578293, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0.730460189919649, 0, 0, 0 },
					{ 0, 0, 0, 0, 0, 81.6326530612245, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0.138408304498270, 0 },
					{ 0, 0, 0, 0, 0, 0, 0, 0.170874449997864 } }) };

	@Test
	public void CMRxTest() {
		DoubleMatrix2D w = DoubleFactory2D.dense.identity(MEANS[0].length).assign(Functions.mult(10));
		DoubleMatrix2D[] weights = { w, w, w };

		CMRxProblem problem = new CMRxProblem(MEANS, weights, null, new DenseDoubleMatrix2D(MODEL_A));
		CMRxSolver solver = new CMRxSolver();
		CMRSolution soln = solver.solve(problem, false);
		assertEquals(0, soln.getFStar(), 1e-15);
		assertEquals(0, soln.getIters().size());

		problem = new CMRxProblem(MEANS, weights, null, new DenseDoubleMatrix2D(MODEL_B));
		for (int i = 0; i < 10; i++)
			soln = solver.solve(problem, false);
		assertEquals(6.04773596513748, soln.getFStar(), 1e-5);
		// assertEquals(67, soln.getIters().size());
	}

	@Test
	public void isFeasible3nTest() {
		TIntHashSet infeas = new TIntHashSet();
		infeas.add(13);
		DoubleMatrix2D tmpVolumes = new DenseDoubleMatrix2D(MEANS[0].length, MEANS[0].length);
		int[][] tmpZoneNumbers = new int[MEANS[0].length][MEANS[0].length];
		int[] if3n = CMRxSolver.isFeasible3n(MEANS, infeas, tmpVolumes, tmpZoneNumbers);
		assertArrayEquals(new int[] { 3, 8, -13 }, if3n);
	}

	@Test
	public void getFeasible6Test() {
		DoubleMatrix2D w = DoubleFactory2D.dense.identity(MEANS[0].length).assign(Functions.mult(10));
		DoubleMatrix2D[] weights = { w, w, w };

		CMRxProblem problem = new CMRxProblem(MEANS, weights, null, new DenseDoubleMatrix2D(MODEL_B));

		CMRxTrial greedy = CMRxSolver.getFeasible6(problem, new MRSolverAJOptimiser(), new TIntObjectHashMap<int[]>());
		assertEquals(6.04773596513748, greedy.getF(), 1e-5);
		
		CMRxTrial greedy7 = CMRxSolver.getFeasible7(problem, new MRSolverAJOptimiser(), new TIntObjectHashMap<int[]>());
		assertEquals(6.04773596513748, greedy7.getF(), 1e-5);
	}

//	@Test
//	public void getFeasible8Test() {
//		DoubleMatrix2D w = DoubleFactory2D.dense.identity(MEANS[0].length).assign(Functions.mult(10));
//		DoubleMatrix2D[] weights = { w, w, w };
//
//		CMRxProblem problem = new CMRxProblem(MEANS, weights, null, new DenseDoubleMatrix2D(MODEL_B));
//
//		CMRxTrial greedy = CMRxSolver.getFeasible8(problem, new MRSolverAJOptimiser(), new TIntObjectHashMap<int[]>());
//		assertEquals(6.97093086063586, greedy.getF(), 1e-5);
//	}

	@Test
	public void sampleDataTest() {
		CMRxProblem problem = new CMRxProblem(MEANS_C, WEIGHTS_C, null, new DenseDoubleMatrix2D(MODEL_C));
		CMRxSolver solver = new CMRxSolver();

		CMRSolution s = solver.solve(problem, false);
		assertEquals(9.54529887661007, s.getFStar(), 1e-5);
	}

	@Test
	public void sampleGetFeasTest() {
		CMRxProblem problem = new CMRxProblem(MEANS_C, WEIGHTS_C, null, new DenseDoubleMatrix2D(MODEL_C));
		CMRxTrial trial = CMRxSolver.getFeasible6(problem, new MRSolverAJOptimiser(), new TIntObjectHashMap<int[]>());
		assertEquals(9.54529898853269, trial.getF(), 1e-5);
		
		trial = CMRxSolver.getFeasible7(problem, new MRSolverAJOptimiser(), new TIntObjectHashMap<int[]>());
		assertEquals(9.54529898853269, trial.getF(), 1e-5);
	}

	// @Test
	public void randomTest() {
		for (int i = 21948; i < 1000000; i++) {
			System.out.println(i);
			try {
				Random rand = new Random(103 * i);
				/// Random rand = new Random(108 * i);

				int nDim = 3;
				int nCond = 4 + rand.nextInt(10);
				// System.out.print(nCond);

				double[][] means = new double[nDim][nCond];
				DoubleMatrix2D[] weights = new DenseDoubleMatrix2D[nDim];

				for (int d = 0; d < nDim; d++) {
					DoubleMatrix2D tmp = new DenseDoubleMatrix2D(nCond, nCond);
					for (int c = 0; c < nCond; c++) {
						means[d][c] = rand.nextDouble() * 100;
						for (int r = 0; r < nCond; r++) {
							tmp.setQuick(r, c, rand.nextDouble() * (r == c ? 10 : 4));
						}
					}

					weights[d] = new DenseDoubleMatrix2D(nCond, nCond);
					weights[d].assign(tmp);
					tmp.zMult(tmp.viewDice(), weights[d]);

					ShrinkDiagonal sd = new ShrinkDiagonal(weights[d]);
					sd.shrink();
					weights[d] = sd.getResult();
				}

				// DoubleMatrix2D model = new DenseDoubleMatrix2D(new double[][]
				// { {
				// 1,
				// 1 }, { 0, 1 }, { 1, 1 }, { 0, 1 } });
				// DoubleMatrix2D model = new DenseDoubleMatrix2D(
				// new double[][] { { 1, 0 }, { 1, 0 }, { 1, 1 }, { 1, 1 }, { 1,
				// 1
				// }, { 1, 0 } });

				DoubleMatrix2D model = new DenseDoubleMatrix2D(new double[][] { { 1, 0 }, { 1, 0 }, { 1, 1 } });
				CMRxProblem problem = new CMRxProblem(means, weights, null, model);

				// System.out.println(problem.toMatlab());
				ParCMRxSolver psolver = new ParCMRxSolver();
				CMRSolution psol = psolver.solve(problem, false);

				CMRxSolver solver = new CMRxSolver();
				CMRSolution sol = solver.solve(problem, false);

				// System.out.println(i + "\t" + sol.getSeconds() + " vs " +
				// psol.getSeconds() + "\t"
				// + sol.getSeconds() / psol.getSeconds());
				//
				// System.out
				// .println("\t" + sol.getFreeMRCalls() + " vs " +
				// psol.getFreeMRCalls() + "\tf=\t" + sol.getFStar());
				//

				assertEquals(sol.getFStar(), psol.getFStar(), 1e-15);
			} catch (Exception e) {
				// We want to know the value of i
				break;
			}
		}
	}

	@Test
	public void evilProblemTest() {
		DoubleMatrix2D model = new DenseDoubleMatrix2D(
				new double[][] { { 1, 1, 1 }, { 1, 0, 0 }, { 1, 1, 0 }, { 1, 0, 1 } });
		DoubleMatrix2D[] weights = new DoubleMatrix2D[4];
		weights[0] = new DenseDoubleMatrix2D(new double[][] {
				{ 66.56026365001046, 40.48796811330176, 45.44202481976386, 60.03610154060311, 49.88393112727002,
						39.64132503754829, 29.64113298383973, 42.58862417332691, 46.55155994389604, 54.54572164461972 },
				{ 40.48796811330176, 42.72996696802675, 38.29836018590341, 40.02524534815892, 39.0824642198973,
						32.70816303458086, 29.02814557608245, 41.09205702856066, 35.90866030255894, 27.81026316601695 },
				{ 45.44202481976386, 38.29836018590341, 80.61817917165084, 80.28191387089844, 60.40428706021969,
						59.63190728450236, 59.87161482854997, 68.73891553029938, 63.50965371123573, 61.72306340468002 },
				{ 60.03610154060311, 40.02524534815892, 80.28191387089844, 108.3403306069031, 68.68778390932324,
						60.1601754460578, 63.43749377924702, 73.79189956779428, 60.6328275744375, 60.16324861624298 },
				{ 49.88393112727002, 39.0824642198973, 60.40428706021969, 68.68778390932324, 83.69318105692609,
						52.82024248053513, 37.13109350233737, 64.17769439149922, 60.56518393007953, 51.16546241904786 },
				{ 39.64132503754829, 32.70816303458086, 59.63190728450236, 60.1601754460578, 52.82024248053513,
						55.96146815353654, 50.96039329918588, 52.24596049262178, 50.24799260149859, 56.09050050045686 },
				{ 29.64113298383973, 29.02814557608245, 59.87161482854997, 63.43749377924702, 37.13109350233737,
						50.96039329918588, 63.6933036552302, 44.13981782952486, 45.40984951192214, 39.25296523759424 },
				{ 42.58862417332691, 41.09205702856066, 68.73891553029938, 73.79189956779428, 64.17769439149922,
						52.24596049262178, 44.13981782952486, 109.6708611452597, 57.53497698931261, 65.97651845631825 },
				{ 46.55155994389604, 35.90866030255894, 63.50965371123573, 60.6328275744375, 60.56518393007953,
						50.24799260149859, 45.40984951192214, 57.53497698931261, 61.21272121504117, 54.75649310706439 },
				{ 54.54572164461972, 27.81026316601695, 61.72306340468002, 60.16324861624298, 51.16546241904786,
						56.09050050045686, 39.25296523759424, 65.97651845631825, 54.75649310706439,
						89.05627354311736 } });
		weights[1] = new DenseDoubleMatrix2D(new double[][] {
				{ 146.2255113184532, 105.1911435544905, 88.69901978243564, 42.7014049843671, 74.31372347100697,
						46.72016685191599, 86.16950401999426, 85.59228509534675, 93.20263004027686, 78.94741746285305 },
				{ 105.1911435544905, 100.8591932761985, 68.43337206115533, 34.31844474636898, 72.22728864599598,
						47.12100275094966, 72.55829654310908, 73.5299315531084, 84.38995231174616, 65.86852572514709 },
				{ 88.69901978243564, 68.43337206115533, 82.79114848663566, 40.37009068684876, 65.01692142282505,
						42.11749397359431, 61.42854779027432, 79.10737006193934, 74.41775238070545, 59.25865037449388 },
				{ 42.7014049843671, 34.31844474636898, 40.37009068684876, 33.44539123845816, 35.45074610592258,
						22.98330203269859, 38.02176757452224, 34.47250334767512, 37.85408212363126, 37.68541190836812 },
				{ 74.31372347100697, 72.22728864599598, 65.01692142282505, 35.45074610592258, 66.59572694330706,
						46.46601312640881, 57.05851104216794, 66.46248755485226, 71.6375378639782, 59.91197240504495 },
				{ 46.72016685191599, 47.12100275094966, 42.11749397359431, 22.98330203269859, 46.46601312640881,
						51.74808502175149, 35.85231425173495, 38.7997747991117, 49.94680070329265, 42.48771781832933 },
				{ 86.16950401999426, 72.55829654310908, 61.42854779027432, 38.02176757452224, 57.05851104216794,
						35.85231425173495, 72.26005429344393, 61.44464083415729, 70.93982162015571, 54.88471326581062 },
				{ 85.59228509534675, 73.5299315531084, 79.10737006193934, 34.47250334767512, 66.46248755485226,
						38.7997747991117, 61.44464083415729, 86.12043464784715, 70.90438318087864, 58.93736021896163 },
				{ 93.20263004027686, 84.38995231174616, 74.41775238070545, 37.85408212363126, 71.6375378639782,
						49.94680070329265, 70.93982162015571, 70.90438318087864, 100.3597116159891, 56.51263779409582 },
				{ 78.94741746285305, 65.86852572514709, 59.25865037449388, 37.68541190836812, 59.91197240504495,
						42.48771781832933, 54.88471326581062, 58.93736021896163, 56.51263779409582,
						66.59500554484548 } });
		weights[2] = new DenseDoubleMatrix2D(new double[][] {
				{ 149.7600290957014, 58.95911605407425, 94.20907889441251, 86.12857766741755, 49.96341550946517,
						50.24890380374627, 86.11866263204963, 91.13441552301141, 53.57522886276648, 49.70804682091084 },
				{ 58.95911605407425, 70.80902264428457, 57.77997208670534, 52.79681132636284, 39.53781338364232,
						36.32386546046854, 64.95208396292362, 60.10769767052074, 54.91502312632731, 65.28635080592821 },
				{ 94.20907889441251, 57.77997208670534, 101.1483869286081, 80.90274584995885, 42.46795708030192,
						41.37491015494496, 80.38783604182214, 76.82812865259049, 57.15518785607539, 73.53944209413524 },
				{ 86.12857766741755, 52.79681132636284, 80.90274584995885, 106.1090541469702, 28.87359800422778,
						49.21278202256198, 85.14893954637449, 67.27478644065276, 55.62124769457762, 69.32795383503027 },
				{ 49.96341550946517, 39.53781338364232, 42.46795708030192, 28.87359800422778, 57.09460129129246,
						27.66391748519222, 30.36219456762311, 40.27174285157558, 46.46164617845989, 32.4228922583656 },
				{ 50.24890380374627, 36.32386546046854, 41.37491015494496, 49.21278202256198, 27.66391748519222,
						32.24690900583943, 40.31760682722386, 40.57631413190061, 35.17377864986883, 35.99237082705593 },
				{ 86.11866263204963, 64.95208396292362, 80.38783604182214, 85.14893954637449, 30.36219456762311,
						40.31760682722386, 96.2258664866953, 72.24619255545431, 52.8935150329863, 63.66912881367284 },
				{ 91.13441552301141, 60.10769767052074, 76.82812865259049, 67.27478644065276, 40.27174285157558,
						40.57631413190061, 72.24619255545431, 78.24623302710491, 48.14634317795563, 66.14646269973171 },
				{ 53.57522886276648, 54.91502312632731, 57.15518785607539, 55.62124769457762, 46.46164617845989,
						35.17377864986883, 52.8935150329863, 48.14634317795563, 62.36281468128915, 58.12305771737548 },
				{ 49.70804682091084, 65.28635080592821, 73.53944209413524, 69.32795383503027, 32.4228922583656,
						35.99237082705593, 63.66912881367284, 66.14646269973171, 58.12305771737548,
						96.64714165459212 } });
		weights[3] = new DenseDoubleMatrix2D(new double[][] {
				{ 94.83499358753383, 75.8975098873823, 81.37960602366442, 58.117258449721, 58.18766900957036,
						49.46296906007306, 59.45926820561287, 83.97407393311707, 41.71226533630053, 49.86502992234135 },
				{ 75.8975098873823, 96.87389009979417, 88.23917384871993, 62.92007655667514, 74.161200346697,
						69.09874630756431, 56.98834133968609, 81.37513188762033, 40.0128237127739, 48.91411613074698 },
				{ 81.37960602366442, 88.23917384871993, 112.3589108898978, 57.84373495028666, 62.95563490411436,
						61.28078700948056, 51.46125198797698, 84.12135185863633, 43.10620635943614, 39.7986378970224 },
				{ 58.117258449721, 62.92007655667514, 57.84373495028666, 55.68163219371755, 53.56809000943632,
						46.85961221725945, 50.9018304526938, 62.58182600978911, 21.26315735517571, 28.39284452199993 },
				{ 58.18766900957036, 74.161200346697, 62.95563490411436, 53.56809000943632, 81.60143866395831,
						56.1437407747348, 48.50954993864221, 66.08879021077152, 29.40163438653731, 42.18259405615963 },
				{ 49.46296906007306, 69.09874630756431, 61.28078700948056, 46.85961221725945, 56.1437407747348,
						56.19851537229361, 39.96046361316065, 55.92229325308152, 29.17464568191973, 33.18302643924641 },
				{ 59.45926820561287, 56.98834133968609, 51.46125198797698, 50.9018304526938, 48.50954993864221,
						39.96046361316065, 67.30231840635734, 62.01817446971955, 18.7910117473262, 42.78363197142505 },
				{ 83.97407393311707, 81.37513188762033, 84.12135185863633, 62.58182600978911, 66.08879021077152,
						55.92229325308152, 62.01817446971955, 86.8098030479389, 37.11807952296604, 57.2257583367211 },
				{ 41.71226533630053, 40.0128237127739, 43.10620635943614, 21.26315735517571, 29.40163438653731,
						29.17464568191973, 18.7910117473262, 37.11807952296604, 29.1151992592851, 22.80454918241622 },
				{ 49.86502992234135, 48.91411613074698, 39.7986378970224, 28.39284452199993, 42.18259405615963,
						33.18302643924641, 42.78363197142505, 57.2257583367211, 22.80454918241622,
						93.36065491860742 } });

		double[][] means = {
				{ 0.9798285438396459, 0.4021927734329805, 0.1435206761271004, 0.6055255149413336, 0.09248416720834529,
						0.9654163711509816, 0.4847120369983434, 0.07570605679334186, 0.2628723777053004,
						0.5747492232270877 },
				{ 0.6650876734694631, 0.03346830589026983, 0.2516508429839454, 0.5733505251109853, 0.8783066555088583,
						0.3480416131449979, 0.8996559426985865, 0.2743314131425028, 0.5934526079714039,
						0.2925381504703778 },
				{ 0.2375011288457212, 0.6606490318311666, 0.1265077464048101, 0.1746023653217372, 0.5036060413390182,
						0.6622222587694857, 0.3461530710520057, 0.7612554964785619, 0.6654920290028329,
						0.3686146898488592 },
				{ 0.9192118613405541, 0.23930548463542, 0.6699545652872549, 0.03739705314017217, 0.8004956914978277,
						0.2404941710632594, 0.7115665316206006, 0.5334601735815366, 0.870030576471695,
						0.2900899670396585 } };

		CMRxProblem p = new CMRxProblem(means, weights, null, model);

		CMRxSolver solver = new CMRxSolver();
		CMRSolution s = solver.solve(p, true);
		// for (HashSet<SimpleLinearConstraint> a : s.getAdjs())
		// System.out.println(a);

		// System.out.println(s.getFStar());
		assertEquals(0.010700599902264537, s.getFStar(), 1e-5);
	}

	@Test
	public void identityTest() {
		DoubleMatrix2D w = DoubleFactory2D.dense.identity(MEANS[0].length).assign(Functions.mult(10));
		DoubleMatrix2D[] weights = { w, w, w };

		CMRxProblem problem = new CMRxProblem(MEANS, weights, null, DoubleFactory2D.dense.identity(3));
		CMRxSolver solver = new CMRxSolver();
		CMRSolution soln = solver.solve(problem, false);
		assertEquals(0, soln.getFStar(), 1e-14);
	}

	@Test
	public void juneDataTest() {
		double[][] means = { { 0.98, 0.95, 0.94, 0.82, 0.89, 0.73, 0.73, 0.96, 0.94, 0.96 },
				{ 0.05, 0.06, 0.05, 0.41, 0.22, 0.06, 0.06, 0.07000000000000001, 0.04, 0.07000000000000001 },
				{ 0.98, 0.94, 0.95, 0.94, 0.9300000000000001, 0.86, 0.86, 0.94, 0.95, 0.94 },
				{ 0.11, 0.17, 0.14, 0.82, 0.55, 0.11, 0.11, 0.26, 0.08, 0.26 } };

		DoubleMatrix2D w = new DenseDoubleMatrix2D(new double[][] { { 10, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 10, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 10, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 10, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 10, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 10, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 10, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 10, 0, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 10, 0 }, { 0, 0, 0, 0, 0, 0, 0, 1, 0, 10 } });
		DoubleMatrix2D[] weights = { w, w, w, w };

		DoubleMatrix2D model1 = new DenseDoubleMatrix2D(new double[][] { { 1 }, { -1 }, { 1 }, { -1 } });
		CMRxProblem problem = new CMRxProblem(means, weights, null, model1);

		ParCMRxSolver solver = new ParCMRxSolver();
		CMRSolution soln = solver.solve(problem, false);
		assertEquals(0.483454217589317, soln.getFStar(), 1e-5);
	}
}
