import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;


/*
* The usage:
* javac GS.java 
* java GS file.fasta 
*
* The objective of this work was to 
* analyze a chip-seq fasta file for a 
* motif sequence that represents a 
* bind zone for a DNA binding protein/
* transcription factor.
* 
*
*/



public class GS {

	// global variable motif indicies
	static int[][] motifIndicies = null;

	public static String readFileAsString(String filename) throws IOException {
		InputStream fasta = new FileInputStream(filename);
		BufferedReader buf = new BufferedReader(new InputStreamReader(fasta));

		String sequence = buf.readLine();
		StringBuilder sb = new StringBuilder();

		sequence = null;
		while ((sequence = buf.readLine()) != null) {
			if (sequence.trim().startsWith(">")) {
				// do nothing
			} else {
				sb.append(sequence).append("\n");
				// sequence = buf.readLine();
			}
		}
		String fileAsString = sb.toString().toUpperCase();
		fileAsString = fileAsString.replaceAll("\r", "").replaceAll("\n", "");
		buf.close();
		// System.out.println(fileAsString);
		return fileAsString;
	}

	public static String[][] readMotifs(String filename) throws IOException {

		InputStream motif_file = new FileInputStream(filename);
		BufferedReader buf = new BufferedReader(new InputStreamReader(motif_file));

		String sequence = buf.readLine();
		StringBuilder sb = new StringBuilder();

		sequence = null;
		while ((sequence = buf.readLine()) != null) {
			if (sequence.trim().startsWith(">")) {
				// do nothing
			} else {
				sb.append(sequence).append("\n");
			}
		}
		String fileAsString = sb.toString().toUpperCase();
		fileAsString = fileAsString.replaceAll("\r", "").replaceAll("\n", "");

		String[] motifs = fileAsString.split("(?<=\\G.{201})");
		String[][] motifMatrix = new String[motifs.length][motifs[0].length()];

		for (int i = 0; i < motifs.length; i++) {
			// System.out.println(motifs[i].toString());
			motifMatrix[i] = motifs[i].split("");

		}

		// System.out.println("Nucleotide motif matrix");
		for (int i = 0; i < motifMatrix.length; i++) {
			// System.out.println(Arrays.toString(motifMatrix[i]));
		}
		// System.out.println("\n");

		buf.close();
		return motifMatrix;
	}

	// get candidate matrix before taking counts

	// pick the windows that we will be to construct a PSSM

	public static String[][] makeCandidatematrix(String[][] motifMatrix, double windowSize, int min, int max) {

		// make matrix to
		String[][] candidateMatrix = new String[motifMatrix.length][(int) windowSize];

		// iterate over the motif matrix and on each line assign a window
		// randomly
		// assign each letter in that window into a spot in the candidateMatrix
		// This matrix will serve as the PSSM

		motifIndicies = new int[motifMatrix.length][2];

		for (int i = 0; i < motifMatrix.length; i++) {

			// make a start position
			int randomStart = new Random().nextInt(183);
			int randomEnd = randomStart + 17;
			// System.out.println(randomStart + " " + randomEnd + "\n");

			if (randomEnd >= 201) {
				randomStart = randomStart - 17;
				randomEnd = randomEnd - 17;

			}

			motifIndicies[i][0] = randomStart;
			motifIndicies[i][1] = randomEnd;

			ArrayList<String> ar = new ArrayList<String>();
			for (int j = 0; j < motifMatrix[i].length; j++) {
				// in this section will assess whether the i and j
				// term is >= the start position or =< the endposition
				// of the random window

				if (j >= randomStart && j <= randomEnd) {
					ar.add(motifMatrix[i][j]);

					// candidateMatrix[i][j] = motifMatrix[i][j];
					// System.out.print(j + " ");
					// System.out.print(motifMatrix[i][j] + " ");
					// candidateMatrix[i][j]=motifMatrix[i][j];
					// candidate j is too big assess the
				}

				if (j > randomEnd && j <= randomEnd + 1) {
					// System.out.print("\n\n");
					// do nothing
				}

				else {
					// do nothing
				}

				// if this is true we will assign that term to a position in
				// the candidate matrix

			}

			for (int k = 0; k < ar.size(); k++) {
				candidateMatrix[i][k] = ar.get(k);
			}

		}

		// System.out.println("Candidates Matrix");
		for (int i = 0; i < candidateMatrix.length; i++) {
			// System.out.print(Arrays.toString(candidateMatrix[i]) + "\n");
		}

		System.out.println("motifIndicies");
		for (int i = 0; i < motifIndicies.length; i++) {
			System.out.print(Arrays.toString(motifIndicies[i]) + "\n");
		}

		// System.out.println("\n");

		return candidateMatrix;
	}

	public static int[][] countResidueoccurences(String[][] motifMatrix, double windowSize, int[][] motifIndicies) {

		int[][] resudieCounts = new int[4][(int) windowSize + 1];

		for (int i = 0; i < resudieCounts.length; i++) {
			for (int j = 0; j < resudieCounts[i].length; j++) {

				resudieCounts[i][j] = 0;
			}

		}

		for (int i = 0; i < motifMatrix.length; i++) {
			for (int j = 0; j < motifMatrix[i].length; j++) {
				if (j < motifIndicies[i][0] || j > motifIndicies[i][1]) {

					if (motifMatrix[i][j].equals("A")) {
						resudieCounts[0][0] = resudieCounts[0][0] + 1;
					}

					else if (motifMatrix[i][j].equals("C")) {
						resudieCounts[1][0] = resudieCounts[1][0] + 1;
					}

					else if (motifMatrix[i][j].equals("T")) {
						resudieCounts[2][0] = resudieCounts[2][0] + 1;
					}

					else if (motifMatrix[i][j].equals("G")) {
						resudieCounts[3][0] = resudieCounts[3][0] + 1;

					}
				}

				else {
					// for (int k = 1; k < windowSize+1; k++) {

					if (motifMatrix[i][j].equals("A")) {
						resudieCounts[0][j - motifIndicies[i][0] + 1]++;
					}

					else if (motifMatrix[i][j].equals("C")) {
						resudieCounts[1][j - motifIndicies[i][0] + 1]++;
					}

					else if (motifMatrix[i][j].equals("T")) {
						resudieCounts[2][j - motifIndicies[i][0] + 1]++;
					}

					else if (motifMatrix[i][j].equals("G")) {
						resudieCounts[3][j - motifIndicies[i][0] + 1]++;

					}

				}

			}
		}

		// System.out.println("residue Counts");
		for (int i = 0; i < resudieCounts.length; i++) {
			// System.out.print(Arrays.toString(resudieCounts[i]) + "\n");
		}
		// System.out.println("\n");

		return resudieCounts;
	}

	// integrate the adjustment

	public static double[][] makeFrequencyCounts(int[][] residueCounts, String[][] motifMatrix) {

		double[][] FreqCounts = new double[residueCounts.length][residueCounts[0].length];

		/*
		 * Algorithm: Window q[i][j]=c[i][j]+b[j]/N-1+B
		 * 
		 * Window q[i][j]=c[i][j]+0.5/4-1+2
		 * 
		 * 
		 * Background q[i][j]=c[0][j]+b[j]/sum(C[0][k])+B
		 * 
		 * 
		 * 
		 */

		// Window

		for (int i = 0; i < residueCounts.length; i++) {
			for (int j = 1; j < residueCounts[i].length; j++) {

				FreqCounts[i][j] = (residueCounts[i][j] + 0.5) / (motifMatrix.length - 1 + 2);

			}
		}

		int sum = 0;

		for (int i = 0; i < residueCounts.length; i++) {
			int j = 0;

			sum = residueCounts[i][j] + sum;
		}

		// Background

		for (int i = 0; i < residueCounts.length; i++) {
			int j = 0;
			FreqCounts[i][j] = (residueCounts[i][j] + 0.5) / (sum + 2);

		}

		// System.out.println("FreqCounts");
		for (int i = 0; i < FreqCounts.length; i++) {
			// System.out.print(Arrays.toString(FreqCounts[i]) + "\n");
		}
		// System.out.println("\n");

		return FreqCounts;
	}

	// make weights for each position in the big matrix

	// remember to input current line at the end

	public static double[] makePositionWeightArray(String[][] motifMatrix, double windowSize, double[][] FreqCounts,
			int currentLine) {

		// double [][] motifsScores = new double
		// [motifMatrix.length][motifMatrix[0].length];
		double[] weights = new double[(int) (motifMatrix[0].length - windowSize)]; // 183
		// System.out.println(weights.length);

		for (int i = 0; i < weights.length; i++) { // 18

			double numerator = 1;
			double denominator = 1;

			for (int j = i; j < windowSize + i; j++) {// 0-17 first round
				// System.out.print(j + ", ");
				if (motifMatrix[currentLine][j].equals("A")) {
					numerator *= FreqCounts[0][j - i + 1];
					denominator *= FreqCounts[0][0];
				}

				else if (motifMatrix[currentLine][j].equals("C")) {
					numerator *= FreqCounts[1][j - i + 1];
					denominator *= FreqCounts[1][0];
				}

				else if (motifMatrix[currentLine][j].equals("T")) {
					numerator *= FreqCounts[2][j - i + 1];
					denominator *= FreqCounts[2][0];
				}

				else if (motifMatrix[currentLine][j].equals("G")) {
					numerator *= FreqCounts[3][j - i + 1];
					denominator *= FreqCounts[3][0];
				}

			}
			// System.out.println();
			// System.out.println(i);
			weights[i] = numerator / denominator;

		}

		// System.out.println("weights");
		for (int i = 0; i < weights.length; i++) {

			// System.out.println(weights[i]);
		}
		// System.out.println("\n");
		return weights;

	}

	// normalize the weights of the weight array

	public static double[] normalizeWeights(double[] weights) {
		double[] normWeights = new double[weights.length];

		double sum = 0;

		for (int i = 0; i < weights.length; i++) {

			sum += weights[i];

		}

		for (int i = 0; i < weights.length; i++) {
			normWeights[i] = weights[i] / sum;
		}
		// System.out.println("Normalized Weights");

		for (int i = 0; i < normWeights.length; i++) {

			// System.out.println(normWeights[i]);
		}
		// System.out.println("\n");
		return normWeights;
	}

	// pick from distribution: gets you value + start index

	// add windowSize +start -1 = Endindex

	// update the motifindicies at the current line

	// make a que with past motifindicies and check after each round for changes
	// (every 5 rounds)

	// iterate until queue

	public static int pickRandomindex(double[] normWeights) {

		// draw random sample from the distribution
		double sum = 0;

		double rand = new Random().nextDouble();

		for (int i = 0; i < normWeights.length; i++) {
			sum += normWeights[i];

			if (sum > rand) {

				// System.out.println();
				// System.out.println(i);
				return i;

			}
		}
		return normWeights.length - 1;

	}

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		// curent line
		int currentLine = 1;

		// window width

		double windowSize = 18;

		// read in the fasta file

		//String path = "/Users/michaelskaro/Desktop/BINF8500/Assignment_5_BINF8500/SalmonellaRpoN-sequences-ChIP.fasta";

		String path = args[0];

		String filename = readFileAsString(path);

		String[][] motifMatrix = readMotifs(path);

		Gibbs(motifMatrix);

		/*
		 * Assign n windows leave one seq out to assess calculate background
		 * frequencies for the search space(left out seq)
		 */

		/*
		 * Use N-1 windows to make PSSM record score for each window
		 */

		/*
		 * Assign probabilities to window locations within seq (evaluate left
		 * out seq) assign a total probability of 0 to length of sequence
		 * between 0 and 1 assign proportional windows
		 */

		/*
		 * Use random number to hit a window that is proportional to the window
		 * probabilities on left out seq repeat for every sequence, incrimenting
		 * over each seq as "forget seq"
		 */

		/*
		 * report areas where highest scores occur report scores.
		 */

	}

	public static void Gibbs(String[][] motifMatrix) {

		// boolean array

		boolean[] hasConverged = new boolean[motifMatrix.length];

		for (int i = 0; i < hasConverged.length; i++) {
			hasConverged[i] = false;
		}

		LinkedList<Integer>[] queue = new LinkedList[motifMatrix.length];
		int historySize = 7;

		for (int i = 0; i < motifMatrix.length; i++) {
			queue[i] = new LinkedList<Integer>();
		}
		double windowSize = 18;
		int min = 0;
		int max = motifMatrix[0].length;
		String[][] candidatematrix = makeCandidatematrix(motifMatrix, windowSize, min, max);

		while (!checkStatus(hasConverged)) {

			for (int i = 0; i < motifMatrix.length; i++) {
				
					// current line
					int currentLine = i;
					
					int[][] residueCounts = countResidueoccurences(motifMatrix, windowSize, motifIndicies);

					double[][] FreqCounts = makeFrequencyCounts(residueCounts, motifMatrix);

					double[] weights = makePositionWeightArray(motifMatrix, windowSize, FreqCounts, currentLine);

					double[] normWeights = normalizeWeights(weights);

					int draw = pickRandomindex(normWeights);

					// draw can be out of bounds
					motifIndicies[i][0] = draw;

					if (draw > motifMatrix[0].length - 18) {
						motifIndicies[i][1] = draw + (201 - draw);
					} else {
						motifIndicies[i][1] = draw + 18 - 1;
					}

					if (queue[i].size() < historySize) {
						queue[i].addLast(draw);
						//System.out.print("hit if statement ");
					}

					else {
						//System.out.print("else statement \n");
						queue[i].remove();
						queue[i].addLast(draw);
						boolean converged = true;
						int value = queue[i].peek();
						for (int j = 0; j < historySize; j++) {
							if (queue[i].get(j) != value) {
								converged = false;
							}
						}


						hasConverged[i] = converged;
					}

				}
			
			
			// print out has converged matricies

			for (int i = 0; i < hasConverged.length; i++) {
				//System.out.println(hasConverged[i]);

			}
			//System.out.println();

			// print queue
			for (int i = 0; i < queue.length; i++) {
				System.out.println(queue[i]);
				for (int j = 0; j < queue[i].size(); j++) {

					System.out.print(queue[i].get(j) + ", ");
				}
				System.out.println();
			}
			System.out.println();
			
			for(int i =0; i<hasConverged.length; i++){
				System.out.println(hasConverged[i]);
			}
			
		}

		

		// System.out.println("motifIndicies");
		for (int i = 0; i < motifIndicies.length; i++) {
			System.out.print(Arrays.toString(motifIndicies[i]) + "\t");
			
			for(int j =motifIndicies[i][0]; j<motifIndicies[i][1]; j++){
				System.out.print(motifMatrix[i][j]);
				
			}
			
			System.out.println();
		}

	}

	public static boolean checkStatus(boolean[] hasConverged) {

		boolean status = true;

		for (int i = 0; i < hasConverged.length; i++) {
			if (hasConverged[i] != true) {
				status = false;
				return status;
			}

		}

		return status;
	}
}
