import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class GS {

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
		System.out.println(fileAsString);
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

		System.out.println("Nucleotide motif matrix");
		for (int i = 0; i < motifMatrix.length; i++) {
			System.out.println(Arrays.toString(motifMatrix[i]));
		}
		System.out.println("\n");

		buf.close();
		return motifMatrix;
	}

	public static double[][] makeBackgroundFreq(String[][] motifMatrix) {

		double[][] backgroundFreq = new double[motifMatrix.length][4];

		for (int i = 0; i < motifMatrix.length; i++) {
			for (int j = 0; j < motifMatrix[i].length; j++) {

				if (motifMatrix[i][j].equals("A")) {
					backgroundFreq[i][0] = backgroundFreq[i][0] + 1;
				}

				else if (motifMatrix[i][j].equals("C")) {
					backgroundFreq[i][1] = backgroundFreq[i][1] + 1;
				}

				else if (motifMatrix[i][j].equals("T")) {
					backgroundFreq[i][2] = backgroundFreq[i][2] + 1;
				}

				else if (motifMatrix[i][j].equals("G")) {
					backgroundFreq[i][3] = backgroundFreq[i][3] + 1;
				}
			}
		}

		System.out.println("Position FreqMatrix");
		for (int i = 0; i < backgroundFreq.length; i++) {
			System.out.println(Arrays.toString(backgroundFreq[i]));
		}
		System.out.println("\n");

		return backgroundFreq;
	}

	// pick the windows that we will be to construct a PSSM

	public static String[][] makeCandidatematrix(String[][] motifMatrix, double windowSize, int min, int max) {

		// make matrix to
		String[][] candidateMatrix = new String[motifMatrix.length][(int) windowSize];

		// iterate over the motif matrix and on each line assign a window
		// randomly
		// assign each letter in that window into a spot in the candidateMatrix
		// This matrix will serve as the PSSM

		for (int i = 0; i < motifMatrix.length; i++) {

			// make a start position
			int randomStart = new Random().nextInt(183);
			int randomEnd = randomStart + 17;
			System.out.println(randomStart + " " + randomEnd + "\n");

			if (randomEnd >= 201) {
				randomStart = randomStart - 17;
				randomEnd = randomEnd - 17;

			}
			ArrayList<String> ar = new ArrayList<String>();
			for (int j = 0; j < motifMatrix[i].length; j++) {
				// in this section will assess whether the i and j
				// term is >= the start position or =< the endposition
				// of the random window

				if (j >= randomStart && j <= randomEnd) {
					ar.add(motifMatrix[i][j]);

					// candidateMatrix[i][j] = motifMatrix[i][j];
					// System.out.print(j + " ");
					System.out.print(motifMatrix[i][j] + " ");
					// candidateMatrix[i][j]=motifMatrix[i][j];
					// candidate j is too big assess the
				}

				if (j > randomEnd && j <= randomEnd + 1) {
					System.out.print("\n\n");
				}

				else {
					// do nothing
				}

				// if this is true we will assign that term to a position in
				// the candidate matrix

			}

			for (int k = 0; k < ar.size(); k++) {
				candidateMatrix[i][k] = ar.get(k);// what the fuck ?
			}

		}

		System.out.println("Candidates Matrix");
		for (int i = 0; i < candidateMatrix.length; i++) {
			System.out.print(Arrays.toString(candidateMatrix[i]) + "\n");
		}
		System.out.println("\n");

		return candidateMatrix;
	}

	public static double[][] makeFreqScoreMatrix(String[][] candidateMatrix) {

		double[][] scoreFreqMatrix = new double[candidateMatrix[0].length][4];

		for (int i = 0; i < candidateMatrix.length; i++) {
			for (int j = 0; j < candidateMatrix[i].length; j++) {

				if (candidateMatrix[i][j].equals("A")) {
					scoreFreqMatrix[j][0] = scoreFreqMatrix[j][0] + 1;
				}

				else if (candidateMatrix[i][j].equals("C")) {
					scoreFreqMatrix[j][1] = scoreFreqMatrix[j][1] + 1;
				}

				else if (candidateMatrix[i][j].equals("T")) {
					scoreFreqMatrix[j][2] = scoreFreqMatrix[j][2] + 1;
				}

				else if (candidateMatrix[i][j].equals("G")) {
					scoreFreqMatrix[j][3] = scoreFreqMatrix[j][3] + 1;
				}
			}
		}

		System.out.println("Position FreqMatrix");
		for (int i = 0; i < scoreFreqMatrix.length; i++) {
			System.out.println(Arrays.toString(scoreFreqMatrix[i]));
		}
		System.out.println("\n");

		// guard against the 0 frequency with pseudocounts
		// convert to a probability
		//
		for (int i = 0; i < scoreFreqMatrix.length; i++) {
			for (int j = 0; j < scoreFreqMatrix[i].length; j++) {
				// guard against 0
				scoreFreqMatrix[i][j] = scoreFreqMatrix[i][j] + 0.25; // try
																		// different
																		// scores
				// probability
				scoreFreqMatrix[i][j] = scoreFreqMatrix[i][j] / scoreFreqMatrix[i].length;
				// log probability,
				scoreFreqMatrix[i][j] = Math.log(scoreFreqMatrix[i][j] / 0.25); // needs
																				// improvement

			}
		}

		// System.out.print("motif scoring matrix made \n\n");

		System.out.println("Scoring position matrix:\n");
		for (int i = 0; i < scoreFreqMatrix.length; i++) {
			System.out.println(Arrays.toString(scoreFreqMatrix[i]));
		}
		System.out.println("\n");

		return scoreFreqMatrix;
	}

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		// window width

		double windowSize = 18;

		// read in the fasta file

		String path = "/Users/michaelskaro/Desktop/BINF8500/Assignment_5_BINF8500/SalmonellaRpoN-sequences-ChIP.fasta";

		// String path = args[1];

		String filename = readFileAsString(path);

		String[][] motifMatrix = readMotifs(path);
		int min = 0;
		int max = motifMatrix[0].length;

		double[][] backgroundFreq = makeBackgroundFreq(motifMatrix);

		String[][] candidatematrix = makeCandidatematrix(motifMatrix, windowSize, min, max);

		double[][] scoreFreqMatrix = makeFreqScoreMatrix(candidatematrix);

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

}
