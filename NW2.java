import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;

public class NW2 {

	// globals

	public static final int gapPenalty = -2;
	public static final int correctMatch = 4;
	public static final int mismatchPenalty = -1;

	// read as string

	public static String readFileAsString(String filename) throws IOException {
		InputStream fasta = new FileInputStream(filename);
		BufferedReader buf = new BufferedReader(new InputStreamReader(fasta));

		String sequence = buf.readLine();
		StringBuilder sb = new StringBuilder();

		sequence = null;
		while ((sequence = buf.readLine()) != null) {
			sb.append(sequence).append("\n");
			// sequence = buf.readLine();
		}

		String fileAsString = sb.toString().toUpperCase();
		fileAsString = fileAsString.replaceAll("\r", "").replaceAll("\n", "");
		buf.close();
		return fileAsString;
	}

	// make score matrix

	public static int[][] makeScoreMatrix(String dnaString1, String dnaString2) {
		int row;
		int col;
		int topLeft;
		int top;
		int left;
		int first;
		int[][] scoreMatrix = new int[dnaString2.length() + 1][dnaString1.length() + 1];

		for (col = 0; col <= dnaString1.length(); col++) {
			scoreMatrix[0][col] = gapPenalty * col;
		}
		for (row = 0; row <= dnaString2.length(); row++) {
			scoreMatrix[row][0] = gapPenalty * row;
		}
		// Now fill in the rest of the array:
		for (row = 1; row <= dnaString2.length(); row++) {
			for (col = 1; col <= dnaString1.length(); col++) {
				if (dnaString1.charAt(col - 1) == dnaString2.charAt(row - 1))
					topLeft = scoreMatrix[row - 1][col - 1] + correctMatch;
				else
					topLeft = scoreMatrix[row - 1][col - 1] + mismatchPenalty;
				left = scoreMatrix[row][col - 1] + gapPenalty;
				top = scoreMatrix[row - 1][col] + gapPenalty;
				first = topLeft;
				if (top > first)
					first = top;
				if (left > first)
					first = left;
				scoreMatrix[row][col] = first;
			}
		}
		return scoreMatrix;
	}

	public static void print3(int x) {

		String s = "" + x;
		if (s.length() == 1)
			System.out.print("  " + s);
		else if (s.length() == 2)
			System.out.print(" " + s);
		else if (s.length() == 3)
			System.out.print(s);
		else
			System.out.print("***");
	}

	public static void printArray(int[][] scoreMatrix) {
		for (int row = 0; row < scoreMatrix.length; row++) {
			for (int col = 0; col < scoreMatrix[row].length; col++)
				print3(scoreMatrix[row][col]);
			System.out.println();
		}
	}

	// make the alignment via traceback on golden path

	public static String[][] alignment(int[][] scoreMatrix, String dnaString1, String dnaString2) {

		// make a matrix that we can iterate over and print out at the end.
		String[][] alignment = new String[3][1];

		StringBuilder outputString1 = new StringBuilder();
		StringBuilder spacer = new StringBuilder();
		StringBuilder outputString2 = new StringBuilder();

		// start at the bottom right corner
		int j = dnaString1.length();
		int i = dnaString2.length();

		while (i > 0 && j > 0) {
			if (scoreMatrix[i - 1][j] >= scoreMatrix[i][j - 1] && scoreMatrix[i - 1][j] >= scoreMatrix[i - 1][j - 1]) {
				

				outputString2.append(dnaString2.charAt(i - 1)); // left
				spacer.append(" ");
				outputString1.append(dnaString1.charAt(j - 1));
				i--;

			} else if (scoreMatrix[i][j - 1] >= scoreMatrix[i - 1][j] && // top
					scoreMatrix[i][j - 1] >= scoreMatrix[i - 1][j - 1]) {
				// take the top

				outputString2.append("-");
				spacer.append("*");
				outputString1.append(dnaString1.charAt(j - 1));
				j--;
			} else if (scoreMatrix[i - 1][j - 1] >= scoreMatrix[i - 1][j] && // match
					scoreMatrix[i - 1][j - 1] >= scoreMatrix[i][j - 1]) {
				// take the diag
				
				outputString2.append(dnaString2.charAt(i - 1));
				spacer.append("|");
				outputString1.append(dnaString1.charAt(j - 1));
				i--;
				j--;	

			}
			
			
		}
		// reverse to prepend order
		String correctOrder1 = outputString1.reverse().toString();
		String correctOrder2 = spacer.reverse().toString();
		String correctOrder3 = outputString2.reverse().toString();

		alignment[0] = new String[] { correctOrder1 };
		alignment[1] = new String[] { correctOrder2 };
		alignment[2] = new String[] { correctOrder3 };
		
		System.out.println(correctOrder1.length());
		System.out.println(correctOrder2.length());
		System.out.println(correctOrder3.length());
		
		return alignment;
	}

	// write out the alignments

	public static void write(String[][] alignment, String output) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		String output1 = Arrays.toString(alignment[0]);
		String spacer = Arrays.toString(alignment[1]);
		String output2 = Arrays.toString(alignment[2]);

		int line = 100;

		for (int i = 0; i < output1.length(); i += line) {
			if (output1.length() - i < line) {
				line = output1.length() - i;
			}

			writer.write(i + ":" + "\n" + output1.substring(i, i + line)); // first
																			// line
																			// of
																			// the
																			// top
			writer.newLine();

			writer.write(spacer.substring(i, i + line));
			writer.newLine();

			writer.write(output2.substring(i, i + line));
			writer.newLine();
			writer.newLine();

		}

	}

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		// read files as strings. Passed as arguements

		String dnaString1 = readFileAsString(args[0]);
		String dnaString2 = readFileAsString(args[1]);

		// make score matrix
		int[][] scoreMatrix = makeScoreMatrix(dnaString1, dnaString2);

		// align the genomes

		String[][] alignment = alignment(scoreMatrix, dnaString1, dnaString2);

		// write the output

		write(alignment, "output.txt");

	}

}
