import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class NW2 {

	// globals for the hell of it
	
	public static final int gapPenalty = -2;
	public static final int correctMatch = 4;
	public static final int mismatchPenalty = -1;
	private static int[][] scoreMatrix;
	
	
	
	// read as string
	
	public static String readFileAsString(String filename) throws IOException { 
		InputStream fasta = new FileInputStream(filename);
		BufferedReader buf = new BufferedReader(new InputStreamReader(fasta));

		String sequence = buf.readLine();
		StringBuilder sb = new StringBuilder();

		sequence = null;
		while ((sequence = buf.readLine()) != null) {
			sb.append(sequence).append("\n");
			sequence = buf.readLine();
		}

		String fileAsString = sb.toString().toUpperCase();
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
		int[][] scoreMatrix = new int[dnaString2.length() +1][dnaString1.length() +1];

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
	
	
	// make the alignment via traceback on golden path
	
	
	public static String[][] alignment(int[][] scoreMatrix, String dnaString1,String dnaString2) {
		
		// make a matrix that we can iterate over and print out at the end.
		String[][] alignments = new String[3][1];
		
		StringBuilder outputString1 = new StringBuilder();
		StringBuilder spacer = new StringBuilder();
		StringBuilder outputString2 = new StringBuilder();

		// start at the bottom right corner
		int j = dnaString1.length();
		int i = dnaString2.length();

		
		while (i > 0 || j > 0) {
			if (scoreMatrix[i - 1][j] >= scoreMatrix[i][j - 1] && scoreMatrix[i - 1][j] >= scoreMatrix[i - 1][j - 1]) {
				// take the left
				
				outputString2.append(dnaString2.charAt(i - 1));
				spacer.append("|");
				outputString1.append(dnaString1.charAt(j - 1));
				i--;

			} else if (scoreMatrix[i][j - 1] >= scoreMatrix[i - 1][j] && // top
					scoreMatrix[i][j - 1] >= scoreMatrix[i - 1][j - 1]) {
				// take the top
				
				outputString2.append("-");
				spacer.append(" ");
				outputString1.append(dnaString1.charAt(j - 1));
				j = j - 1;
			} else if (scoreMatrix[i - 1][j - 1] >= scoreMatrix[i - 1][j] && // match
					scoreMatrix[i - 1][j - 1] >= scoreMatrix[i][j - 1]) {
				// take the diag
				
				outputString2.append(dnaString2.charAt(i - 1));
				spacer.append("*");
				outputString1.append(dnaString1.charAt(j - 1));
				i--;
				j--;

			}
		}
		return alignments;
	}
	
	public void print3(int x) {

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
	
	public void printArray() {
		for (int row = 0; row < scoreMatrix.length; row++) {
			for (int col = 0; col < scoreMatrix[row].length; col++)
				print3(scoreMatrix[row][col]);
			System.out.println();
		}
	}
	
	
	
	

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		// read file as string

		String dnaString1 = readFileAsString(
				"/Users/michaelskaro/Desktop/Assignment_3_BINF8500/FastaSampleFiles/HIV1a.fasta");
		String dnaString2 = readFileAsString(
				"/Users/michaelskaro/Desktop/Assignment_3_BINF8500/FastaSampleFiles/HIV1b.fasta");

		// make score matrix
		
		int [][] scoreMatrix = makeScoreMatrix(dnaString1, dnaString2);

	}

}
