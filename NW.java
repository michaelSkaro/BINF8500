
public class NW {

	// scores and strings
	public static final int gapPenalty = -2;
	public static final int correctMatch = 4;
	public static final int mismatchPenalty = -1;

	private String dnaStringOne; // First string
	private String dnaStrnigTwo; // Second string
	private int dnaString1Len, dnaString2Len; // their lengths
	private int[][] scoreMatrix;

	// read in file as string arg 1 and arg 2

	public NW(String i, String j) {
		dnaStringOne = i;
		dnaStrnigTwo = j;
		dnaString1Len = dnaStringOne.length();
		dnaString2Len = dnaStrnigTwo.length();
		scoreMatrix = new int[dnaString2Len + 1][dnaString1Len + 1];
	}

	public void makeScoreMatrix() {
		int row;
		int col; 
		int topLeft;
		int top;
		int left; 
		int first; 
		
		for (col = 0; col <= dnaString1Len; col++) {
			scoreMatrix[0][col] = gapPenalty * col;
		}
		for (row = 0; row <= dnaString2Len; row++) {
			scoreMatrix[row][0] = gapPenalty * row;
		}
		// Now fill in the rest of the array:
		for (row = 1; row <= dnaString2Len; row++) {
			for (col = 1; col <= dnaString1Len; col++) {
				if (dnaStringOne.charAt(col - 1) == dnaStrnigTwo.charAt(row - 1))
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
	}
	
	public void print3(int x) {
	    
	    String s = ""+x;
	    if (s.length() == 1) System.out.print("  "+s);
	    else if (s.length() == 2) System.out.print(" "+s);
	    else if (s.length() == 3) System.out.print(s);
	    else System.out.print("***");
	  }

	  public void printArray() {
	    for (int row=0; row < scoreMatrix.length; row++) {
	      for (int col=0; col < scoreMatrix[row].length; col++) 
	        print3(scoreMatrix[row][col]);
	      System.out.println();
	    }
	  }
	
	
	
	

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		NW nw = new NW(args[0], args[1]);
	    nw.makeScoreMatrix();
	    nw.printArray();
	}

}
