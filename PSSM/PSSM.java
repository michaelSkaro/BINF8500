import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class PSSM {
	 
	/*
	 * usage: javac PSSM.java 
	 * 			java PSSM arg0 arg1 arg2
	 *
	 * This script will create a position specific scoring matrix in order to 
	 * assess the likelyhood that a biological sequence contains a TSS motif.
	 * 
	 * This script will provide the user with the start, end, and strand the 
	 * motif occurs on. If the motif occurs on the reverse strand the motfi will 
	 * string will always be printed in 5' to 3' direction. Finally this script 
	 * will provide the matching strings of targets.   
	 *
	 *
	 */
	
	 // read search file into String, string new lines
	
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
	
	// read the motifs were are looking to match into an array
	
	public static String [][] readMotifs(String filename) throws IOException{
	
		InputStream motif_file = new FileInputStream(filename);
		BufferedReader buf = new BufferedReader(new InputStreamReader(motif_file));

		String sequence = buf.readLine();
		StringBuilder sb = new StringBuilder();

		sequence = null;
		while ((sequence = buf.readLine()) != null) {
			sb.append(sequence).append("\n");
			// sequence = buf.readLine();
		}

		String fileAsString = sb.toString().toUpperCase();
		fileAsString = fileAsString.replaceAll("\r", "").replaceAll("\n>[A-Z]+", ">").replaceAll("[0-9]+\n", "");
		
		String [] motifs = fileAsString.split(">");
		String [][] motifMatrix = new String [motifs.length][motifs[0].length()];
		
		for(int i =0; i <motifs.length; i++){
			//System.out.println(motifs[i].toString());
			motifMatrix[i] = motifs[i].split("");
			
			
		}

//		System.out.println(motifMatrix[0].length);
//		System.out.println(motifMatrix.length);
//		for(int i=0; i<motifMatrix[0].length; i++){
//			System.out.println(motifMatrix[0][i].getClass());
//		}
//		
		
		
		buf.close();
		return motifMatrix;
	}
	
	
	// make a two dimensional array of scores based on each of the frequency of the 
	// As, Cs, Ts, and Gs in a position.
	
	// make a matrix that is as long a the number of motifs
	// make it as wide as the motifs are
	
	
	// score the frequency of each letter  in each position 
	// move down the letters matrix and count the As, C, Ts, and Gs
	// Fill and update the score for each position 
	
	
	
	public static double [][] makeFreqScoreMatrix(String [][] motifMatrix){
		
		double [][] scoreFreqMatrix = new double [motifMatrix.length][4];		
		
		
		for(int i=0; i<motifMatrix.length; i++){
			for(int j =0; j<motifMatrix[i].length; j++){
				
				if(motifMatrix[i][j].equals("A")){
					scoreFreqMatrix[i][0] =scoreFreqMatrix[i][0] + 1;
				}
				
				else if(motifMatrix[i][j].equals("C")){
					scoreFreqMatrix[i][1]= scoreFreqMatrix[i][1] + 1;
				}
				
				else if(motifMatrix[i][j].equals("T")){
					scoreFreqMatrix[i][2] =scoreFreqMatrix[i][2] + 1;
				}
				
				else if(motifMatrix[i][j].equals("G")){
					scoreFreqMatrix[i][3] =scoreFreqMatrix[i][3] + 1;	
				}

			}
		}
		
		// guard against the 0 frequency with pseudocounts
		// convert to a probability
		// 
		for(int i=0; i<scoreFreqMatrix.length; i++){
			for(int j =0; j<scoreFreqMatrix[i].length; j++){
				// guard against 0
				scoreFreqMatrix[i][j] = scoreFreqMatrix[i][j] +0.01; // try different scores
				// make probability matrix where the max score would be a 100% chance of having 
				// the letter in that position
				scoreFreqMatrix[i][j] = scoreFreqMatrix[i][j]/scoreFreqMatrix[i].length;
				// log probability, 
				//  divide by 0.25 for the As, C, Ts and Gs to include background model
				scoreFreqMatrix[i][j] = Math.log(scoreFreqMatrix[i][j]/0.25); // needs improvement 
				
			}
		}
		
		System.out.print("score matrix made \n\n");
		
		return scoreFreqMatrix;
	}
	

	
	private static void printScoreArray(double[][] scoreFreqMatrix) {
		for(int i=0; i<scoreFreqMatrix.length; i++)
			for(int j=0; j<scoreFreqMatrix[i].length; j++){
			System.out.println(scoreFreqMatrix[i][j]);
			//System.out.println();
			}
		}	
	
	
	
	// okay so I can slide through the string but i am 
	// not sure how to score it against the matrix
	
	
	
	public static ArrayList<String> slidingWindow(String fileAsString, String [][] motifMatrix, double [][] scoreFreqmatrix){
		//ACTG
		// make and array list
		
		// slide in windows of n and increment by 1 until i find a window without sufficient space and then stop
		
		
		// assign each base a score in window based on the position it is in
		
		// score the entire seq
		
		// filter the seqs,
			//if the seq gets a score high enough: keep the position, the score and the string of the window.
		
		ArrayList<String> keeps = new ArrayList<String>();
		
		ArrayList<Integer> positions = new ArrayList<Integer>();
		
		int number =11;
		
		
		
		for(int i=0; i<fileAsString.length(); i++){
			String window = fileAsString.substring(i, number);
			for(int j=0; j< scoreFreqmatrix.length; j++){
				for(int k=0; k<scoreFreqmatrix[j].length; j++){
					// slide, score and compare
					// my window should get score based on the identity of the letter at the position
					
				}
			}
		}
		
		
		
		
		
		return null;
		
	}
	
	
	
	
	public static void main(String[] args) throws IOException {
		
		// read fasta seq into a string
		String search_seq = readFileAsString("/Users/michaelskaro/Desktop/BINF8500/Assignment_4_BINF8500/ecoK12-MG1655.fasta"); 
		// String search_seq = readFileAsString(args[0]); // for after ;)
		
		// read in motifs
		String [][] motifMatrix = readMotifs("/Users/michaelskaro/Desktop/BINF8500/Assignment_4_BINF8500/FruR.fasta");
		
		// String search_seq = readFileAsString(args[1]); // for after ;)
		
		
		// provide the alpha score to keep the motifs
		double alpha = 0.9; //  args[2]); // for after ;)
		
		double [][] scoreFreqMatrix = makeFreqScoreMatrix(motifMatrix);
		
		printScoreArray(scoreFreqMatrix); // all of them are being filled with 0.01's :(
		
		

	}

}
