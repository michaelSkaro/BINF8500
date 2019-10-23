import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;

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

		System.out.println("Nucleotide motif matrix");
		for (int i = 0; i < motifMatrix.length; i++){
		System.out.println(Arrays.toString(motifMatrix[i]));
		}
		System.out.println("\n");
		
		
		buf.close();
		return motifMatrix;
	}
	
	
	
	
	public static double [][] makeFreqScoreMatrix(String [][] motifMatrix){
		
		double [][] scoreFreqMatrix = new double [motifMatrix[0].length][4];

			for(int i=0; i<motifMatrix.length; i++){
				for(int j =0; j<motifMatrix[i].length; j++){
					
					
					if(motifMatrix[i][j].equals("A")){
						scoreFreqMatrix[j][0] =scoreFreqMatrix[j][0] + 1;
					}
					
					else if(motifMatrix[i][j].equals("C")){
						scoreFreqMatrix[j][1]= scoreFreqMatrix[j][1] + 1;
					}
					
					else if(motifMatrix[i][j].equals("T")){
						scoreFreqMatrix[j][2] =scoreFreqMatrix[j][2] + 1;
					}
					
					else if(motifMatrix[i][j].equals("G")){
						scoreFreqMatrix[j][3] =scoreFreqMatrix[j][3] + 1;	
						}
					}
				}
		
		System.out.println("Position FreqMatrix");
		for (int i = 0; i < scoreFreqMatrix.length; i++){
		System.out.println(Arrays.toString(scoreFreqMatrix[i]));
		}
		System.out.println("\n");
		
		// guard against the 0 frequency with pseudocounts
		// convert to a probability
		// 
		for(int i=0; i<scoreFreqMatrix.length; i++){
			for(int j =0; j<scoreFreqMatrix[i].length; j++){
				// guard against 0
				scoreFreqMatrix[i][j] = scoreFreqMatrix[i][j] +0.25; // try different scores
				// probability
				scoreFreqMatrix[i][j] = scoreFreqMatrix[i][j]/scoreFreqMatrix[i].length;
				// log probability, 
				scoreFreqMatrix[i][j] = Math.log(scoreFreqMatrix[i][j]/0.25); // needs improvement 
				
			}
		}
		
		//System.out.print("motif scoring matrix made \n\n");
		
		
		System.out.println("Scoring position matrix:\n");
		for (int i = 0; i < scoreFreqMatrix.length; i++){
		System.out.println(Arrays.toString(scoreFreqMatrix[i]));
		}
		System.out.println("\n");
		
		return scoreFreqMatrix;
	}
	
	// score each of the motifs and determine a filter based on the lowest score
	
	
	public static double scoreMotifs(double [][] scoreFreqMatrix, String [][] motifMatrix){
		
		double [][] motifsScores = new double [motifMatrix.length][motifMatrix[0].length];
		double [] rowScores = new double[motifsScores.length];
		double lowestMotifScore =0;
		
		for(int i=0; i<motifMatrix.length; i++){
			for(int j=0; j< motifMatrix[0].length; j++){
					if(motifMatrix[i][j].equals("A")){
						motifsScores[i][j] =scoreFreqMatrix[j][0];
					}
					
					else if(motifMatrix[i][j].equals("C")){
						motifsScores[i][j] =scoreFreqMatrix[j][1];
					}
					
					else if(motifMatrix[i][j].equals("T")){
						motifsScores[i][j] =scoreFreqMatrix[j][2];
					}
					
					else if(motifMatrix[i][j].equals("G")){
						motifsScores[i][j] =scoreFreqMatrix[j][3];	
					}
			}	
		}
		
		// Scores of motifs
		System.out.println("Scoring motifs matrix:\n");
		for (int i = 0; i < motifsScores.length; i++){
		System.out.println(Arrays.toString(motifsScores[i]));
		}
		System.out.println("\n");
		
		
		
		// score each row
		
		for(int i=0; i<motifsScores.length; i++){
			double sum=0;
			for(double j: motifsScores[i]){
				sum +=j;
				rowScores[(int) i] = sum;
				
			}
		}
		
		
		double minValue = rowScores[0]; 
	    for(int i=1;i<rowScores.length;i++){ 
	      if(rowScores[i] < minValue){ 
	        minValue = rowScores[i]; 
	      } 
	      
	    }
	    lowestMotifScore =minValue;
	    System.out.println("alpha value");
	    System.out.println(lowestMotifScore);
		return lowestMotifScore;
	}
	
	// slide through the genome and score each
	
	public static void scoreWindows(String search_seq, double lowestMotifScore, double [][] scoreFreqmatrix, String [][] motifMatrix){
		double [][] motifsScores = new double [motifMatrix.length][motifMatrix[0].length];
		double windowWidth = motifMatrix[0].length;
		String [] window = new String[(int) windowWidth];
		String candidate;
		
		double [] candidateScores = new double [motifMatrix.length];
		
		for(int i=0; i<search_seq.length() - windowWidth; i++){
			candidate = (String) search_seq.subSequence(i, i+(int) windowWidth);
			//System.out.print(candidate);
			//System.out.println(candidate);
			
			window= candidate.split("");
			
			
			
			
			for(int j=0; j< window.length; j++){
				if(window[j].equals("A")){
					candidateScores[j] =scoreFreqmatrix[j][0];
				}
				
				else if(window[j].equals("C")){
					candidateScores[j] =scoreFreqmatrix[j][1];
				}
				
				else if(window[j].equals("T")){
					candidateScores[j] =scoreFreqmatrix[j][2];
				}
				
				else if(window[j].equals("G")){
					candidateScores[j] =scoreFreqmatrix[j][3];	
				}
			}
		
			
			double sum=0;
			
			for(double k: candidateScores){
					sum +=k;
					
			}	
			
			// make an array of strings of the keeps
			//System.out.println("Keeps:\n");
			if(sum >=lowestMotifScore){
				
				//System.out.println(window.toString() + i + (i+window.length));
				
				System.out.println(candidate + " Start: " + i + " " + " End: " + (i+window.length) + " Strand +");
				
				
				
			}
				
		}
		
		
		for(int i=search_seq.length(); i > windowWidth; i--){
			candidate = (String) search_seq.subSequence(i-(int) windowWidth,i);
			//System.out.print(candidate);
			//System.out.println(candidate);
			
			window= candidate.split("");
			
			for(int j=0; j< window.length; j++){
				if(window[j].equals("A")){
					candidateScores[j] =scoreFreqmatrix[j][0];
				}
				
				else if(window[j].equals("C")){
					candidateScores[j] =scoreFreqmatrix[j][1];
				}
				
				else if(window[j].equals("T")){
					candidateScores[j] =scoreFreqmatrix[j][2];
				}
				
				else if(window[j].equals("G")){
					candidateScores[j] =scoreFreqmatrix[j][3];	
				}
			}
		
			
			double sum=0;
			
			for(double k: candidateScores){
					sum +=k;
					
			}	
			
			// make an array of strings of the keeps
			//System.out.println("Keeps:\n");
			if(sum >=lowestMotifScore){
				
				//System.out.println(window.toString() + i + (i+window.length));
				String temp= "";
				for(int p=candidate.length() -1; p>0; p--){
					temp += candidate.charAt(p);
					
				}
				candidate =temp;
				System.out.println(temp + " Start: " + i + " " + " End: " + (i+window.length) + " Strand -");
				
				
				
			}
				
		}
			
		
			

	}
	
	
	
	
	
	
	
	public static void main(String[] args) throws IOException {
		
		// read fasta seq into a string
		//String search_seq = readFileAsString("/Users/michaelskaro/Desktop/BINF8500/Assignment_4_BINF8500/Scel-So0157-2.fasta"); 
		String search_seq = readFileAsString(args[0]); // for after ;)
		
		// read in motifs
		//String [][] motifMatrix = readMotifs("/Users/michaelskaro/Desktop/BINF8500/Assignment_4_BINF8500/sigma54.fasta");
		
		String [][] motifMatrix = readMotifs(args[1]); // for after ;)
		
		
		// provide the alpha score to keep the motifs
		
		
		double [][] motifScoreMatrix = makeFreqScoreMatrix(motifMatrix);
		
		
		// score each of the known motifs
		
		double lowestMotifScore = scoreMotifs(motifScoreMatrix, motifMatrix);
		
		// score each window
		
		scoreWindows(search_seq, lowestMotifScore, motifScoreMatrix, motifMatrix);
		
		

	}

}
