
/**
 * QuickSort on a Fastq in Java
 *
 * Initialize the program by compiling the code.
 * $ javac QuickSort_v3.java  
 *
 * Run the code the compiled java code
 * $ java QuickSort_v3 /path_to_fastq_to_Quicksort/
 *
 *
 * Please note that you will get an index out of bounds if you 
 * do not provide the fastq file
 *
 *
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.regex.*;

public class QuickSort_v3 {
	
	
	// read the file into a string
	
	public static String readFileAsString(String fastq_file) {
		String fastq = "";
		try {
			fastq = new String (Files.readAllBytes(Paths.get(fastq_file)));
			//System.out.println(fastq);
		} catch (IOException e) {
			e.printStackTrace();
		}

		return fastq;
	}
	
	// delimit the string into an array
	// may want to make this  an array of arrays, maybe?
	
	public static String[][] delimitFileinotArray(String fastq) throws Exception{
		
		// split the file
		
		String delimiter = "\n@ERR";
		
		String[] Fastq_seqs = fastq.split(delimiter);
		
		
		
		//System.out.println(Fastq_seqs.length);
		
		
		String [][] saa = new String [Fastq_seqs.length][2];
		
		Pattern p = Pattern.compile("[A-Z]{44}");
		Matcher m;
		for(int i = 0; i<saa.length; i++){
			saa[i][0]=Fastq_seqs[i];
			//print 2d array
			m = p.matcher(Fastq_seqs[i]);
			
			if(m.find()){
				saa[i][1] = m.group();
				//System.out.println(saa[i][1]);
			}
			
		
		}
		
		for(int i= 0; i<saa.length; i++){
			if(saa[i][1] == null){
				//System.out.println(i);
			}
		}
		
		
		for (int i=0; i<3; i++) {
		    //System.out.println(saa[i][1] + "\n");
		    //System.out.println(saa[i][0] + "\n");
		
		}
		
		//return Fastq_seqs;
		return saa;
		
	}
	
	
	public static int partition(String[][] saa, int indexStart, int indexEnd){ 
		String pivotValue = saa[indexEnd][1];
		int i = indexStart - 1;
		int j = indexEnd + 1;
		
		while(true){
			do {
				i++;
			}
			while (saa[i][1].compareTo(pivotValue) < 0);
		
			do {
				j--;
			} while(saa[j][1].compareTo(pivotValue) > 0);
				
			//System.out.println(saa[indexStart][0] + "\n\n\n" + "before swap");
			if(i<j){
			String [] temp = saa[i];
	        //swap
	        // am i ever actually doing anything to these values?
	        //System.out.println("made it to swap");
	        saa[i] = saa[j];
	        //System.out.println(saa[indexStart][0] + "\n\n\n" + "after swap");
	        saa[j] = temp;
			} else {
				return i;
			}
		}
		
	}
	
	
	public static void QuicksortFastqReads(String[][] saa, int indexStart, int indexEnd) {
		if (indexStart >= indexEnd) {
			return;
		}
		int pivotIndex = partition(saa, indexStart, indexEnd);
		//System.out.println(pivotIndex);
		QuicksortFastqReads(saa, indexStart, pivotIndex - 1);
		QuicksortFastqReads(saa, pivotIndex, indexEnd);
		//We have not made it
		//System.out.print("\n" + "We have made it out of QuicksortFastqReads");
		
	}
	
	
	public static String [][] QuicksortString(String[][] saa) {
		QuicksortFastqReads(saa, 0, saa.length - 1);
		return saa;
	}
	
	
	private static void printArray(String[][] saa) {
		for(int i=0; i<saa.length; i++)
			System.out.println("@ERR" + saa[i][0]);
			
		}	

	
	public static void write (String output, String [][]saa ) throws IOException{
		  BufferedWriter outputWriter = null;
		  outputWriter = new BufferedWriter(new FileWriter(output));
		  for (int i = 0; i < saa.length; i++) {
		    // Maybe:
		    outputWriter.write("@ERR" + saa[i][0]+"\n");
		  }
		  outputWriter.flush();  
		  outputWriter.close();  
		}
	
	
	public static void main(String[] args) throws Exception {
		/*
		 * String[][] saa = { {"ZD"}, {"DC"}, {"AB"}, {"AR"}, {"TB"}, {"YY"}, {"UC"}, {"VX"} };
		 * printArray(saa);
		 *QuicksortString(saa);
		 *System.out.print("\n");
		 *printArray(saa);
		 */
		
		
		String fastq_file = args[0];
		fastq_file = readFileAsString(fastq_file);
		String [][] saa = delimitFileinotArray(fastq_file);
		QuicksortString(saa);
		write("output.txt", saa);
		//System.out.print("done");
		//printArray(saa);
		
		//System.out.print("\n" + fastq_file.substring(1,150));
		
	}
}
