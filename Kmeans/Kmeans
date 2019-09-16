import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class KMeans {

	// read in data file into string array of arrays

	public static ArrayList<String> readFile(String data) throws IOException{ 
		ArrayList<String> thefile = new ArrayList<>();
		 
		try (BufferedReader br = new BufferedReader(new FileReader(data))) {
		    br.readLine(); //skip header line
			while (br.ready()) {
		        thefile.add(br.readLine());
		    }
		}
		return thefile;
	}
	
	
	public static double[][] makeLargeMatrix(ArrayList<String> thefile){
		//convert array list to 2d matrix with floats
		double[][] matrix = new double[thefile.size()][20]; // change this twenty 1 or 
			for(int i = 0; i < thefile.size(); i++){ 
				String temp = thefile.get(i); // okay
				//System.out.println(temp);
				String[] temparr =  temp.split("\t");
				
				//System.out.println(Arrays.toString(temparr));
				//System.out.println(temparr[18]);
				for(int j = 0; j < 20; j++){	
					
					matrix[i][j] = Double.parseDouble(temparr[j+1]); // oh this is where you do it
					//System.out.println(matrix[i][0]);
					
					
					//}
				}
			}
		return matrix;
		
	}
	
	
	public static String [] rememberNames(ArrayList<String> thefile){
		String [] names = new String [thefile.size()];
		
		for(int i =0; i<names.length; i++){
			names[i]= thefile.get(i).split("\t")[0];
			
		}
		
		return names;
	}
	
	
	
	// normalize the data
	
	private static double [][] Normalizedata(double[][] matrix){
		
		// iterate of col array and make means 
		
		int size = matrix[0].length;
	    double col[] = new double[size];

	    for (int i = 0; i < matrix.length; i++){
	        for (int j = 0; j < matrix[i].length; j++){
	            col[j] += matrix[i][j]/matrix.length;  //sum step
	        }
	    }
		
		// stdv
	    
	    double [] var = new double [matrix[0].length]; 
	    
	    for (int i = 0; i < matrix.length; i++){
	        for (int j = 0; j < matrix[i].length; j++){
	            var[j] += Math.pow(matrix[i][j]-col[j],2)/matrix.length;  //sum step
	        }
	    }
	    
	    double [] stdv = new double [var.length];
	    
	    for(int i=0;i<var.length; i++){
	    	stdv[i] = Math.sqrt(var[i]);
	    }
		
		
	    
	    // Normalize the data
	    
	    
	    double [][] NormMat = new double[matrix.length][matrix[0].length];
	    
	    for (int i = 0; i < matrix.length; i++){
	        for (int j = 0; j < matrix[i].length; j++){
	            NormMat[i][j] = (matrix[i][j]-col[j])/stdv[j];  //sum step
	        }
	    }

		return NormMat;
		
	}
	
	// kmeans stuff
	
	// the number of integers
		public static int Clusters(int clusters) {
			return clusters;
		}
		
	
	
	
	// randomly pick a center of each data within each range
	
	// measure euclidean distance away from each center
	
	// assign 
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	// iterate over each line in the array and only get the int
	
	private static void printArray(Double [][] matrix) {
		for(int i=0; i<matrix.length; i++){
			for(int j=0; j< matrix[0].length; j++ ){
			System.out.println(matrix[i][j]);
			}
		}
	}	
	
	
	

	

	// pick the centroids within the data, must operate on each i and j

	// measure the distance of points from arbitrary centroids

	// cluster the data based on distance

	// update the center of the points
	// Based on what though? cluster mean maybe? i.e. Kmeans?!

	// recurse

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub

		try {
			ArrayList<String> data = readFile("/Users/michaelskaro/Desktop/Assignment_2_BINF8500/Bacteria+Archaea.txt");
			double[][] LargeMatrix = makeLargeMatrix(data);
			String []names = rememberNames(data);
			double [][] datMatrix = Normalizedata(LargeMatrix);
			
			// pick centers of the data
			//  
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} // args[0];
		
		
		//String myFile =  readFileAsString(data);
		//Double[][] matrix = delimitFileInto2Darray(myFile);	
		
		 // this will be the number of clustersint clusters = Integer.parseInt(args[1]);
	
		//Double [] Columnmatrix = columnSum(matrix);
		//double [][] IAA = NormalizeMultidimensionalArray(matrix);
		//printArray(destinationArray);

	}


}
