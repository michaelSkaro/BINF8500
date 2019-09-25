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
import java.util.Random;


/*

I have attached both files in the email.

PSA:
The data in this file is passed in the main. I am not sure what was happening on the cluster but it 
would not accept the file or the clusters as arguments today. sorry


********************
Usage: javac KMeans.java 
		java KMeans
clusters:10
file: Bacteria+Archea.txt
********************


I have decided to use the AIC as the model predictor and it is rrported to the user after 
every iteration for each cluster. The stabilized cluster WCSS is reported and the 
cluster assignment is reported. I ran the code 300 times to find the global minimum. 





*/






public class KMeans {

	// read in data file into string array of arrays

	public static ArrayList<String> readFile(String data) throws IOException {
		ArrayList<String> thefile = new ArrayList<>();

		try (BufferedReader br = new BufferedReader(new FileReader(data))) {
			br.readLine(); // skip header line
			while (br.ready()) {
				thefile.add(br.readLine());
			}
		}
		return thefile;
	}

	public static double[][] makeLargeMatrix(ArrayList<String> thefile) {
		// convert array list to 2d matrix with floats
		double[][] matrix = new double[thefile.size()][20]; // change this
															// twenty 1 or
		for (int i = 0; i < thefile.size(); i++) {
			String temp = thefile.get(i); // okay
			// System.out.println(temp);
			String[] temparr = temp.split("\t");

			// System.out.println(Arrays.toString(temparr));
			// System.out.println(temparr[18]);
			for (int j = 0; j < 20; j++) {

				matrix[i][j] = Double.parseDouble(temparr[j + 1]); // oh this is
																	// where you
																	// do it
				// System.out.println(matrix[i][0]);

				// }
			}
		}
		return matrix;

	}

	public static String[] rememberNames(ArrayList<String> thefile) {
		String[] names = new String[thefile.size()];

		for (int i = 0; i < names.length; i++) {
			names[i] = thefile.get(i).split("\t")[0];

		}

		return names;
	}

	// normalize the data

	private static double[][] Normalizedata(double[][] matrix) {

		// iterate of col array and make means

		int size = matrix[0].length;
		double col[] = new double[size];

		for (int i = 0; i < matrix.length; i++) { // ++130
			for (int j = 0; j < matrix[i].length; j++) { // ++20
				col[j] += matrix[i][j] / matrix.length; // mean step
			}
		}

		// stdv

		double[] var = new double[matrix[0].length]; // ++20

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				var[j] += Math.pow(matrix[i][j] - col[j], 2) / matrix.length; // sum
																				// step
			}
		}

		double[] stdv = new double[var.length];

		for (int i = 0; i < var.length; i++) { // ++20
			stdv[i] = Math.sqrt(var[i]);
		}

		// Normalize the data

		double[][] NormMat = new double[matrix.length][matrix[0].length];

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				NormMat[i][j] = (matrix[i][j] - col[j]) / stdv[j]; // sum step
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

	// k by 20 matrix
	private static double[][] makeCentroids(double[][] NormMat, double clusters) {
		Random rand = new Random();
		// int size = NormMat[0].length;

		double[][] centroids = new double[(int) clusters][NormMat[0].length]; // think
																				// maybe
																				// reversed
		for (int i = 0; i < clusters; i++) {
			for (int j = 0; j < NormMat[0].length; j++) {
				centroids[i][j] = 4 * rand.nextDouble() - 2;

				// System.out.println(centroids[i][j]);// -2:2
			}
			// System.out.println(centroids[i]);
		}

		return centroids;

	}

	// k distances
	// measure euclidean distance away from each center

	private static double[][] measureDistances(double[][] NormMat, double[][] centroids, double clusters) {

		double[][] distanceArray = new double[NormMat.length][(int) clusters]; // 130x10
		// System.out.println(NormMat.length);

		for (int i = 0; i < NormMat.length; i++) { // ++ 130
			for (int j = 0; j < distanceArray[0].length; j++) {// ++ 3
				double sum = 0;
				for (int k = 0; k < NormMat[i].length; k++) { // ++20
					sum = Math.pow(NormMat[i][k] - centroids[j][k], 2);
				}
				distanceArray[i][j] = Math.sqrt(sum);
			}
		}
		return distanceArray;
	}

	// assign clusters based on distances away from centers

	private static int[] assignClusters(double[][] distanceArray) {

		// iterate over the distance matrix and assign names of genomes into
		// clusters based on their index and
		// line distance away from centroid.
		int[] clusterAsssignments = new int[distanceArray.length];

		for (int i = 0; i < distanceArray.length; i++) { // ++10 // this is the
															// problem!!!
			int min = 0;
			for (int j = 1; j < distanceArray[i].length; j++) {
				if (distanceArray[i][j] < distanceArray[i][min]) {
					min = j;
				}
				for (double val : distanceArray[i]) {
					// System.out.print(val + ", ");
				}
				// System.out.println("\n" + min);

			}
			clusterAsssignments[i] = min;
		}

		return clusterAsssignments;

	}

	// pick centroids within clusters: mean of cluster

	private static double[][] assignCentroids(double[][] NormMat, int[] clusterAsssignments, double clusters) {

		double[][] centroids = new double[(int) clusters][NormMat[0].length];

		for (int i = 0; i < centroids.length; i++) { // 10 deep here
			int clusterSize = 0; // calculate size of the cluster
			for (int j = 0; j < NormMat.length; j++) {
				if (clusterAsssignments[j] == i) {
					clusterSize++;
				}
			}
			for (int j = 0; j < NormMat.length; j++) {

				for (int k = 0; k < NormMat[0].length; k++) {
					if (clusterAsssignments[j] == i) {
						centroids[i][k] += NormMat[j][k] / clusterSize;
					}
				}

			}

		}
		return centroids;
	}

	private static double[] computerSumofSquares(double[][] NormMat, double[][] centroids, double clusters,
			int[] clusterAsssignments) {

		double[] wcss = new double[(int) clusters]; // 3x1

		for (int i = 0; i < NormMat.length; i++) { // ++ 130
			for (int j = 0; j < NormMat[i].length; j++) {

				wcss[clusterAsssignments[i]] += Math.pow(centroids[clusterAsssignments[i]][j], 2);
			}

		}
		return wcss;
	}

	
	private static double [] computeAkikeinfromationCriterion(double[] wcss, double clusters){
		double [] aic = new double[wcss.length];
			
		for(int i=0; i<wcss.length; i++){
			aic[i] = (2*clusters)+wcss[i];
		}
		
		return aic;
		
	}
	
	
	private static void reportAIC(double [] aic, double clusters) {
		for (int i = 0; i < aic.length; i++) {
			System.out.println("The number of clusters: " + clusters + "\n" +
								"The aic: " +aic[i]);
		}
	}
	
	private static void printArray(double[] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			System.out.print("This is the stable cluster WCSS for cluster: " +i + "\t" + matrix[i] + ", " + "\n");
		}
		System.out.println();

	}

	private static void reportClusters(String[] names, int[] clusterAsssignments) {
		System.out.println("Genome:\t" + "Cluster Assignent:\n");
		for (int i = 0; i < names.length; i++) {
			System.out.println(names[i] + "\t" + clusterAsssignments[i]);
		}
	}

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub

		try {
			ArrayList<String> data = readFile("Bacteria+Archaea.txt");
			double[][] LargeMatrix = makeLargeMatrix(data);
			String[] names = rememberNames(data);
			double[][] datMatrix = Normalizedata(LargeMatrix);
			int clusters = 10; // Integer.parseInt(args[0]);
			for(int i=0; i<300; i++){
			KMeansAlg(datMatrix, clusters, names);
			}
			// pick centers of the data
			//

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} // args[0];

	}

	public static void KMeansAlg(double[][] datMatrix, int clusters, String[] names) {
		// Clusters

		Clusters(clusters);

		// pick centroids

		double[][] centroids = makeCentroids(datMatrix, clusters);

		double[] oldWcss = new double[(int) clusters];

		for (double x : oldWcss) {
			x = 0.0;
		}

		boolean keepGoing = true;
		int[] clustersAssignments;
		

		do {

			// measure distances

			double[][] distances = measureDistances(datMatrix, centroids, clusters);

			// assign clusters

			clustersAssignments = assignClusters(distances);
			
			/// calculate the wcss

			centroids = assignCentroids(datMatrix, clustersAssignments, clusters);


			double[] wcss = computerSumofSquares(datMatrix, centroids, clusters, clustersAssignments);
			
			double [] aic = computeAkikeinfromationCriterion(wcss, clusters);
			
			System.out.println("\n" + "New cycle to stablize cluster"+ "\n");
			reportAIC(aic, clusters);
			
			// keep going

			keepGoing = false;

			for (int i = 0; i < wcss.length; i++) {
				if (wcss[i] - oldWcss[i] > 1.0 * 10e-30) {
					keepGoing = true;
				}
			}
			oldWcss = wcss;

			//printArray(oldWcss);

		} while (keepGoing);

		reportClusters(names,clustersAssignments);
		
		
		printArray(oldWcss);

	}

}
