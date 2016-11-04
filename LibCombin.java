public class LibCombin {
	static int[][] choose; // choose[n][k] = nCk
	/*
	 * Computes all combinations of the form aCb for 0 <= b <= a <= n.
	 * Takes advantage of arrays in Java being initialized to contain zeroes.
	 */
	static void initChooseDP(int n) {
		choose = new int[n+1][n+1];
		for (int i=0; i<=n; i++)
			choose[i][0] = 1;
		for (int i=1; i<=n; i++)
			for (int j=1; j<=i; j++) 
				choose[i][j] = choose[i-1][j]+choose[i-1][j-1];
	}
	/*
	 * Computes combinations mod m.
	 */
	static void initChooseDP(int n, int m) {
		choose = new int[n+1][n+1];
		for (int i=0; i<=n; i++)
			choose[i][0] = 1;
		for (int i=1; i<=n; i++)
			for (int j=1; j<=i; j++) 
				choose[i][j] = (choose[i-1][j]+choose[i-1][j-1])%m;
	}
	/*
	 * Returns the nth Catalan number, C_n. 
	 * We also have C_n = (2n)Cn-(2n)C(n+1) = (2nCn)/(n+1).
	 */
	static int catalan(int n) {
		if (n == 0)
			return 1;
		return catalan(n-1)*(4*n+2)/(n+2);
	}
	// stand in cyclic function used for implementing Floyd's cycle finding algorithm.
	static int f(int n) {
		return (n+1)%10;
	}
	/*
	 * Floyd's cycle finding algorithm, adapted from the C++ code given in CP1 to Java.
	 * Returns an int[] whose first entry is the first position after x0 where there's a cycle,
	 * and whose second entry is the length of the cycle.
	 */
	static int[] floyd(int x0) { // The main phase of the algorithm, finding a repetition x_i = x_2i, hare speed is 2x tortoise’s
		int a = f(x0), b = f(f(x0)); // f(x0) is the element/node next to x0
		while (a != b) { a = f(a); b = f(f(b)); }
		// Find the position of mu, the hare and tortoise move at the same speeds
		int mu = 0; b = a; a = x0;
		while (a != b) {
			a = f(a); 
			b = f(b); 
			mu++; 
		}
		// Find the length of the shortest cycle starting from x_mu, hare moves, tortoise stays
		int lambda = 1; 
		b = f(a);
		while (a != b) { 
			b = f(b); 
			lambda++; 
		}
		return new int[]{mu,lambda};
	}
}
