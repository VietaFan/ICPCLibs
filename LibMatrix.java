import java.io.*;
import java.util.*;
class Matrix {
	public double[][] grid;
	public double[][] basis;
	public boolean[] isPivot;
	public int[] pivotCol;
	public int[] freeNum;
	static double EPS = 1e-6;
	public int m, n, nrows, ncols;
	public Matrix L, U, P;
	public int[] perm;
	static boolean rounding = true;
	public boolean factored;
	public double det;
	public Matrix(int n) {
		this(n,n);
	}
	public Matrix(int m, int n) {
		this(new double[m][n]);
	}
	public Matrix(double[][] A) {
		factored = false;
		nrows = m = A.length;
		ncols = n = A[0].length;
		grid = new double[m][n];
		for (int i=0; i<m; i++) 
			for (int j=0; j<n; j++)
				grid[i][j] = A[i][j];
	}
	public Matrix(Matrix m) {
		this(m.grid);
	}
	public Matrix add(Matrix M) {
		Matrix ans = new Matrix(this);
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				ans.grid[i][j] += M.grid[i][j];
		return ans;
	}
	public Matrix scale(double c) {
		Matrix ans = new Matrix(this);
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				ans.grid[i][j] *= c;
		return ans;
	}
	public Matrix sub(Matrix M) {
		Matrix ans = new Matrix(this);
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				ans.grid[i][j] -= M.grid[i][j];
		return ans;
	}
	public Matrix mul(Matrix M) {
		if (n != M.m) {
			System.err.printf("Attempted matrix multiplication with invalid sizes: %d X %d and %d X %d.\n", m, n, M.m, M.n);
			return null;
		}
		Matrix ans = new Matrix(m, M.n);
		for (int i=0; i<m; i++) 
			for (int j=0; j<M.n; j++) 
				for (int k=0; k<n; k++)
					ans.grid[i][j] += this.grid[i][k]*M.grid[k][j];
		if (rounding)
			ans.intRound();
		return ans;			
	}
	public static Matrix Identity(int n) {
		double[][] L = new double[n][n];
		for (int i=0; i<n; i++)
			L[i][i] = 1.0;
		return new Matrix(L);
	}
	public static Matrix Zero(int m, int n) {
		return new Matrix(new double[m][n]);
	}
	public static Matrix Zero(int n) {
		return new Matrix(new double[n][n]);
	}
	public Matrix transpose() {
		double[][] L = new double[n][m];
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				L[j][i] = grid[i][j];
		return new Matrix(L);
	}
	public void rowSwap(int i, int j) {
		double[] temp = grid[i];
		grid[i] = grid[j];
		grid[j] = temp;
	}
	public void rowScale(int row, double c) {
		for (int i=0; i<n; i++)
			grid[row][i] *= c;
	}
	public void rowAdd(int sourceRow, int targetRow, double c) {
		for (int i=0; i<n; i++)
			grid[targetRow][i] += grid[sourceRow][i]*c;
	}
	public Matrix echelonForm() {
		if (!factored)
			this.LUfactorize();
		return new Matrix(U);
	}
	public void LUfactorize() {
		U = new Matrix(this);
		P = Matrix.Identity(m);
		int r = 0, temp;
		// find permutation matrix P
		for (int c=0; r<m && c<n; c++) {
			if (Math.abs(U.grid[r][c]) < EPS)  {
				for (int i=r+1; i<m; i++) {
					if (Math.abs(U.grid[i][c]) < EPS) {
						continue;
					}
					U.rowSwap(r,i);
					P.rowSwap(r,i);
					break;
				}
				if (Math.abs(U.grid[r][c]) < EPS)
					continue;
			}
			for (int i=r+1; i<m; i++) {
				U.rowAdd(r, i, -U.grid[i][c]/U.grid[r][c]);
			}
			r++;
		}
		// find L and U
		U = P.mul(this);
		L = Matrix.Identity(m);
		r = 0;
		for (int c=0; r<m && c<n; c++) {
			if (Math.abs(U.grid[r][c]) < EPS)
				continue;
			for (int i=r+1; i<m; i++) {
				L.grid[i][r] = U.grid[i][c]/U.grid[r][c];
				U.rowAdd(r, i, -U.grid[i][c]/U.grid[r][c]);
			}
			r++;
		}	
		perm = new int[m];
		for (int i=0; i<m; i++) {
			for (int j=0; j<m; j++)
				if (P.grid[i][j] > 0.5) // if there's a 1
					perm[i] = j;
		}
		if (rounding) {
			L.intRound();
			U.intRound();
		}
		//compute determinant (it's cheap like O(N) so might as well)
		if (m == n) {
			det = 1.0;
			for (int i=0; i<m; i++)
				det *= grid[i][i];
		}
	}
	public void findBasis() {
		// find the basis for the null space of this matrix	- useful for generating additional solutions			
		int row = 0, col = 0, oldcol = -1;
		isPivot = new boolean[n];
		freeNum = new int[n];
		pivotCol = new int[m];
		while (row < m && col < n) {
			if (Math.abs(U.grid[row][col]) > EPS) {
				if (col != oldcol) {
					isPivot[col] = true;
					pivotCol[row] = col;
				} else {
					pivotCol[row] = -1;
				}
				row++;
			} else {
				col++;
			}
		}
		for (; row<m; row++)
			pivotCol[row] = -1;
		int nfree=0;
		for (int i=0; i<n; i++)
			if (!isPivot[i])
				freeNum[i] = nfree++;
		basis = new double[nfree][n];
		for (int i=0; i<n; i++) {
			if (!isPivot[i]) {
				for (int j=0; j<nfree; j++)
					basis[j][i] = 0;
				basis[freeNum[i]][i] = 1;
			}
		}
		for (int i=m-1; i>-1; i--) {
			if (pivotCol[i] == -1)
				continue;
			for (int j=0; j<nfree; j++) {
				basis[j][pivotCol[i]] = 0;
				for (int k=pivotCol[i]+1; k<n; k++) 
					basis[j][pivotCol[i]] -= basis[j][k]*U.grid[i][k];
				basis[j][pivotCol[i]] /= U.grid[i][pivotCol[i]];
			}
		}
		factored = true;		
	}
	public double determinant() {
		if (n != m)
			return Double.NaN;
		if (!factored)
			this.LUfactorize();
		return det;
	}
	// returns null if not invertible
	public Matrix inverse() {
		if (n != m || Math.abs(det) < EPS)
			return null;
		Matrix ans = new Matrix(this);
		Matrix inv = Matrix.Identity(n);
		int r = 0,j;
		for (int c=0; r<m && c<n; c++) {
			if (Math.abs(ans.grid[r][c]) < EPS)  {
				for (int i=r+1; i<m; i++) {
					if (Math.abs(ans.grid[i][c]) < EPS) {
						continue;
					}
					inv.rowSwap(r,i);
					ans.rowSwap(r,i);
					break;
				}
				if (Math.abs(ans.grid[r][c]) < EPS)
					continue;
			}
			inv.rowScale(r, 1.0/ans.grid[r][c]);
			ans.rowScale(r, 1.0/ans.grid[r][c]);
			for (int i=r+1; i<m; i++) {
				inv.rowAdd(r, i, -ans.grid[i][c]);
				ans.rowAdd(r, i, -ans.grid[i][c]);
			}
			r++;
		}
		for (int i=n-1; i>-1; i--) {
			for (j=m-1; j>-1 && Math.abs(ans.grid[j][i]) < EPS; j--);
			for (int k=0; k<j; k++) {
				inv.rowAdd(j, k, -ans.grid[k][i]/ans.grid[j][i]);
				ans.rowAdd(j, k, -ans.grid[k][i]/ans.grid[j][i]);
			}
		}
		return inv;
	}
	// returns the reduced echelon form of the matrix
	public Matrix getREF() {
		Matrix M = this.echelonForm();
		int j;
		for (int i=n-1; i>-1; i--) {
			for (j=m-1; j>-1 && Math.abs(M.grid[j][i]) < EPS; j--);
			for (int k=0; k<j; k++)
				M.rowAdd(j, k, -M.grid[k][i]/M.grid[j][i]);
		}
		if (rounding)
			M.intRound();
		return M;
	}
	// rounds entries off to integers - especially useful for making sure zeroes are actually zero
	public void intRound() {
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				if (Math.abs(grid[i][j]-Math.round(grid[i][j])) < EPS)
					grid[i][j] = Math.round(grid[i][j]);
	}
	// returns this * x
	public double[] vecMul(double[] x) {
		double[] b = new double[m];
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++)
				b[i] += grid[i][j]*x[j];
		}
		return b;
	}
	//assumes augmented matrix in echelon form, tells whether the associated system is consistent
	public boolean augConsistent() {
		for (int i=m-1; i>-1; i--) {
			for (int j=0; j<n-1; j++) {
				if (Math.abs(grid[i][j]) > EPS)
					return true;
			}
			if (Math.abs(grid[i][n-1]) > EPS)
				return false;
		}
		return true;
	}
	public Matrix submatrix(int rowstart, int colstart, int numrows, int numcols) {
		Matrix M = new Matrix(numrows, numcols);
		for (int i=rowstart; i<rowstart+numrows; i++)
			for (int j=colstart; j<colstart+numcols; j++)
				M.grid[i-rowstart][j-colstart] = grid[i][j];
		return M;
	}
	public Matrix cofactor(int row, int col) {
		double[][] L = new double[m-1][n-1];
		int a=0,b;
		for (int i=0; i<m; i++) {
			if (i == row) continue;
			b=0;
			for (int j=0; j<n; j++)
				if (j != col)
					L[a][b++] = grid[i][j];
			a++;
		}
		return new Matrix(L);
	}
	// assumes this is an augmented matrix in echelon form
	public int augNumPivots() {
		int r = 0;
		for (int c=0; r<m && c<n-1; c++) {
			if (Math.abs(grid[r][c]) > EPS)
				r++;
		}
		return r;
	}
	public String toString() {
		String s = "----------\n";
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++)
				s += String.format("%f\t", grid[i][j]);
			s += '\n';
		}
		s += "----------\n";
		return s;
	}
	public void printVec(double[] vec) {
		System.out.print(new Matrix(new double[][]{vec}));
	}
	// returns a particular solution to the system this*x=b0. use this.basis to get general solution
	public double[] solve(double[] b0) {
		if (!factored)
			LUfactorize();
		if (basis == null)
			findBasis();
		double[] b = new double[m];
		for (int i=0; i<m; i++)
			b[i] = b0[perm[i]];
		double[] x = new double[n];
		double[] y = new double[m];
		double[] sol = new double[n];
		for (int i=0; i<m; i++) {
			y[i] = b[i];
			for (int j=0; j<i; j++) 
				y[i] -= L.grid[i][j]*y[j];
			y[i] /= L.grid[i][i];
		}
		for (int i=0; i<n; i++)
			if (!isPivot[i])
				sol[i] = 0;
		double test;
		for (int i=m-1; i>-1; i--) {
			// b[i] better be equal to the correct linear combination or it's inconsistent
			if (pivotCol[i] == -1) {
				continue;
			}
			sol[pivotCol[i]] = y[i];
			for (int j=pivotCol[i]+1; j<n; j++) 
				sol[pivotCol[i]] -= sol[j]*U.grid[i][j];
			sol[pivotCol[i]] /= U.grid[i][pivotCol[i]];
		}
		// if this doesn't work, it's inconsistent
		for (int i=0; i<m; i++) {
			test = 0.0;
			for (int j=0; j<n; j++)
				test += grid[i][j]*sol[j];
			if (Math.abs(test-b0[i]) > EPS) {
				//System.out.printf("%f != %f\n",test,b0[i]);
				//printVec(sol);
				return null;
			}
		}
		return sol;
	}
	// says whether two matrices are equal and gives debugging print
	public boolean equals(Matrix other) {
		if (m != other.m || n != other.n)
			return false;
		for (int i=0; i<other.m; i++)
			for (int j=0; j<other.n; j++)
				if (Math.abs(grid[i][j]-other.grid[i][j]) > EPS) {
					System.out.printf("Unequal matrices: entries in row %d and column %d differ by %.12f.\n",i,j,grid[i][j]-other.grid[i][j]);
					return false;
				}
		return true;
	}
}

public class LibMatrix {
	// false = fail, true = it worked
	static boolean testSolver(Matrix A, double[] b) {
		A.LUfactorize();
		if (!A.L.mul(A.U).equals(A.P.mul(A))) {
			System.out.println("bad factorization");
			return false;
		}
		double[] sol = A.solve(b);
		if (sol == null) {
			// make sure it actually is inconsistent
			double[][] L = new double[A.m][A.n+1];
			for (int i=0; i<A.m; i++) {
				for (int j=0; j<A.n; j++)
					L[i][j] = A.grid[i][j];
				L[i][A.n] = b[i];
			}
			Matrix aug = new Matrix(L);
			aug = aug.getREF();
			if (aug.augConsistent()) {
				System.out.println("fail: it's consistent but no solution found.");
				return false;
			}
			return true;
		}	
		// make sure the solution actually works
		double[] x = new double[A.n];
		double[] c = new double[A.basis.length];
		double[] bprime = new double[A.m];
		for (int i=0; i<16; i++) {
			for (int j=0; j<c.length; j++)
				c[j] = (Math.random()-0.5)*200.0;
			for (int j=0; j<A.n; j++) {
				x[j] = sol[j];
				for (int t=0; t<c.length; t++)
					x[j] += c[t] * A.basis[t][j];
			}
			bprime = A.vecMul(x);
			for (int j=0; j<A.m; j++) {
				if (Math.abs(b[j]-bprime[j]) > 1e-4) {
					System.out.println("fail:");
					System.out.println(new Matrix(new double[][] {bprime}));
					System.out.println(new Matrix(new double[][] {b}));
					System.out.println(new Matrix(new double[][] {sol}));
					return false;
				}
			}
		}
		return true;
	}
	public static void main(String[] args) {	
		System.out.println(testSolver(new Matrix(new double[][] { new double[] {2,3}, new double[] {1,0}, new double[] {0,1}}), new double[] {3,0,0}));
		double[][] L;
		double[] b;
		Random r = new Random();
		for (int i=0; i<100000; i++) {
			L = new double[r.nextInt(6)+1][r.nextInt(6)+1];
			System.out.printf("Testing matrix with dimensions %d x %d ...", L.length, L[0].length);
			for (int a=0; a<L.length; a++)
				for (int c=0; c<L[0].length; c++) {
					//L[a][c] = r.nextInt(3)-1;
					L[a][c] = (Math.random()-0.5)*200.0;
				}
			Matrix M = new Matrix(L);
			b = new double[M.m];
			for (int a=0; a<M.m; a++)
				b[a] = (Math.random()-0.5)*200.0;
				//b[a] = r.nextInt(3)-1;
			if (testSolver(M,b))
				System.out.printf("test %d ok!\n",i+1);
			else {
				System.out.printf("failed!\n");
				System.out.print(M);
				System.out.print(new Matrix(new double[][] {b}));
				break;
			}
		}
		
	}
}