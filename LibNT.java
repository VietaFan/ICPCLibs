import java.util.*;
import java.io.*;
public class LibNT {
	public static int smallPrimeMax = 10;
	public static int[] smallPrimes = {2,3,5,7};
	
	public static void initPrimes(int n) {
		smallPrimeMax=n;
		ArrayList L = getPrimesTo(n);
		smallPrimes = new int[L.size()];
		for (int i=0; i<smallPrimes.length; i++)
			smallPrimes[i] = (int)L.get(i);
	}
	public static boolean isPrime(long n) {
		if (n < smallPrimeMax) {
			return Arrays.binarySearch(smallPrimes,(int)n) >= 0;
		} else {
			HashMap<Long,Integer> L = primeFact(n);
			return L.size() == 1 && L.containsValue(1);
		}
	}
	/*
	 * Returns an ArrayList containing all prime numbers less than or equal to n.
	 * Uses the sieve of Eratosthenes.
	 */
	public static ArrayList<Integer> getPrimesTo(int n) {
		boolean[] isprime = new boolean[n+1];
		Arrays.fill(isprime, true);
		ArrayList<Integer> primes = new ArrayList<Integer>();
		primes.add(2);
		int k=0,p,pos=2;
		while (pos*pos<=n) {
			p = primes.get(k++);
			while (++pos<2*p)
				if (isprime[pos])
					primes.add(pos);
			for (int j=2*p; j<=n; j+=p)
				isprime[j]=false;
		}
		for (;k<primes.size();k++) {
			p = primes.get(k);
			for (int j=2*p; j<=n; j+=p)
				isprime[j] = false;
		}
		for (;pos<=n;pos++)
			if (isprime[pos])
				primes.add(pos);
		return primes;
	}
	/*
	 * Returns a HashMap whose keys are the prime factors of n and whose values 
	 * are the exponents in the prime factorization of n. Uses naive algorithm, starting
	 * with precomputed primes and then using 1 and 5 mod 6.
	 */
	public static HashMap<Long, Integer> primeFact(long n) {
		int ct;
		long q;
		HashMap<Long, Integer> ans = new HashMap<Long, Integer>();
		for (long p: smallPrimes) {
			ct = 0;
			while (n%p == 0) {
				n /= p;
				ct++;
			}
			if (ct > 0)
				ans.put(p, ct);
		}
		for (long p=(smallPrimeMax/6)*6; p*p<=n; p+=6) {
			q = p+1;
			ct = 0;
			while (n%q == 0) {
				n /= q;
				ct++;
			}
			if (ct > 0)
				ans.put(q, ct);
			q = p+5;
			ct = 0;
			while (n%q == 0) {
				n /= q;
				ct++;
			}
			if (ct > 0)
				ans.put(q, ct);
		}
		if (n > 1) {
			ans.put(n,1);
		}
		return ans;
	}
	public static long product(HashMap<Long, Integer> factors) {
		long ans = 1;
		int e;
		for (long f: factors.keySet()) {
			e = factors.get(f);
			for (int i=0; i<e; i++)
				ans *= e;
		}
		return ans;
	}
	/*
	 * Returns the smallest j such that L[j] >= x.
	 */
	public static int lower_bound(int[] L, int x) {
		int a = 0, b = L.length-1, c;
		while (a < b) {
			c = (a+b)>>1;
			if (L[c] < x) {
				a = c+1;
			} else {
				b = c;
			}
		}
		return a;
	}
	public static int lower_bound(long[] L, long x) {
		int a = 0, b = L.length-1, c;
		while (a < b) {
			c = (a+b)>>1;
			if (L[c] < x) {
				a = c+1;
			} else {
				b = c;
			}
		}
		return a;
	}
	/*
	 * Returns the smallest j such that L[j] > x.
	 */
	public static int upper_bound(int[] L, int x) {
		int a = 0, b = L.length-1, c;
		while (a < b) {
			c = (a+b)>>1;
			if (L[c] <= x) {
				a = c+1;
			} else {
				b = c;
			}
		}
		return a;
	}
	/*
	 * gcd function from CP1
	 */
	static int gcd(int a, int b) {
		return (b == 0 ? a : gcd(b, a % b)); 
	}
	static long gcd(long a, long b) {
		return (b == 0 ? a : gcd(b, a % b)); 
	}
	/*
	 * lcm function from CP1
	 */
	static int lcm(int a, int b) { 
		return (a * (b / gcd(a, b)));
	} // divide before multiply!
	/*
	 * Returns the number of positive integers k such that k < n and gcd(k,n) = 1.
	 */
	static long phi(long n) {
		return phi(primeFact(n));
	}
	static int phi(int n) {
		return (int)phi((long)n);
	}
	/*
	 * Returns the number of divisors of n.
	 */
	static int tau(long n) {
		return tau(primeFact(n));
	}
	/*
	 * Returns the sum of the divisors of n.
	 */
	static long sigma(long n) {
		return sigma(primeFact(n));
	}
	static int sigma(int n) {
		return (int)sigma((long)n);
	}
	static long phi(HashMap<Long, Integer> factors) {
		long ans = 1;
		int e;
		for (long p: factors.keySet()) {
			e = factors.get(p);
			for (int i=1; i<e; i++) {
				ans *= p;
			}
			ans *= (p-1);
		}
		return ans;
	}
	static int tau(HashMap<Long, Integer> factors) {
		int ans = 1;
		for (int e: factors.values()) {
			ans *= (e+1);
		}
		return ans;
	}
	static long sigma(HashMap<Long, Integer> factors) {
		long ans = 1,d;
		int e;
		for (long p: factors.keySet()) {
			d=0;
			e = factors.get(p);
			for (int i=0; i<=e; i++)
				d *= p;
			d--;
			d /= (p-1);
			ans *= d;
		}
		return ans;
	}
	/*
	 * Returns the multiplicative inverse of a modulo m assuming gcd(a,m) = 1.
	 * If gcd(a,m) != 1, the output is undefined.
	 */
	static long inverse(long a, long m) {
		if (a == 0) return 1;
		if (a > m)
			return inverse(a%m, m);
		return ((1-m*inverse(m%a, a))/a+m)%m;
	}
	/*
	 * Given coprime integers a,b, this returns an integer y such that
	 * if x is the inverse of a mod b, then ax+by = 1. 
	 */
	static long secondCoeff(long a, long b) {
		return (1-a*inverse(a,b))/b;
	}
	/*
	 * Returns b^p mod m. Use (BigInt).modPow(BigInt e, BigInt m) if you need larger integers.
	 */
	static int powerMod(int b, int p, int m) {
		long pow2mods[] = new long[31];
		pow2mods[0] = b;
		for (int i=1; i<31; i++) {
			pow2mods[i] = pow2mods[i-1]*pow2mods[i-1];
			pow2mods[i] %= m;
		}
		long ans = 1;
		for (int i=0; i<31; i++) {
			if ((p&(1<<i))>0) {
				ans *= pow2mods[i];
				ans %= m;
			}
		}
		return (int)ans;
	}
}
