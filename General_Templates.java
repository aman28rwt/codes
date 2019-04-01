import java.math.*;
import java.util.*;

public class General_Templates {
	// Modular inverse of a w.r.t. b
	long inv(long a, long b) {
		return a > 1 ? b - inv(b % a, a) * b / a : 1; // This will not work when (a * b) overflows
	}
	// Modular inverse of 2 w.r.t. to m : (m + 1) / 2

	long pow(long a, long x, long m) {
		long rv = 1;
		while (x != 0) {
			if ((x & 1) == 1) {
				rv = (rv * a) % m;
			}
			a = (a * a) % m;
			x = x >> 1;
		}
		return rv;
	}

	// This function gives (a * b) % m where a and b are 64 bits integers
	// and m = (1 << 61) - 1 ONLY
	long mul(long a, long b, long mod) {
		long l1 = (int) a, h1 = a >> 32;
		long l2 = (int) b, h2 = b >> 32;

		long l = l1 * l2, m = l1 * h2 + l2 * h1, h = h1 * h2;

		long rv = (l & mod) + (l >> 61) + (h << 3) + (m >> 29) + (m << 35 >> 3) + 1;
		rv = (rv & mod) + (rv >> 61);
		rv = (rv & mod) + (rv >> 61);

		return rv - 1;
	}

	// 1D-BIT
	void add(long[] bit, int ind, long val) {
		while (ind < bit.length) {
			bit[ind] += val;
			ind += ind & (-ind);
		}
	}

	long sum(long[] bit, int ind) {
		int rv = 0;
		while (ind > 0) {
			rv += bit[ind];
			ind -= ind & (-ind);
		}
		return rv;
	}
	// 1D-BIT

	// 2D-BIT
	// Copy above + below
	long[] b1, b2;

	void rangeUpdate(int l, int r, long val) {
		add(b1, l, val);
		add(b1, r + 1, -val);
		add(b2, l, val * (l - 1));
		add(b2, r + 1, -val * r);
	}

	long prefixSum(int ind) {
		return sum(b1, ind) * ind - sum(b2, ind);
	}

	long rangeSum(int l, int r) {
		return prefixSum(r) - prefixSum(l - 1);
	}
	// 2D-BIT

	// SegTree of Sum
	int[] tree, lazy;

	void updateSingle(int ind, int si, int li, int pos, int val) {
		if (si > li || si > pos || li < pos) {
			return;
		}

		if (si == li && si == pos) {
			tree[ind] = val;
			return;
		}

		int mid = (si + li) >> 1;

		updateSingle(2 * ind, si, mid, pos, val);
		updateSingle(2 * ind + 1, mid + 1, li, pos, val);

		tree[ind] = tree[2 * ind] + tree[2 * ind + 1];
	}

	void updateRange(int ind, int si, int li, int l, int r, int val) {
		if (lazy[ind] != 0) {
			tree[ind] = lazy[ind] * (li - si + 1);
			if (si != li) {
				lazy[2 * ind] = lazy[ind];
				lazy[2 * ind + 1] = lazy[ind];
			}
			lazy[ind] = 0;
		}

		if (si > li || si > r || li < l) {
			return;
		}

		if (si >= l && li <= r) {
			tree[ind] = val;
			if (si != li) {
				lazy[2 * ind] = val;
				lazy[2 * ind + 1] = val;
			}
			return;
		}

		int mid = (si + li) >> 1;

		updateRange(2 * ind, si, mid, l, r, val);
		updateRange(2 * ind + 1, mid + 1, li, l, r, val);

		tree[ind] = tree[2 * ind] + tree[2 * ind + 1];
	}

	int query(int ind, int si, int li, int l, int r) {
		if (lazy[ind] != 0) {
			tree[ind] = lazy[ind];
			if (si != li) {
				lazy[2 * ind] = lazy[ind];
				lazy[2 * ind + 1] = lazy[ind];
			}
			lazy[ind] = 0;
		}

		if (si > li || si > r || li < l) {
			return 0;
		}

		if (si >= l && li <= r) {
			return tree[ind];
		}

		int mid = (si + li) >> 1;

		int left = query(2 * ind, si, mid, l, r);
		int right = query(2 * ind + 1, mid + 1, li, l, r);

		return left + right;
	}
	// SegTree

	void Hashing() {
		String str = "";
		char[] arr = str.toCharArray();
		int n = arr.length;
		long[] f = new long[n], r = new long[n]; // f: forward hash | r: reverse hash

		long mod = BigInteger.probablePrime(30, new Random()).longValue();
		long mul = BigInteger.probablePrime(30, new Random()).longValue();

		long[] power = new long[n], invPow = new long[n];

		power[0] = 1;
		for (int i = 1; i < n; i++) {
			power[i] = (mul * power[i - 1]) % mod;
		}

		// Use this for hashing
		for (int i = 0; i < n; i++) {
			invPow[i] = pow(power[i], mod - 2, mod);
		}
		// Use this for hashing

		// this method is for inverse factorial
		// invPow[n - 1] = pow(power[n - 1], mod - 2, mod);
		// for (int i = n - 2; i >= 0; i--) {
		// invPow[i] = (invPow[i + 1] * (i + 1)) % mod;
		// }
		// this method is for inverse factorial

		f[0] = arr[0];
		for (int i = 1; i < n; i++) {
			f[i] = ((arr[i] * power[i]) % mod + f[i - 1]) % mod;
		}

		r[n - 1] = arr[n - 1];
		for (int i = n - 2; i >= 0; i--) {
			r[i] = ((arr[i] * power[n - 1 - i]) % mod + r[i + 1]) % mod;
		}

		for (int i = 0, j = str.length() - 1; j < n; i++, j++) {
			long k1 = fh(f, i, j, invPow, mod); // hash of forward string
			long k2 = rh(r, i, j, invPow, mod, n); // hash of reverse string
			System.out.println(k1 + k2);
		}
	}

	long fh(long[] f, int i, int j, long[] invPow, Long M) {
		long jh = f[j]; // hash for S[0..j]
		long ih = (i > 0) ? f[i - 1] : 0; // hash for S[0..i-1]
		long sb = ((jh + M - ih) * invPow[i]) % M;
		return sb;
	}

	long rh(long[] r, int i, int j, long[] invPow, Long M, int n) {
		long ih = r[i]; // hash for reverse S[i..n-1]
		long jh = (j < n - 1) ? r[j + 1] : 0; // hash for reverse S[j+1..n-1]
		long sb = ((ih + M - jh) * invPow[n - 1 - j]) % M;
		return sb;
	}

	long[] radixSort(long[] f) {
		int n = f.length;
		long[] to = new long[n];
		{
			int[] b = new int[65537];
			for (int i = 0; i < n; i++)
				b[1 + (int) (f[i] & 0xffff)]++;
			for (int i = 1; i <= 65536; i++)
				b[i] += b[i - 1];
			for (int i = 0; i < n; i++)
				to[b[(int) (f[i] & 0xffff)]++] = f[i];
			long[] d = f;
			f = to;
			to = d;
		}
		{
			int[] b = new int[65537];
			for (int i = 0; i < n; i++)
				b[1 + (int) (f[i] >>> 16 & 0xffff)]++;
			for (int i = 1; i <= 65536; i++)
				b[i] += b[i - 1];
			for (int i = 0; i < n; i++)
				to[b[(int) (f[i] >>> 16 & 0xffff)]++] = f[i];
			long[] d = f;
			f = to;
			to = d;
		}
		{
			int[] b = new int[65537];
			for (int i = 0; i < n; i++)
				b[1 + (int) (f[i] >>> 32 & 0xffff)]++;
			for (int i = 1; i <= 65536; i++)
				b[i] += b[i - 1];
			for (int i = 0; i < n; i++)
				to[b[(int) (f[i] >>> 32 & 0xffff)]++] = f[i];
			long[] d = f;
			f = to;
			to = d;
		}
		{
			int[] b = new int[65537];
			for (int i = 0; i < n; i++)
				b[1 + (int) (f[i] >>> 48 & 0xffff)]++;
			for (int i = 1; i <= 65536; i++)
				b[i] += b[i - 1];
			for (int i = 0; i < n; i++)
				to[b[(int) (f[i] >>> 48 & 0xffff)]++] = f[i];
			long[] d = f;
			f = to;
			to = d;
		}
		return f;
	}

	class DisJointSet {
		int[] table, rank, count;
		int size;

		DisJointSet(int size) {
			this.table = new int[size];
			this.rank = new int[size];
			this.count = new int[size];
			this.size = size;
			for (int i = 0; i < size; i++) {
				this.table[i] = i;
				this.rank[i] = 1;
				this.count[i] = 1;
			}
		}

		boolean isSame(int x, int y) {
			return root(x) == root(y);
		}

		int root(int node) {
			if (table[node] == node) {
				return node;
			} else {
				return table[node] = root(table[node]);
			}
		}

		void union(int x, int y) {
			x = root(x);
			y = root(y);
			if (x != y)
				this.size--;
			if (rank[x] < rank[y]) {
				table[x] = y;
				count[y] += count[x];
			} else if (rank[x] > rank[y]) {
				table[y] = x;
				count[x] += count[y];
			} else if (x != y) {
				table[y] = x;
				count[x] += count[y];
				rank[x]++;
			}
		}
	}

	// Sqrt
	long sqrt(long n) {
		long i = Math.max(0, (long) Math.sqrt(n) - 2);
		while (i * i <= n)
			i++;
		return i - 1;
	}

	// Mobius O(N) : Multiplicative fn
	// mob(p^k) = [k == 0] - [k == 1]
	int[] mobius(int n) {
		boolean[] composites = new boolean[n + 1];
		int[] mob = new int[n + 1], primes = new int[n + 1];
		mob[1] = 1;
		for (int i = 2, k = 0; i <= n; i++) {
			if (!composites[i]) {
				primes[k++] = i;
				mob[i] = -1;
			}
			for (int j = 2; j < k && i * 1.0 * primes[j] <= n; j++) {
				composites[i * primes[j]] = true;
				if (i % primes[j] != 0) {
					mob[i * primes[j]] = mob[i] * mob[primes[j]];
				} else {
					mob[i * primes[j]] = 0;
					break;
				}
			}
		}
		return mob;
	}

	// Euler Totient's O(N) : Multiplicative fn
	// phi(p^k) = p^k - p^(k - 1)
	int[] eulerphi(int n) {
		boolean[] composites = new boolean[n + 1];
		int[] phi = new int[n + 1], primes = new int[n + 1];
		phi[1] = 1;
		for (int i = 2, k = 0; i <= n; i++) {
			if (!composites[i]) {
				primes[k++] = i;
				phi[i] = i - 1;
			}
			for (int j = 0; j < k && i * 1.0 * primes[j] <= n; j++) {
				composites[i * primes[j]] = true;
				if (i % primes[j] != 0) {
					phi[i * primes[j]] = phi[i] * phi[primes[j]];
				} else {
					phi[i * primes[j]] = phi[i] * primes[j];
					break;
				}
			}
		}
		return phi;
	}

	int[] eulerPhilgN(int n) {
		int[] arr = new int[n + 1];
		arr[1] = 1;
		for (int i = 1; i <= n; i++) {
			if (arr[i] == 0) {
				continue;
			}
			for (int j = 2 * i; j <= n && j > 0; j += i) {
				arr[j] -= arr[i];
			}
		}
		return arr;
	}

	// Linear Sieve: 1e8 in 1000ms
	int[] linearSieve(int n) {
		boolean[] composites = new boolean[n + 1];
		int[] prime = new int[n + 1];
		int k = 0;
		for (int i = 2; i <= n; i++) {
			if (!composites[i]) {
				prime[k++] = i;
			}
			for (int j = 0; j < k && i * 1.0 * prime[j] <= n; j++) {
				composites[i * prime[j]] = true;
				if (i % prime[j] == 0) {
					break;
				}
			}
		}
		return Arrays.copyOf(prime, k);
	}

	// Fastest Sieve by uwi: 1e8 in 600ms
	int[] uwiSieve(int n) {
		if (n <= 32) {
			int[] primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };
			for (int i = 0; i < primes.length; i++) {
				if (n < primes[i]) {
					return Arrays.copyOf(primes, i);
				}
			}
			return primes;
		}

		int u = n + 32;
		double lu = Math.log(u);
		int[] ret = new int[(int) (u / lu + u / lu / lu * 1.5)];
		ret[0] = 2;
		int pos = 1;

		int[] isp = new int[(n + 1) / 32 / 2 + 1];
		int sup = (n + 1) / 32 / 2 + 1;

		int[] tprimes = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };
		for (int tp : tprimes) {
			ret[pos++] = tp;
			int[] ptn = new int[tp];
			for (int i = (tp - 3) / 2; i < tp << 5; i += tp)
				ptn[i >> 5] |= 1 << (i & 31);
			for (int i = 0; i < tp; i++) {
				for (int j = i; j < sup; j += tp)
					isp[j] |= ptn[i];
			}
		}

		// 3,5,7
		// 2x+3=n
		int[] magic = { 0, 1, 23, 2, 29, 24, 19, 3, 30, 27, 25, 11, 20, 8, 4, 13, 31, 22, 28, 18, 26, 10, 7, 12, 21, 17,
				9, 6, 16, 5, 15, 14 };
		int h = n / 2;
		for (int i = 0; i < sup; i++) {
			for (int j = ~isp[i]; j != 0; j &= j - 1) {
				int pp = i << 5 | magic[(j & -j) * 0x076be629 >>> 27];
				int p = 2 * pp + 3;
				if (p > n)
					break;
				ret[pos++] = p;
				for (int q = pp; q <= h; q += p)
					isp[q >> 5] |= 1 << (q & 31);
			}
		}

		return Arrays.copyOf(ret, pos);
	}

	// Matrix Multiplication : Stressen
	long[][] multiply(long[][] a, long[][] b, long mod) {
		if (a[0].length != b.length) {
			return null;
		}

		long[][] rv = new long[a.length][b[0].length];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b[0].length; j++) {
				for (int k = 0; k < b.length; k++) {
					rv[i][j] += (a[i][k] * b[k][j]) % mod;
					rv[i][j] %= mod;
				}
			}
		}

		return rv;
	}

	// LCA
	int lca(int u, int v, int[] lvl, int[][] ancestor, int log) { // Here log is lg(nNoOfNodes)
		if (lvl[u] > lvl[v]) {
			u = u ^ v;
			v = u ^ v;
			u = u ^ v;
		}
		for (int i = log - 1; i >= 0; i--) {
			if (lvl[ancestor[v][i]] >= lvl[u]) {
				v = ancestor[v][i];
			}
		}
		if (u == v) {
			return u;
		}
		for (int i = log - 1; i >= 0; i--) {
			if (ancestor[u][i] != ancestor[v][i]) {
				u = ancestor[u][i];
				v = ancestor[v][i];
			}
		}
		return ancestor[u][0];
	}

	void dfs_lca(int[][] g, int u, int p, int l, int[][] ancestor, int[] lvl, int log) {
		ancestor[u][0] = p;
		lvl[u] = l;
		for (int i = 1; i < log; i++) {
			ancestor[u][i] = ancestor[ancestor[u][i - 1]][i - 1];
		}
		for (int v : g[u]) {
			if (v != p) {
				dfs_lca(g, v, u, l + 1, ancestor, lvl, log);
			}
		}
	}
	// LCA

	// FFT without MOD
	long[] convolute(long[] a, long[] b, int p) {
		int m = Integer.highestOneBit(Math.max(Math.max(a.length, b.length) - 1, 1)) << 2;
		double[][] fa = fft(a, m, false);
		double[][] fb = a == b ? fa : fft(b, m, false);
		for (int i = 0; i < m; i++) {
			double nfa0 = fa[0][i] * fb[0][i] - fa[1][i] * fb[1][i];
			double nfa1 = fa[0][i] * fb[1][i] + fa[1][i] * fb[0][i];
			fa[0][i] = nfa0;
			fa[1][i] = nfa1;
		}
		fft(fa[0], fa[1], true);
		long[] ret = new long[m];
		for (int i = 0; i < m; i++) {
			ret[i] = Math.round(fa[0][i]);
		}
		return ret;
	}

	double[][] fft(long[] srcRe, int n, boolean inverse) {
		int m = srcRe.length;
		double[] dstRe = new double[n];
		double[] dstIm = new double[n];
		for (int i = 0; i < m; i++) {
			dstRe[i] = srcRe[i];
		}

		int h = Integer.numberOfTrailingZeros(n);
		for (int i = 0; i < n; i++) {
			int rev = Integer.reverse(i) >>> 32 - h;
			if (i < rev) {
				double d = dstRe[i];
				dstRe[i] = dstRe[rev];
				dstRe[rev] = d;
			}
		}

		for (int s = 2; s <= n; s <<= 1) {
			int nt = s >>> 1;
			double theta = inverse ? -2 * Math.PI / s : 2 * Math.PI / s;
			double wRe = Math.cos(theta);
			double wIm = Math.sin(theta);
			for (int j = 0; j < n; j += s) {
				double wr = 1, wi = 0;
				for (int t = j; t < j + nt; t++) {
					int jp = t + nt;
					double re = dstRe[jp] * wr - dstIm[jp] * wi;
					double im = dstRe[jp] * wi + dstIm[jp] * wr;
					dstRe[jp] = dstRe[t] - re;
					dstIm[jp] = dstIm[t] - im;
					dstRe[t] += re;
					dstIm[t] += im;
					double nwre = wr * wRe - wi * wIm;
					double nwim = wr * wIm + wi * wRe;
					wr = nwre;
					wi = nwim;
				}
			}
		}

		if (inverse) {
			for (int i = 0; i < n; i++) {
				dstRe[i] /= n;
				dstIm[i] /= n;
			}
		}

		return new double[][] { dstRe, dstIm };
	}

	void fft(double[] re, double[] im, boolean inverse) {
		int n = re.length;
		int h = Integer.numberOfTrailingZeros(n);
		for (int i = 0; i < n; i++) {
			int rev = Integer.reverse(i) >>> 32 - h;
			if (i < rev) {
				double d = re[i];
				re[i] = re[rev];
				re[rev] = d;
				d = im[i];
				im[i] = im[rev];
				im[rev] = d;
			}
		}

		for (int s = 2; s <= n; s <<= 1) {
			int nt = s >>> 1;
			double theta = inverse ? -2 * Math.PI / s : 2 * Math.PI / s;
			double wRe = Math.cos(theta);
			double wIm = Math.sin(theta);
			for (int j = 0; j < n; j += s) {
				double wr = 1, wi = 0;
				for (int t = j; t < j + nt; t++) {
					int jp = t + nt;
					double lre = re[jp] * wr - im[jp] * wi;
					double lim = re[jp] * wi + im[jp] * wr;
					re[jp] = re[t] - lre;
					im[jp] = im[t] - lim;
					re[t] += lre;
					im[t] += lim;
					double nwre = wr * wRe - wi * wIm;
					double nwim = wr * wIm + wi * wRe;
					wr = nwre;
					wi = nwim;
				}
			}
		}

		if (inverse) {
			for (int i = 0; i < n; i++) {
				re[i] /= n;
				im[i] /= n;
			}
		}
	}
	// FFT without MOD

	// KMP : O(n + m)
	int[] build_KMP_Table(char[] small) {
		int n = small.length;
		int[] mp = new int[n + 1];
		mp[0] = -1;
		for (int i = 1, j = 0; i < n; i++) {
			while (j >= 0 && small[i] != small[j]) {
				j = mp[j];
			}
			mp[i + 1] = ++j;
		}
		return mp;
	}

	int findFirstInd(char[] big, char[] small, int[] KMP) {
		for (int i = 0, j = 0; i < big.length; i++) {
			while (j >= 0 && big[i] != small[j]) {
				j = KMP[j];
			}
			if (++j == small.length) {
				return i - small.length + 1;
			}
		}
		return -1; // pattern not found
	}
	// KMP : O(n + m)
	
	// nextPermutation
	boolean nextPermutation(int[] nums) {
        // Find the last location where an element is smaller than its right neighbour.
        int left = -1, first = nums[0];
        for (int i = nums.length - 1; i > 0; i--) {
            if (nums[i] > nums[i - 1]) {
                left = i - 1;
                break;
            }
        }
        
        // If no element smaller than right neighbour was found, sort entire array and exit.
        if (left == -1) {
            Arrays.sort(nums);
            return nums[0] >= first;
        }
        
        // Find the *smallest* element to the right of "left" that is larger than "left"
        int right = -1;
        for (int i = left + 1; i < nums.length; i++) {
            int n = nums[i];
            int l = nums[left];
            if ((n > l && right == -1) || n > l && n < nums[right]){
                right = i;
            } 
        }
        
        // Swap the values of left and right
        int temp = nums[left];
        nums[left] = nums[right];
        nums[right] = temp;
        
        // Sort all elements to the right of lefts index
        Arrays.sort(nums, left + 1, nums.length);
        
        return nums[0] >= first;
    }
}