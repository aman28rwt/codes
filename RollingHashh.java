import java.math.BigInteger;
import java.util.Random;

public class RollingHashh {
	public class RollingHash {
		public RollingHashFactory rhf;
		public long[][] buf;
		public int p;

		public RollingHash(int bufsize, RollingHashFactory rhf) {
			buf = new long[rhf.deg][bufsize + 1];
			this.rhf = rhf;
			this.p = 1;
		}

		public void add(int c) {
			for (int i = 0; i < rhf.deg; i++)
				buf[i][p] = (buf[i][p - 1] * rhf.muls[i] + c) % rhf.mods[i];
			p++;
		}

		public void addr(int c) {
			for (int i = 0; i < rhf.deg; i++)
				buf[i][p] = (buf[i][p - 1] + rhf.powers[i][p - 1] * c) % rhf.mods[i];
			p++;
		}

		public long queryTwin(int r) {
			return buf[0][r] << 32 | buf[1][r];
		}

		public long queryTwin(int l, int r) {
			assert l <= r;
			assert rhf.deg == 2;
			long h = 0;
			for (int i = 0; i < rhf.deg; i++) {
				long v = (buf[i][r] - buf[i][l] * rhf.powers[i][r - l]) % rhf.mods[i];
				if (v < 0)
					v += rhf.mods[i];
				h = h << 32 | v;
			}
			return h;
		}

		public long[] query(int l, int r) {
			assert l <= r;
			long[] h = new long[rhf.deg];
			for (int i = 0; i < rhf.deg; i++) {
				h[i] = (buf[i][r] - buf[i][l] * rhf.powers[i][r - l]) % rhf.mods[i];
				if (h[i] < 0)
					h[i] += rhf.mods[i];
			}
			return h;
		}

		public  long add(long a, long b, int w, RollingHashFactory rhf) {
			assert rhf.deg == 2;
			long high = ((a >>> 32) * rhf.powers[0][w] + (b >>> 32)) % rhf.mods[0];
			long low = ((long) (int) a * rhf.powers[1][w] + (int) b) % rhf.mods[1];
			return high << 32 | low;
		}
	}

	public  class RollingHashFactory {
		public int[] mods;
		public int[] muls;
		public long[][] powers;
		public int deg;

		public RollingHashFactory(int deg, int n, Random gen) {
			this.deg = deg;
			mods = new int[deg];
			muls = new int[deg];
			for (int i = 0; i < deg; i++) {
				mods[i] = BigInteger.probablePrime(30, gen).intValue();
				muls[i] = BigInteger.probablePrime(30, gen).intValue();
			}
			muls[0] = 100;
			mods[0] = 1000000007;
			powers = new long[deg][n + 1];
			for (int i = 0; i < deg; i++) {
				powers[i][0] = 1;
				for (int j = 1; j <= n; j++) {
					powers[i][j] = powers[i][j - 1] * muls[i] % mods[i];
				}
			}
		}
	}
}