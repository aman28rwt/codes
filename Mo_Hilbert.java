import java.util.Arrays;

public class Mo_Hilbert {
	void func() {
		int Q = 0; // inputs
		query[] qr = new query[Q]; // inputs // qr.l can be input or input - 1
		Arrays.sort(qr);
		int[] ans = new int[Q];
		int curL = 0, curR = -1; // indexes[0, N)
		for (query q : qr) {
			while (curR < q.r) {
				curR++;
				// WORK: add
			}
			while (curL > q.l) {
				curL--;
				// WORK: add
			}
			while (curR > q.r) {
				// WORK: remove
				curR--;
			}
			while (curL < q.l) {
				// WORK: remove
				curL++;
			}
			ans[q.ind] = 0; // req ans from work
		}
	}
	
	final int size = 21;
	public class query implements Comparable<query> {
		int l, r, ind;
		long hilbert;

		query(int x, int y, int i) {
			l = x;
			r = y;
			ind = i;
			hilbert = gilbertOrder(l, r, size, 0); // size = ceil(lg(n))
		}

		@Override
		public int compareTo(query o) {
			if (this.hilbert < o.hilbert) {
				return -1;
			} else if (this.hilbert > o.hilbert) {
				return 1;
			} else {
				return 0;
			}
		}

		long gilbertOrder(int x, int y, int pow, int rotate) {
			if (pow == 0) {
				return 0;
			}
			int hpow = 1 << (pow - 1);
			int seg = (x < hpow) ? ((y < hpow) ? 0 : 3) : ((y < hpow) ? 1 : 2);
			seg = (seg + rotate) & 3;
			int[] rotateDelta = { 3, 0, 0, 1 };
			int nx = x & (x ^ hpow), ny = y & (y ^ hpow);
			int nrot = (rotate + rotateDelta[seg]) & 3;
			long subSquareSize = 1L << (2 * pow - 2);
			long ans = seg * subSquareSize;
			long add = gilbertOrder(nx, ny, pow - 1, nrot);
			ans += (seg == 1 || seg == 2) ? add : (subSquareSize - add - 1);
			return ans;
		}
	}
}