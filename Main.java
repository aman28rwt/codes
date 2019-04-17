import java.io.*;
import java.util.*;

public class Main implements Runnable {
	FastReader scn;
	PrintWriter out;
	String INPUT = "5\n" + 
			"2 4 8 3 6\n" + 
			"";

	void solve() {
		int n = scn.nextInt(), max = (int) 1e7 + 1;
		int[] arr = scn.nextIntArray(n);

		int[] lol = new int[max];
		for (int a : arr) {
			lol[a]++;
		}

		long lcm = Long.MAX_VALUE;
		int a = -1, b = -1;

		for (int i = 1; i < max; i++) {
			long z = Long.MAX_VALUE;
			int p = -1, q = -1;
			if (lol[i] == 1) {
				p = i;
			} else if (lol[i] > 1) {
				q = i;
				z = i;
				if (q != -1 && z < lcm) {
					lcm = z;
					a = p;
					b = q;
				}
				continue;
			}

			if (lol[i] != 0) {
				for (int j = 2 * i; j < max; j += i) {
					if (lol[j] > 0) {
						if (j < z) {
							z = j;
							q = j;
						}
						break;
					}
				}
			} else {
				p = i;
				while(p < max && lol[p] == 0) {
					p += i;
				}
				if(p > max) {
					continue;
				}
				for (int j = p + i; j < max; j += i) {
					if (lol[j] > 0) {
						long lc = (p * 1L * j) / gcd(p, j);
						if (lc < z) {
							z = lc;
							q = j;
						}
					}
				}
			}
			
			if (q != -1 && z < lcm) {
				lcm = z;
				a = p;
				b = q;
			}
		}

		if (a == -1) {
			int p = 1;
			while (lol[p] == 0) {
				p++;
			}
			for (int q = p + 1; q < max; q++) {
				if (lol[q] > 0) {
					long lc = (p * 1L * q) / gcd(p, q);
					if (lc < lcm) {
						lcm = lc;
						a = p;
						b = q;
					}
					p = q;
				}
			}
		}

		int i1 = -1, i2 = -1;
		for (int i = 0; i < n; i++) {
			if (arr[i] == a) {
				i1 = i + 1;
			} else if (arr[i] == b) {
				i2 = i + 1;
			}
		}

		out.println(i1 + " " + i2);
	}

	long gcd(long a, long b) {
		return b == 0 ? a : gcd(b, a % b);
	}

	public void run() {
		long time = System.currentTimeMillis();
		boolean oj = System.getProperty("ONLINE_JUDGE") != null;
		out = new PrintWriter(System.out);
		scn = new FastReader(oj);
		solve();
		out.flush();
		if (!oj) {
			System.out.println(Arrays.deepToString(new Object[] { System.currentTimeMillis() - time + " ms" }));
		}
	}

	public static void main(String[] args) {
		new Thread(null, new Main(), "Main", 1 << 26).start();
	}

	class FastReader {
		InputStream is;

		public FastReader(boolean onlineJudge) {
			is = onlineJudge ? System.in : new ByteArrayInputStream(INPUT.getBytes());
		}

		byte[] inbuf = new byte[1024];
		public int lenbuf = 0, ptrbuf = 0;

		int readByte() {
			if (lenbuf == -1)
				throw new InputMismatchException();
			if (ptrbuf >= lenbuf) {
				ptrbuf = 0;
				try {
					lenbuf = is.read(inbuf);
				} catch (IOException e) {
					throw new InputMismatchException();
				}
				if (lenbuf <= 0)
					return -1;
			}
			return inbuf[ptrbuf++];
		}

		boolean isSpaceChar(int c) {
			return !(c >= 33 && c <= 126);
		}

		int skip() {
			int b;
			while ((b = readByte()) != -1 && isSpaceChar(b))
				;
			return b;
		}

		double nextDouble() {
			return Double.parseDouble(next());
		}

		char nextChar() {
			return (char) skip();
		}

		String next() {
			int b = skip();
			StringBuilder sb = new StringBuilder();
			while (!(isSpaceChar(b))) { // when nextLine, (isSpaceChar(b) && b != ' ')
				sb.appendCodePoint(b);
				b = readByte();
			}
			return sb.toString();
		}

		String nextLine() {
			int b = skip();
			StringBuilder sb = new StringBuilder();
			while ((!isSpaceChar(b) || b == ' ')) { // when nextLine, (isSpaceChar(b) && b != ' ')
				sb.appendCodePoint(b);
				b = readByte();
			}
			return sb.toString();
		}

		char[] next(int n) {
			char[] buf = new char[n];
			int b = skip(), p = 0;
			while (p < n && !(isSpaceChar(b))) {
				buf[p++] = (char) b;
				b = readByte();
			}
			return n == p ? buf : Arrays.copyOf(buf, p);
		}

		int nextInt() {
			int num = 0, b;
			boolean minus = false;
			while ((b = readByte()) != -1 && !((b >= '0' && b <= '9') || b == '-'))
				;
			if (b == '-') {
				minus = true;
				b = readByte();
			}

			while (true) {
				if (b >= '0' && b <= '9') {
					num = num * 10 + (b - '0');
				} else {
					return minus ? -num : num;
				}
				b = readByte();
			}
		}

		long nextLong() {
			long num = 0;
			int b;
			boolean minus = false;
			while ((b = readByte()) != -1 && !((b >= '0' && b <= '9') || b == '-'))
				;
			if (b == '-') {
				minus = true;
				b = readByte();
			}

			while (true) {
				if (b >= '0' && b <= '9') {
					num = num * 10 + (b - '0');
				} else {
					return minus ? -num : num;
				}
				b = readByte();
			}
		}

		char[][] nextMatrix(int n, int m) {
			char[][] map = new char[n][];
			for (int i = 0; i < n; i++)
				map[i] = next(m);
			return map;
		}

		int[] nextIntArray(int n) {
			int[] a = new int[n];
			for (int i = 0; i < n; i++)
				a[i] = nextInt();
			return a;
		}

		long[] nextLongArray(int n) {
			long[] a = new long[n];
			for (int i = 0; i < n; i++)
				a[i] = nextLong();
			return a;
		}

		int[][] next2DInt(int n, int m) {
			int[][] arr = new int[n][];
			for (int i = 0; i < n; i++) {
				arr[i] = nextIntArray(m);
			}
			return arr;
		}

		long[][] next2DLong(int n, int m) {
			long[][] arr = new long[n][];
			for (int i = 0; i < n; i++) {
				arr[i] = nextLongArray(m);
			}
			return arr;
		}

		int[] shuffle(int[] arr) {
			Random r = new Random();
			for (int i = 1, j; i < arr.length; i++) {
				j = r.nextInt(i);
				arr[i] = arr[i] ^ arr[j];
				arr[j] = arr[i] ^ arr[j];
				arr[i] = arr[i] ^ arr[j];
			}
			return arr;
		}

		int[] uniq(int[] arr) {
			Arrays.sort(arr);
			int[] rv = new int[arr.length];
			int pos = 0;
			rv[pos++] = arr[0];
			for (int i = 1; i < arr.length; i++) {
				if (arr[i] != arr[i - 1]) {
					rv[pos++] = arr[i];
				}
			}
			return Arrays.copyOf(rv, pos);
		}

		int[] reverse(int[] arr) {
			int l = 0, r = arr.length - 1;
			while (l < r) {
				arr[l] = arr[l] ^ arr[r];
				arr[r] = arr[l] ^ arr[r];
				arr[l] = arr[l] ^ arr[r];
				l++;
				r--;
			}
			return arr;
		}
	}
}