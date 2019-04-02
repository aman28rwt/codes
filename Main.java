import java.io.*;
import java.util.*;

public class Main implements Runnable {
	FastReader scn;
	StringBuilder out;
	String INPUT = "";

	void solve() {
		int n = scn.nextInt(), m = scn.nextInt();
		int[] from = new int[m], to = new int[m], w = new int[m];
		for (int i = 0; i < m; i++) {
			from[i] = scn.nextInt() - 1;
			to[i] = scn.nextInt() - 1;
			w[i] = scn.nextInt();
		}
		int[][][] g = packWU(n, from, to, w);
		if (!dfs(g, 0, n - 1, new boolean[n])) {
			out.append(-1 + "\n");
			return;
		}
		
		int[] par = dijkstra(g, 0);
		int[][] gg = parentToG(par);
		dfs(gg, 0, -1, n - 1);
		for (int i = ans.size() - 1; i >= 0; i--) {
			out.append(ans.get(i) + 1 + " ");
		}
		out.append("\n");
	}
	
	int[] dijkstra(int[][][] g, int from) {
		int n = g.length;
		
		long[] dist = new long[n]; // Contains the distance
		int[] par = new int[n]; // Contains the tree formed (not necessary MST). It contains parent of node I

		Arrays.fill(dist, Long.MAX_VALUE / 10);
		Arrays.fill(par, -1);
		dist[from] = 0;
		
		PriorityQueue<long[]> pq = new PriorityQueue<>((o1, o2) -> Long.compare(o1[1], o2[1]));
		pq.add(new long[] {from, 0});
		
		while(!pq.isEmpty()) {
			int cur = (int)pq.poll()[0];
			
			for (int[] e : g[cur]) {
				int next = e[0];
				long newDist = dist[cur] + e[1];
				if (newDist < dist[next]) {
					dist[next] = newDist;
					par[next] = cur;
					pq.add(new long[] {next, newDist});
				} else if (newDist == dist[next]) {
					if (dist[cur] > dist[par[next]]) {
						par[next] = cur;
					}
				}
			}
		}
		
		return par;
	}
	
	ArrayList<Integer> ans = new ArrayList<>();

	boolean dfs(int[][] g, int u, int p, int des) {
		if (u == des) {
			ans.add(des);
			return true;
		}

		for (int v : g[u]) {
			if (v != p) {
				if (dfs(g, v, u, des)) {
					ans.add(u);
					return true;
				}
			}
		}

		return false;
	}
	
	int[][] parentToG(int[] par) {
		int n = par.length;
		int[] ct = new int[n];
		for (int i = 0; i < n; i++) {
			if (par[i] >= 0) {
				ct[i]++;
				ct[par[i]]++;
			}
		}
		int[][] g = new int[n][];
		for (int i = 0; i < n; i++) {
			g[i] = new int[ct[i]];
		}
		for (int i = 0; i < n; i++) {
			if (par[i] >= 0) {
				g[par[i]][--ct[par[i]]] = i;
				g[i][--ct[i]] = par[i];
			}
		}
		return g;
	}
	
	boolean dfs(int[][][] g, int u, int des, boolean[] vis) {
		if (vis[u]) {
			return false;
		}
		vis[u] = true;
		if (u == des) {
			return true;
		}

		for (int[] v : g[u]) {
			if (dfs(g, v[0], des, vis)) {
				return true;
			}
		}

		return false;
	}
	
	int[][][] packWU(int n, int[] from, int[] to, int[] w) {
		int[][][] g = new int[n][][];
		int[] p = new int[n];
		for (int f : from)
			p[f]++;
		for (int t : to)
			p[t]++;
		for (int i = 0; i < n; i++)
			g[i] = new int[p[i]][2];
		for (int i = 0; i < from.length; i++) {
			--p[from[i]];
			g[from[i]][p[from[i]]][0] = to[i];
			g[from[i]][p[from[i]]][1] = w[i];
			--p[to[i]];
			g[to[i]][p[to[i]]][0] = from[i];
			g[to[i]][p[to[i]]][1] = w[i];
		}
		return g;
	}

	public void run() {
		long time = System.currentTimeMillis();
		boolean oj = System.getProperty("ONLINE_JUDGE") != null;
		out = new StringBuilder();
		scn = new FastReader(oj);
		solve();
		System.out.print(out);
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