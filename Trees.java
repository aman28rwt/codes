import java.util.*;

public class Trees {
	int[][] packU(int n, int[] from, int[] to) {
		int[][] g = new int[n][];
		int[] p = new int[n];
		for (int f : from)
			p[f]++;
		for (int t : to)
			p[t]++;
		for (int i = 0; i < n; i++)
			g[i] = new int[p[i]];
		for (int i = 0; i < from.length; i++) {
			g[from[i]][--p[from[i]]] = to[i];
			g[to[i]][--p[to[i]]] = from[i];
		}
		return g;
	}

	int[][] packD(int n, int[] from, int[] to) {
		int[][] g = new int[n][];
		int[] p = new int[n];
		for (int f : from)
			p[f]++;
		for (int i = 0; i < n; i++)
			g[i] = new int[p[i]];
		for (int i = 0; i < from.length; i++) {
			g[from[i]][--p[from[i]]] = to[i];
		}
		return g;
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

	int[][][] packWD(int n, int[] from, int[] to, int[] w) {
		int[][][] g = new int[n][][];
		int[] p = new int[n];
		for (int f : from)
			p[f]++;
		for (int i = 0; i < n; i++)
			g[i] = new int[p[i]][2];
		for (int i = 0; i < from.length; i++) {
			--p[from[i]];
			g[from[i]][p[from[i]]][0] = to[i];
			g[from[i]][p[from[i]]][1] = w[i];
		}
		return g;
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
	
	int[] sortTopologically(int[][] g) {
		int n = g.length;
		int[] ec = new int[n];
		for (int i = 0; i < n; i++) {
			for (int to : g[i])
				ec[to]++;
		}
		int[] ret = new int[n];
		int q = 0;

		// sources
		for (int i = 0; i < n; i++) {
			if (ec[i] == 0)
				ret[q++] = i;
		}

		for (int p = 0; p < q; p++) {
			for (int to : g[ret[p]]) {
				if (--ec[to] == 0)
					ret[q++] = to;
			}
		}
		// loop
		for (int i = 0; i < n; i++) {
			if (ec[i] > 0)
				return null;
		}
		return ret;
	}

	// KOSARAJU
	int[] decomposeToSCC(int[][] g, int[][] ig) {
		int n = g.length;
		boolean[] visited = new boolean[n];
		int[] po = new int[n];
		int pop = 0;
		int[] stack = new int[n];
		int[] sinds = new int[n];
		int sp = 0;
		for (int i = 0; i < n; i++) {
			if (!visited[i]) {
				sinds[sp] = 0;
				stack[sp++] = i;
				while (sp > 0) {
					int cur = stack[sp - 1];
					visited[cur] = true;
					while (sinds[sp - 1] < g[cur].length && visited[g[cur][sinds[sp - 1]]])
						sinds[sp - 1]++;
					if (sinds[sp - 1] == g[cur].length) {
						po[pop++] = cur;
						sp--;
					} else {
						stack[sp] = g[cur][sinds[sp - 1]];
						sinds[sp] = 0;
						sp++;
					}
				}
			}
		}
		int[] ret = new int[n];
		Arrays.fill(visited, false);
		int clus = 0;
		Queue<Integer> q = new ArrayDeque<Integer>();
		for (int i = n - 1; i >= 0; i--) {
			if (!visited[po[i]]) {
				q.add(po[i]);
				visited[po[i]] = true;
				while (!q.isEmpty()) {
					int cur = q.poll();
					ret[cur] = clus;
					for (int k : ig[cur]) {
						if (!visited[k]) {
							q.add(k);
							visited[k] = true;
						}
					}
				}
				clus++;
			}
		}
		return ret;
	}

	// CENTROID
	void centroidT(int[][] g, int[] cnt, int u, int p, int[] par, boolean[] vis) {
		int max = 0, v = -1;
		for (int x : g[u]) {
			if (!vis[x] && cnt[x] > max) {
				max = cnt[x];
				v = x;
			}
		}
		while (max > cnt[u] / 2) {
			int t = cnt[v];
			cnt[v] = cnt[u];
			cnt[u] = cnt[v] - t;
			u = v;
			max = 0;
			for (int x : g[u]) {
				if (!vis[x] && cnt[x] > max) {
					max = cnt[x];
					v = x;
				}
			}
		}
		par[u] = p;
		vis[u] = true;
		for (int x : g[u]) {
			if (!vis[x]) {
				centroidT(g, cnt, x, u, par, vis);
			}
		}
	}

	int dfs_count(int[][] g, int u, int p, int[] cnt) {
		cnt[u] = 1; // can be weight of edge
		for (int v : g[u]) {
			if (v != p) {
				cnt[u] += dfs_count(g, v, u, cnt);
			}
		}
		return cnt[u];
	}
	// CENTROID

	// [0] gives parents, [1] gives bfs order, [2] gives depth
	int[][] parents3(int[][] g, int root) {
		int n = g.length;
		int[] par = new int[n];
		Arrays.fill(par, -1);

		int[] depth = new int[n];
		depth[0] = 0;

		int[] q = new int[n];
		q[0] = root;
		for (int p = 0, r = 1; p < r; p++) {
			int cur = q[p];
			for (int nex : g[cur]) {
				if (par[cur] != nex) {
					q[r++] = nex;
					par[nex] = cur;
					depth[nex] = depth[cur] + 1;
				}
			}
		}
		return new int[][] { par, q, depth };
	}

	// Used for Weighted graph
	int[][] parents(int[][][] g, int root) {
		int n = g.length;
		int[] par = new int[n];
		Arrays.fill(par, -1);
		int[] dw = new int[n];
		int[] pw = new int[n];
		int[] dep = new int[n];

		int[] q = new int[n];
		q[0] = root;
		for (int p = 0, r = 1; p < r; p++) {
			int cur = q[p];
			for (int[] nex : g[cur]) {
				if (par[cur] != nex[0]) {
					q[r++] = nex[0];
					par[nex[0]] = cur;
					dep[nex[0]] = dep[cur] + 1;
					dw[nex[0]] = Math.max(dw[cur], nex[1]);
					pw[nex[0]] = nex[1];
				}
			}
		}
		return new int[][] { par, q, dep, dw, pw };
	}

	// Sparse Table
	int[][] logstepParents(int[] parent) {
		int n = parent.length;
		int m = Integer.numberOfTrailingZeros(Integer.highestOneBit(n - 1)) + 1;
		int[][] sparse = new int[m][n];
		sparse[0] = parent;
		for (int j = 1; j < m; j++) {
			for (int i = 0; i < n; i++) {
				sparse[j][i] = sparse[j - 1][i] == -1 ? -1 : sparse[j - 1][sparse[j - 1][i]];
			}
		}
		return sparse;
	}
	// Sparse Table

	// Weighted Tree Sparse : pw can be obtained from above functions
	long[][] logstepSums(int[][] sp, int[] pw) {
		int m = sp.length;
		int n = sp[0].length;
		long[][] sums = new long[m][n];
		for (int i = 0; i < n; i++) {
			sums[0][i] = pw[i];
		}
		for (int j = 1; j < m; j++) {
			for (int i = 0; i < n; i++) {
				if (sp[j - 1][i] != -1) {
					sums[j][i] = sums[j - 1][i] + sums[j - 1][sp[j - 1][i]];
				} else {
					sums[j][i] = 0;
				}
			}
		}
		return sums;
	}
	// Weighted Tree Sparse

	// LCA
	int lca2(int a, int b, int[][] sparse, int[] depth) {
		if (depth[a] < depth[b]) {
			b = ancestor(b, depth[b] - depth[a], sparse);
		} else if (depth[a] > depth[b]) {
			a = ancestor(a, depth[a] - depth[b], sparse);
		}
		if (a == b) {
			return a;
		}
		int sa = a, sb = b;
		for (int low = 0, high = depth[a], t = Integer.highestOneBit(high), k = Integer
				.numberOfTrailingZeros(t); t > 0; t >>>= 1, k--) {
			if ((low ^ high) >= t) {
				if (sparse[k][sa] != sparse[k][sb]) {
					low |= t;
					sa = sparse[k][sa];
					sb = sparse[k][sb];
				} else {
					high = low | t - 1;
				}
			}
		}
		return sparse[0][sa];
	}

	int ancestor(int u, int dep, int[][] sparse) {
		for (int i = 0; dep > 0 && u != -1; dep >>>= 1, i++) {
			if ((dep & 1) == 1)
				u = sparse[i][u];
		}
		return u;
	}
	// LCA

	// DJIKSTRA
	// For weighted graph
	// If graph is not weighted assume weights to be 1
	long[] dijkstra(int[][][] g, int from) {
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
		
		return dist;
	}

	// Heavy Light Decomposition : HLD

	// Few tips and tricks ;)

	// in[i] : stores the time at which node i comes
	// size[i] : stores the size of subtree rooted at i
	// order : stores DFS order

	// <<<<<<<<< IMPORTANT <<<<<<<<<< IMPORTANT >>>>>>>>>> IMPORTANT >>>>>>>>>

	// next[i] : stores the node at which chain from i to root breaks and remaining
	// chain starts from parent[next[i]] (i.e. i = par[next[i]]). This goes on until
	// parent[i] is not equal -1

	// Subtree of node u is from [in[u], in[u] + size[u]) in array : order

	// Heavy path from node u to root is :
	// while u != -1
	// [in[next[u], in[u]]
	// u = par[next[u]]

	// Path from u to v is path from u to LCA(u, v) and LCA(u, v) to v

	// Path from u to v (where v is in same branch which goes to root
	// and depth[u] >= depth[v])
	// while depth[u] >= depth[v]
	// [max(in[next[u], in[v]), in[u]]
	// u = par[in[next[u]]

	int time = 0;

	void dfs_hld(int[][] g, int[] in, int[] order, int[] next, int u, int p) {
		in[u] = time++;
		order[in[u]] = u;
		for (int v : g[u]) {
			if (v != p) {
				next[v] = (v == g[u][0] ? next[u] : v);
				dfs_hld(g, in, order, next, v, u);
			}
		}
	}

	void dfs_size(int[][] g, int[] size, int[] parent, int[] depth, int u, int p, int d) {
		size[u] = 1;
		parent[u] = p;
		depth[u] = d;

		if (g[u].length > 1 && g[u][0] == p) {
			int temp = g[u][0];
			g[u][0] = g[u][1];
			g[u][1] = temp;
		}

		for (int i = 0; i < g[u].length; i++) {
			if (g[u][i] != p) {
				dfs_size(g, size, parent, depth, g[u][i], u, d + 1);
				size[u] += size[g[u][i]];
				if (size[g[u][i]] > size[g[u][0]]) {
					int temp = g[u][0];
					g[u][0] = g[u][i];
					g[u][i] = temp;
				}
			}
		}
	}
	// Heavy Light Decomposition

	// Minimum Spanning Tree
	// ----------------------------------------------------------------------
	// 1. To make MST from weighted graph use DSU
	// 2. Sort edges based on their W
	// 3. For each edge e, check their vertices are connected or not in DSU
	// 3.1 if connected, continue
	// 3.2 else, connect them, add that edge to required graph(to, from, w)
	// 4. DONE!
	// ----------------------------------------------------------------------
	// Minimum Spanning Tree
}