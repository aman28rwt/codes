import java.util.ArrayList;
import java.util.HashSet;

public class LCA_RMQ {
	
	public static class lcaHelp{
		int ind, val = Integer.MAX_VALUE;
	}
	
	public static ArrayList<Integer>[] tree;
	public static int[] level, euler, firstInd;
	public static lcaHelp[] segTree;
	public static int root, move;
	
	public static void intitLCA(int n) {
		level = new int[2 * n - 1];
		euler = new int[2 * n - 1];
		firstInd = new int[n + 1];
		segTree = new lcaHelp[8 * n + 1];
		move = 0;
		HashSet<Integer> vis = new HashSet<>();
		dfs(root, 0, vis);
		for(int i = 0; i < level.length; i++) {
			buildST(0, 0, level.length, i, level[i]);
		}
	}

	public static int LCA(int u, int v, int n) {
		if(firstInd[u] > firstInd[v]) {
			u = u ^ v;
			v = u ^ v;
			u = u ^ v;
		}
		
		lcaHelp rv = segQuery(0, 0, 2 * n - 1, firstInd[u], firstInd[v]);
		return euler[rv.ind];
	}
	
	public static int dfs(int node, int lvl, HashSet<Integer> vis) {
		euler[move] = node;
		level[move] = lvl;
		if(vis.contains(node)) {
			return -1;
		}
		vis.add(node);
		firstInd[node] = move;
		move++;
		
		for (Integer u : tree[node]) {
			int x = dfs(u, lvl + 1, vis);
			if(x == -1) {
				continue;
			}
			euler[move] = node;
			level[move] = lvl;
			move++;
		}
		
		return 0;
	}
	
	public static void buildST(int ind, int si, int li, int pos, int val) {
		if(segTree[ind] == null) {
			segTree[ind] = new lcaHelp();
		}
		
		if(si > li || si > pos || li < pos) {
			return;
		}
		
		if(si == li && si == pos) {
			segTree[ind].val = val;
			segTree[ind].ind = pos;
			return;
		}
		
		int mid = (si + li) >> 1;
		
		buildST(2 * ind + 1, si, mid, pos, val);
		buildST(2 * ind + 2, mid + 1, li, pos, val);
		
		if(segTree[2 * ind + 1].val < segTree[2 * ind + 2].val) {
			segTree[ind].val = segTree[2 * ind + 1].val;
			segTree[ind].ind = segTree[2 * ind + 1].ind;
		} else {
			segTree[ind].val = segTree[2 * ind + 2].val;
			segTree[ind].ind = segTree[2 * ind + 2].ind;
		}
	}

	public static lcaHelp segQuery(int ind, int si, int li, int l, int r) {
		if(l > li || r < si || si > li) {
			return new lcaHelp();
		}
		
		if(l <= si && r >= li) {
			return segTree[ind];
		}
		
		int mid = (si + li) >> 1;
		
		lcaHelp left = segQuery(2 * ind + 1, si, mid, l, r);
		lcaHelp right = segQuery(2 * ind + 2, mid + 1, li, l, r);
			
		lcaHelp rv = new lcaHelp();
		
		if(left.val < right.val) {
			rv.val = left.val;
			rv.ind = left.ind;
		} else {
			rv.val = right.val;
			rv.ind = right.ind;
		}
		
		return rv;
	}
}