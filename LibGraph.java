import java.util.*;

/*
 * Warning: This code does not work yet and still needs to be debugged.
 * Hopefully I'll get around to finishing it tomorrow.
 */

class Tarjan {
	static ArrayList<Integer>[] adjLists;
	static ArrayList<Integer> S;
	static int indexCount;
	static int[] index, lowlink;
	static boolean[] onStack;
	public static ArrayList<ArrayList<Integer> > strongConnect(int v) {
		index[v] = indexCount;
		lowlink[v] = indexCount;
		++indexCount;
		S.add(v);
		onStack[v] = true;
		ArrayList<ArrayList<Integer> > components = new ArrayList<ArrayList<Integer> >();
		for (int w: adjLists[v]) {
			if (index[w] == -1) {
				strongConnect(w);
				lowlink[w] = Math.min(lowlink[v], lowlink[w]);
			} else if (onStack[w]) {
				lowlink[v] = Math.min(lowlink[v], index[w]);
			}
		}
		if (lowlink[v] == index[v]) {
			int w;
			ArrayList<Integer> comp = new ArrayList<Integer>();
			do {
				w = S.get(S.size()-1);
				S.remove(S.size()-1);
				onStack[w] = false;
				comp.add(w);
			} while (w != v);
			components.add(comp);
		}
		return components;
	}
	public static ArrayList<ArrayList<Integer> > getSCCs(ArrayList<Integer>[] G) {
		adjLists = G;
		indexCount = 0;
		S = new ArrayList<Integer>();
		index = new int[adjLists.length];
		lowlink = new int[adjLists.length];
		onStack = new boolean[adjLists.length];
		for (int i=0; i<adjLists.length; ++i) {
			index[i] = -1;
		}
		ArrayList<ArrayList<Integer> > sccList = new ArrayList<ArrayList<Integer> >();
		for (int v=0; v<adjLists.length; ++v) {
			if (index[v] == -1) {
				for (ArrayList<Integer> arr: strongConnect(v)) {
					sccList.add(arr);
				}
			}
		}
		return sccList;
	}
}

class Kosaraju {
	static boolean[] visited;
	static ArrayList<Integer> L;
	static ArrayList<Integer>[] adjLists, revAdjLists;
	static boolean[] assigned;
	public static ArrayList<ArrayList<Integer> > getSCCs(ArrayList<Integer>[] G) {
		adjLists = G;
		revAdjLists = new ArrayList[adjLists.length];
		for (int i=0; i<adjLists.length; ++i) {
			revAdjLists[i] = new ArrayList<Integer>();
		}
		for (int i=0; i<adjLists.length; ++i) {
			for (int j=0; j<adjLists[i].size(); ++j) {
				revAdjLists[adjLists[i].get(j)].add(i);
			}
		}
		visited = new boolean[adjLists.length];
		L = new ArrayList<Integer>();
		for (int u=0; u<adjLists.length; ++u) {
			Visit(u);
		}
		ArrayList<ArrayList<Integer> > comps = new ArrayList<ArrayList<Integer> >();
		Collections.reverse(L);
		assigned = new boolean[L.size()];
		for (int u=0; u<L.size(); ++u) {
			ArrayList<Integer> comp = new ArrayList<Integer>();
			Assign(u, comp);
			if (comp.size() > 0) {
				comps.add(comp);
			}
		}
		return comps;
	}
	public static void Visit(int u) {
		if (!visited[u]) {
			visited[u] = true;
			for (int v: adjLists[u]) {
				Visit(v);
				L.add(v);
			}
		}
	}
	public static void Assign(int u, ArrayList<Integer> comp) {
		if (!assigned[u]) {
			assigned[u] = true;
			comp.add(u);
			for (int v: revAdjLists[u]) {
				Assign(v, comp);
			}
		}
	}
}

class LibGraph {
	// given adjacency matrix of distances returns matrix of minimum distances
	public static int[][] FloydWarshall(int[][] distMat) {
		int n = distMat.length;
		int spaths[][] = new int[n][n];
		for (int i=0; i<n; ++i) {
			for (int j=0; j<n; ++j) {
				spaths[i][j] = distMat[i][j];
			}
		}
		for (int k=0; k<n; ++k) {
			for (int i=0; i<n; ++i) {
				for (int j=0; j<n; ++j) {
					if (spaths[i][k]+spaths[k][j] < spaths[i][j]) {
						spaths[i][j] = spaths[i][k]+spaths[k][j];
					}
				}
			}
		}
		return spaths;
	}
	// distances stored in distance, predecessors stored in pred after algorithm runs - these should be declared before calling
	public static void BellmanFord(ArrayList<Integer>[] adjLists, ArrayList<Integer>[] edgeWeights, int src, int[] distance, int[] pred) {
		for (int v=0; v<adjLists.length; ++v) {
			distance[v] = 1000000000;
			pred[v] = -1;
		}
		distance[src] = 0;
		
		for (int i=1; i<adjLists.length; ++i) {
			for (int u=0; u<adjLists.length; ++u) {
				for (int j=0; j<adjLists[u].size(); ++j) {
					int v = adjLists[u].get(j), w = edgeWeights[u].get(j);
					if (distance[u]+w < distance[v]) {
						distance[v] = distance[u]+w;
						pred[v] = u;
					}
				}
			}
		}
		
		for (int u=0; u<adjLists.length; ++u) {
			for (int j=0; j<adjLists[u].size(); ++j) {
				int v = adjLists[u].get(j), w = edgeWeights[u].get(j);
				if (distance[u]+w < distance[v]) {
					System.err.printf("negative cycle");
				}
			}
		}
	}
	
	public static void main(String[] args) {
		// testing code
		ArrayList<Integer>[] G = new ArrayList[6];
		for (int i=0; i<6; ++i) {
			G[i] = new ArrayList<Integer>();
		}
		G[0].add(1);
		G[1].add(2);
		G[2].add(0);
		G[3].add(4);
		G[4].add(3);
		G[1].add(3);
		G[5].add(3);
		System.out.println(Kosaraju.getSCCs(G));
	}
}
