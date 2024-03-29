template <class T = int>
struct RMQ {
  vector<vector<T>> jmp;  // 0-idx
  RMQ(const vector<T>& V) : jmp(1, V) {
    for (int pw = 1, k = 1; pw * 2 <= V.size(); pw *= 2, ++k) {
      jmp.emplace_back(int(V.size()) - pw * 2 + 1);
      for (int j = 0; j < jmp[k].size(); j++)
        jmp[k][j] = min(jmp[k - 1][j], jmp[k - 1][j + pw]);
    }
  }
  T query(int a, int b) { // [a, b)
    assert(a < b);  // or return inf if a == b
    int dep = 31 - __builtin_clz(b - a);
    return min(jmp[dep][a], jmp[dep][b - (1 << dep)]);
  }
};
// 0-idx
// make vector<vector<int>> adj
// then build LCA lca(adj, root)
// then l = lca.lca(u, v)
struct LCA {
  int tim = 0;
  vector<int> time, path, ret;
  RMQ<int> rmq;
  LCA(vector<vector<int>>& g, int rt) : time(g.size()), rmq((dfs(g, rt, -1), ret)) {}
  void dfs(vector<vector<int>>& g, int u, int fa) {
    time[u] = tim++;
    for (int v : g[u]) {
      if (v == fa) continue;
      path.push_back(u), ret.push_back(time[u]);
      dfs(g, v, u);
    }
  }
  int lca(int a, int b) {
    if (a == b) return a;
    tie(a, b) = minmax(time[a], time[b]);
    return path[rmq.query(a, b)];
  }
  // dist(a, b) { return depth[a] + depth[b] - 2 * depth[lca(a, b)]; }
};
