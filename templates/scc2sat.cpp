struct SCC {
  int N, comps = 0;
  vector<int> scc, topo;
  vector<vector<int>> adj, radj;
  vector<bool> vis;
  SCC(int n) : N(n), scc(n), adj(n), radj(n), vis(n) {}
  void toposort(int cur) {
    if (vis[cur]) return;
    vis[cur] = true;
    for (int nbr : adj[cur]) toposort(nbr);
    topo.push_back(cur);
  }
  void rdfs(int cur) {
    if (vis[cur]) return;
    vis[cur] = true;
    scc[cur] = comps;
    for (int nbr : radj[cur]) rdfs(nbr);
  }
  // use these:
  void add_edge(int from, int to) {
    adj[from].push_back(to);
    radj[to].push_back(from);
  }
  void solve() {
    for (int i = 0; i < N; i++) {
      toposort(i);
    }
    vis.assign(N, 0);
    for (int i = int(topo.size()) - 1; i >= 0; i--) {
      if (!vis[topo[i]]) {
        rdfs(topo[i]);
        comps++;
      }
    }
  }
  int get(int node) {
    return scc[node];
  }
};

struct TWOSAT {
  SCC scc;
  int N;
  TWOSAT(int n) : N(2 * n), scc(2 * n) {}
  void addOR(int a, bool aa, int b, bool bb) {
    a *= 2, b *= 2;
    if (!aa) a ^= 1;
    if (!bb) b ^= 1;
    scc.add_edge(a ^ 1, b);
    scc.add_edge(b ^ 1, a);
  }
  bool solve() {
    scc.solve();
    for (int i = 0; i < N; i += 2) {
      if (scc.get(i) == scc.get(i + 1)) {
        return false;
      }
    }
    return true;
  }
  bool get(int node) {
    node *= 2;
    return scc.get(node) > scc.get(node + 1);
  }
};

int main() {
  cin.tie(0)->sync_with_stdio(0);
  int n, m;
  cin >> n >> m;
  TWOSAT sat(m);
  for (int i = 0; i < n; i++) {
    int a, b;
    bool aa = 1, bb = 1;
    char x, y;
    cin >> x >> a >> y >> b;
    a--, b--;
    if (x == '-') aa = 0;
    if (y == '-') bb = 0;
    sat.addOR(a, aa, b, bb);
  }
  if (!sat.solve()) {
    cout << "IMPOSSIBLE\n";
    return 0;
  }
  for (int i = 0; i < m; i++) {
    cout << (sat.get(i) ? '+' : '-') << ' ';
  }
}
