#include <bits/stdc++.h>
using namespace std;

/////////////////////////////////////////////////////////////////

const int MAXN = 1e6;
const int MOD = 1e9 + 7;
const int64_t INF64 = 1e18;
int n, a[MAXN + 1];


// x to the y-th power mod MOD in O(log y)
int64_t powmod(int64_t x, int64_t y) {
  if (y == 0) return 1;
  int64_t res = powmod(x * x % MOD, y / 2);
  if (y & 1) res *= x;
  return res % MOD;
}


// n choose r mod MOD in O(1) after O(n) precomputation
int64_t fac[MAXN + 1], finv[MAXN + 1], inv[MAXN + 1];
int64_t ncr(int n, int r) {
  if (fac[0] == 0) {
    fac[0] = finv[0] = inv[1] = 1;
    for (int i = 2; i <= MAXN; i++) {
      inv[i] = inv[MOD % i] * (-MOD / i) % MOD;
      if (inv[i] < 0) inv[i] += MOD;
    }
    for (int i = 1; i <= MAXN; i++) {
      fac[i] = fac[i - 1] * i % MOD;
      finv[i] = finv[i - 1] * inv[i] % MOD;
    }
  }
  return (n < 0 || r < 0 || r > n ? 0 : fac[n] * finv[r] % MOD * finv[n - r] % MOD);
}


// seg_tree<int> st([](int a, int b) { return a + b; });
// vector<int> a(100, 5);
// st.build(a);
// st.upd(0, 10);
// st.query(0, 99);
template <class T>
struct seg_tree {
  const static int MAXN = 100000;
  function<T(T, T)> combine;
  T t[2 * MAXN];
  int n = MAXN;
  seg_tree(function<T(T, T)> f) { combine = f; }
  void build(vector<T>& a) {
    n = a.size();
    for (int i = 0; i < n; i++)
      t[i + n] = a[i];
    for (int i = n - 1; i > 0; --i)
      t[i] = combine(t[i << 1], t[i << 1 | 1]);
  }
  void upd(int p, T v) {
    for (t[p += n] = v; p >>= 1; )
      t[p] = combine(t[p << 1], t[p << 1 | 1]);
  }
  T query(int l, int r) {
    T resl = 0, resr = 0;
    for (++r, l += n, r += n; l < r; l >>= 1, r >>= 1) {
      if (l & 1) resl = combine(resl, t[l++]);
      if (r & 1) resr = combine(t[--r], resr);
    }
    return combine(resl, resr);
  }
};


// vector<int64_t> min_dists = dijkstra(1, edges);
vector<int64_t> dijkstra(int root, vector<vector<array<int64_t, 2>>> &edges) {
  vector<int64_t> minDist(n + 1, INF64);
  priority_queue<array<int64_t, 2>> pq;
  minDist[root] = 0;
  pq.push({0, root});
  while (!pq.empty()) {
    auto [dist, node] = pq.top(); dist *= -1; pq.pop();
    if (minDist[node] != dist) continue;
    for (auto [nbr, weight] : edges[node]) {
      if (dist + weight < minDist[nbr]) {
        minDist[nbr] = dist + weight;
        pq.push({-minDist[nbr], nbr});
      }
    }
  }
  return minDist;
}

// disjoint_set<100000> ds;
// ds.join(1, 2);
// -ds.par[ds.find(1)];
template <int N>
struct disjoint_set {
  int par[N + 1];
  disjoint_set() { memset(par, -1, sizeof(par)); }
  int find(int v) { return par[v] < 0 ? v : par[v] = find(par[v]); }
  void join(int u, int v) {
    if ((u = find(u)) != (v = find(v))) {
      if (par[u] > par[v]) swap(u, v);
      par[u] += par[v], par[v] = u;
    }
  }
};

/////////////////////////////////////////////////////////////////
int main() {
  ios::sync_with_stdio(0), cin.tie(0);
  cin >> n;


  // sieve[i] = true if i is prime
  vector<bool> sieve(MAXN + 1, true);
  sieve[0] = sieve[1] = false;
  for (int i = 2; i * i <= MAXN; i++)
    for (int j = i * i; sieve[i] && j <= MAXN; j += i)
      sieve[j] = false;
  vector<int> primes = {2};
  for (int i = 3; i <= MAXN; i += 3)
    if (sieve[i])
      primes.push_back(i);


  // floyd-warshall: get all shortest paths in O(n^3)
  vector<vector<array<int64_t, 2>>> edges;
  vector dist(n + 1, vector<int64_t>(n + 1, INF64));
  for (int i = 1; i <= n; i++) {
    dist[i][i] = 0;
    for (auto [j, weight] : edges[i])
      dist[i][j] = weight;
  }
  for (int k = 1; k <= n; k++)
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= n; j++)
        dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);


  // topological sort with cycle detection in O(n)
  vector<vector<int>> out_edges(n + 1);
  vector<bool> visited(n + 1);
  vector<int> topo, topoIdx(n + 1);
  function<void(int)> dfs = [&](int cur) {
    if (visited[cur]) return;
    visited[cur] = true;
    for (int nbr : out_edges[cur]) dfs(nbr);
    topo.push_back(cur);
  };
  for (int i = 1; i <= n; i++) if (!visited[i]) dfs(i);
  reverse(topo.begin(), topo.end());
  for (int i = 0; i < n; i++) topoIdx[topo[i]] = i;
  bool cycle = false;
  for (int i = 1; !cycle && i <= n; i++)
    for (int nbr : out_edges[i])
      if (topoIdx[i] > topoIdx[nbr]) {
        cycle = true;
        break;
      }

  ///////////////////////////////////////////////////////////////
  return 0;
}
