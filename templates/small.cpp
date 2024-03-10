#include <bits/stdc++.h>
using namespace std;

// indexed_set v;
// v.insert(5);
// *v.find_by_order(0);
// v.order_of_key(7);
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> indexed_set;

/////////////////////////////////////////////////////////////////

const int MAXN = 1e6;
const int MOD = 1e9 + 7;
const int64_t INF64 = 1e18;
int n, a[MAXN + 1];

// parsing utility
vector<string> split_string(string s, string spl) { vector<string> ret; int pos = 0; while (pos < s.size()) { size_t nxt = s.find(spl, pos); if (nxt == string::npos) { ret.push_back(s.substr(pos)); break; } else if (nxt == pos) { ret.push_back(""); pos += spl.size(); continue; } else { ret.push_back(s.substr(pos, int(nxt) - pos)); pos = nxt + spl.size(); } } return ret; }

vector<int> parse_ints(string s) { vector<int> ret; string temp = ""; bool neg = false; s += "#"; for (char c : s) { if (c == '-') { neg = !neg; } else if ('0' <= c && c <= '9') { temp += c; } else { if (temp.size()) { int x = stoi(temp); if (neg) x *= -1; ret.push_back(x); } temp = ""; neg = false; } } return ret; }


// x to the y-th power mod MOD in O(log y)
int powmod(int x, int y) {
  int res = 1;
  for (; y; y >>= 1, x = 1ll * x * x % MOD)
    if (y & 1) res = 1ll * res * x % MOD;
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


// RMQ + LCA
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


// vector<int64_t> min_dists = dijkstra(1, edges);
vector<int64_t> dijkstra(int root, vector<vector<array<int64_t, 2>>> &edges) {
  vector<int64_t> minDist(n + 1, INF64);
  priority_queue<array<int64_t, 2>, vector<array<int64_t, 2>>, greater<array<int64_t, 2>>> pq;
  minDist[root] = 0;
  pq.push({0, root});
  while (!pq.empty()) {
    auto [dist, node] = pq.top(); pq.pop();
    if (minDist[node] != dist) continue;
    for (auto [nbr, weight] : edges[node]) {
      if (dist + weight < minDist[nbr]) {
        minDist[nbr] = dist + weight;
        pq.push({minDist[nbr], nbr});
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


// vector<int> zarr = z(b + "#" + a);
// everywhere where zarr[i] == b.size() means b appears
vector<int> z(string s) {
  int n = s.size();
  vector<int> z(n);
  int x = 0, y = 0;
  for (int i = 1; i < n; i++) {
    z[i] = max(0, min(z[i - x], y - i + 1));
    while (i + z[i] < n && s[z[i]] == s[i + z[i]]) {
      x = i; y = i + z[i]++;
    }
  }
  return z;
}


// vector<int> m = manacher(s);
// m[2 * i] for even length palindrome centered between i - 1 and i
// m[2 * i + 1] for odd length palindrome centered at i
vector<int> manacher(string &s) {
  string odd = "#";
  for (char c : s) odd += c, odd += '#';
  odd = "$" + odd + "^";
  vector<int> m(odd.size());
  for (int i = 1, l = 1, r = 1; i < odd.size(); i++) {
    if (i < r) m[i] = max(0, min(r - i, m[l + r - i]));
    if (i + m[i] >= r) while (odd[i - m[i] - 1] == odd[i + m[i] + 1])
      m[i]++, l = i - m[i], r = i + m[i];
  }
  return vector<int>(m.begin() + 1, m.end() - 2);
}


struct custom_hash {
  static uint64_t splitmix64(uint64_t x) {
    // http://xorshift.di.unimi.it/splitmix64.c
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }

  size_t operator()(uint64_t x) const {
    static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
    return splitmix64(x + FIXED_RANDOM);
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
  for (int i = 3; i <= MAXN; i += 2)
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
}
