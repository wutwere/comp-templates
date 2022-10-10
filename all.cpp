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


// segment tree: 1-indexed, builds from array a, upd to set
struct {
  int t[2 * MAXN] = {0};
  inline int combine(int l, int r) { return l + r; }
  void build() {
    for (int i = 0; i < n; i++)
      t[i + n] = a[i + 1];
    for (int i = n - 1; i > 0; --i)
      t[i] = combine(t[i << 1], t[i << 1 | 1]);
  }
  void upd(int p, int v) {
    for (--p, t[p += n] = v; p >>= 1; )
      t[p] = combine(t[p << 1], t[p << 1 | 1]);
  }
  int query(int l, int r) {
    int resl = 0, resr = 0;
    for (--l, l += n, r += n; l < r; l >>= 1, r >>= 1) {
      if (l & 1) resl = combine(resl, t[l++]);
      if (r & 1) resr = combine(t[--r], resr);
    }
    return combine(resl, resr);
  }
} seg_tree;


// get all shortest paths from a root in O(m log n)
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


/////////////////////////////////////////////////////////////////
int main() {
  ios::sync_with_stdio(0), cin.tie(0);
  cin >> n;
  vector<vector<array<int64_t, 2>>> edges;

  // floyd-warshall: get all shortest paths in O(n^3)
  vector dist(n + 1, vector<int64_t>(n + 1, INF64));
  for (int i = 1; i <= n; i++) {
    dist[i][i] = 0;
    for (auto [j, weight] : edges[i])
      dist[i][j] = weight;
  }
  for (int k = 1; k <= n; k++) {
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
      }
    }
  }

  ///////////////////////////////////////////////////////////////
  return 0;
}
