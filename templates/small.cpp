#include <bits/stdc++.h>
using namespace std;

/////////////////////////////////////////////////////////////////

const int MAXN = 1e6;
const int MOD = 1e9 + 7;
const int64_t INF64 = 1e18;
int n, a[MAXN + 1];

// parsing utility
vector<string> split_string(string s, string spl) { vector<string> ret; int pos = 0; while (pos < s.size()) { size_t nxt = s.find(spl, pos); if (nxt == string::npos) { ret.push_back(s.substr(pos)); break; } else if (nxt == pos) { ret.push_back(""); pos += spl.size(); continue; } else { ret.push_back(s.substr(pos, int(nxt) - pos)); pos = nxt + spl.size(); } } return ret; }

vector<int> parse_ints(string s) { vector<int> ret; string temp = ""; bool neg = false; s += "#"; for (char c : s) { if (c == '-') { neg = !neg; } else if ('0' <= c && c <= '9') { temp += c; } else { if (temp.size()) { int x = stoi(temp); if (neg) x *= -1; ret.push_back(x); } temp = ""; neg = false; } } return ret; }


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
