#include <bits/stdc++.h>
using namespace std;

// BEGIN


struct MCMF {
  static const int MAXN = 1010, MAXM = 2010, INF = 1e9;
  struct edge {
    int node, to, f, c;
  };
  edge e[MAXM * 2];
  int S, T, m = 1, hd[MAXN], d[MAXN], pre[MAXN], pid[MAXN];
  vector<int> cur_path;
  bool vis[MAXN], in[MAXN];
  void add_edge(int u, int v, int f, int c) {
    e[++m] = {v, hd[u], f, c};
    e[++m] = {u, hd[v], 0, -c};
    hd[u] = m - 1;
    hd[v] = m;
  }
  void print_path(int cur, int goal) {
    cur_path.push_back(cur);
    if (cur == goal) {
      cout << cur_path.size() << '\n';
      for (int i = 0; i < cur_path.size(); i++) {
        cout << cur_path[i] << " \n"[i == cur_path.size() - 1];
      }
    } else {
      for (int i = 2; i <= m; i++) {
        if (e[i].node == cur && i % 2 && e[i ^ 1].f == 0) {
          e[i ^ 1].f = 1;
          print_path(e[i ^ 1].node, goal);
          return;
        }
      }
    }
  }
  bool SPFA() {
    memset(d, 0x3f, sizeof(d));
    queue<int> q;
    in[S] = 1;
    d[S] = 0;
    q.push(S);
    while (!q.empty()) {
      int u = q.front();
      q.pop();
      in[u] = 0;
      for (int i = hd[u]; i; i = e[i].to) {
        if (!e[i].f || d[e[i].node] <= d[u] + e[i].c) continue;
        d[e[i].node] = d[u] + e[i].c;
        pre[e[i].node] = u;
        pid[e[i].node] = i;
        if (!in[e[i].node]) {
          in[e[i].node] = 1;
          q.push(e[i].node);
        }
      }
    }
    return d[T] < INF;
  }
  array<int, 2> expath() {
    int cur = T, flow = INF, cost = 0;
    while (cur != S) {
      flow = min(flow, e[pid[cur]].f);
      cur = pre[cur];
    }
    cur = T;
    while (cur != S) {
      int nxt = pid[cur];
      cur = pre[cur];
      e[nxt].f -= flow;
      e[nxt ^ 1].f += flow;
      cost += flow * e[nxt].c;
    }
    return {flow, cost};
  }
  array<int, 2> mcmf() { // return {flow, cost}
    array<int, 2> ans;
    while (SPFA()) {
      auto [f, c] = expath();
      ans[0] += f;
      ans[1] += c;
    }
    return ans;
  }
} G;

// END

int main() {
	cin.tie(0)->sync_with_stdio(0);
	for (int __, _ = (cin >> __, 0); ++_ <= __;) {
		cout << "Case #" << _ << ": ";
	}
}
