// seg2d st;
// st.build(vector<vector<int>>(5, vector<int>(5)));
// st.upd(0, 0, 23);
// st.query(0, 4, 2, 2); // row1, row2, col1, col2
struct seg2d {
  const static int MAXR = 1000, MAXC = 1000;
  using T = int;
  const T start = 0;
  T combine(T a, T b) {
    return a + b;
  }
  T t[2 * MAXR][2 * MAXC];
  int n = MAXR, m = MAXC;
  void build(vector<vector<T>>& a) {
    n = a.size();
    m = a[0].size();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < a[i].size(); j++) {
        t[i + n][j + m] = a[i][j];
      }
    }
    for (int i = n - 1; i >= 0; i--) {
      for (int j = m - 1; j >= 0; j--) {
        t[i + n][j] = combine(t[i + n][j << 1], t[i + n][j << 1 | 1]);
        t[i][j + m] = combine(t[i << 1][j + m], t[i << 1 | 1][j + m]);
      }
    }
    for (int i = n - 1; i > 0; i--) {
      for (int j = m - 1; j > 0; j--) {
        t[i][j] = combine(t[i << 1][j], t[i << 1 | 1][j]);
      }
    }
  }
  void upd(int r, int c, T v) {
    r += n, c += m;
    t[r][c] = v;
    for (int i = c; i >>= 1; ) {
      t[r][i] = combine(t[r][i << 1], t[r][i << 1 | 1]);
    }
    for (int j = r; j >>= 1; ) {
      t[j][c] = combine(t[j << 1][c], t[j << 1 | 1][c]);
      for (int i = c; i >>= 1; ) {
        t[j][i] = combine(t[j << 1][i], t[j << 1 | 1][i]);
      }
    }
  }
  T query1(int r, int a, int b) {
    T resl = start, resr = start;
    for (++b, a += m, b += m; a < b; a >>= 1, b >>= 1) {
      if (a & 1) resl = combine(resl, t[r][a++]);
      if (b & 1) resr = combine(t[r][--b], resr);
    }
    return combine(resl, resr);
  }
  T query(int r1, int r2, int a, int b) {
    T resl = start, resr = start;
    for (++r2, r1 += n, r2 += n; r1 < r2; r1 >>= 1, r2 >>= 1) {
      if (r1 & 1) resl = combine(resl, query1(r1++, a, b));
      if (r2 & 1) resr = combine(query1(--r2, a, b), resr);
    }
    return combine(resl, resr);
  }
};
