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
