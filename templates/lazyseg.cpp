// seg_tree<int,int> sumupd_sumqry(
//  [](int val1, int val2) { return val1 + val2; }, 
//  [](int qry1, int qry2) { return qry1 + qry2; },
//  [](int qry, int val, int width) { return val + qry * width; });
// vector<int> a(100, 5);
// st.build(a, 0); // Second arg is default lazy node
// st.upd(0, 10, 3);
// st.query(0, 99);
template <class T, class Q>
struct seg_tree {
  const static int MAXN = 2e5;
  function<T(T,T)> combT;
  function<Q(Q,Q)> combQ;
  function<T(Q,T,int)> apply;
  T t[2 * MAXN];
  Q q[MAXN];
  int w[MAXN];
  int n, h;
  seg_tree(
    function<T(T,T)> combT,
    function<Q(Q,Q)> combQ,
    function<T(Q,T,int)> apply
  ) : combT(combT), combQ(combQ), apply(apply) {
  }
  inline void build(vector<T>& a, Q def) {
    n = a.size();
    h = (n>1)?32-__builtin_clz(n-1):1;
    for (int i = 0; i < n; i++) {
      t[i + n] = a[i], q[i] = def;
      int o = h+__builtin_clz(i)-32;
      w[i] = (1<<o)+max(0, min(1<<o,n-(i<<o)));
    }
    for (int i = n - 1; i > 0; --i)
      t[i] = combT(t[i << 1], t[i << 1 | 1]);
  }
  inline int cw(int i) {
    return (i<n)?w[i]:1;
  }
  inline void prpd(int i) {
    int j = i<<1;
    if (j<n)q[j]=combQ(q[i],q[j]),q[j|1]=combQ(q[i],q[j|1]); 
    t[j]=apply(q[i],t[j],cw(j)),t[j|1]=apply(q[i],t[j|1],cw(j|1));
    q[i]=q[0];
  }
  inline void prpu(int i) {
    if(i>1)t[i>>1]=apply(q[i>>1],combT(t[i&~1],t[i|1]),cw(i>>1));
  }
  void prp(int l, int r) {
    for (int i = h; i >= 1; i--) {
      prpd(r>>i);
      if ((r>>i)^(l>>i))prpd(l>>i);
    }
  }
  inline void ins(int i, Q v) {
    if(i<n) q[i]=combQ(v,q[i]);
    t[i]=apply(v,t[i],cw(i));
    prpu(i);
  }
  void upd(int l, int r, Q v) {
    l += n, r += n;
    prp(l,r);
    int el = l, er = r++;
    for (;l < r; l >>= 1, r >>= 1) {
      if (l & 1) ins(l++, v);
      if (r & 1) ins(--r, v);
    }
    for (;el|er;el>>=1,er >>=1) {
      prpu(el);
      if (el^er) prpu(er);
    }
  }
  T query(int l, int r) {
    l += n, r += n;
    prp(l,r);
    if (l==r) return t[l];
    T resl = t[l++], resr = t[r];
    for (; l < r; l >>= 1, r >>= 1) {
      if (l & 1) resl = combT(resl, t[l++]);
      if (r & 1) resr = combT(t[--r], resr);
    }
    return combT(resl, resr);
  }
};

// Example code for https://cses.fi/problemset/task/1735/
typedef tuple<int,L,bool> Q;
seg_tree<L, Q> seg(
  [](L a, L b){return a+b;},
  [](Q a, Q b) {
    if (b>a)swap(a,b);
    auto [t1,v1,op1] = a;
    auto [t2,v2,op2] = b;
    return (op1)?a:Q{t1,v1+v2,op2};
  },
  [](Q a, L b, int w){return get<0>(a)?get<2>(a)?get<1>(a)*w:b+get<1>(a)*w:b;}
);

int main() {
  int n, q;
  cin >> n >> q;
  vector<L> nums(n,0);
  for (int i = 0; i < n; i++) {
    cin >> nums[i];
  }
  seg.build(nums,{});
  for (int i = 1; i <= q; i++) {
    int t;
    cin >> t;
    if (t == 1) {
      int a, b; L x;
      cin >> a >> b >> x; a--; b--;
      seg.upd(a,b,{i,x,false});
    } else if (t == 2) {
      int a, b; L x;
      cin >> a >> b >> x; a--; b--;
      seg.upd(a,b,{i,x,true});
    } else if (t == 3) {
      int a, b;
      cin >> a >> b; a--; b--;
      cout << seg.query(a,b) << '\n';
    }
  }
}
