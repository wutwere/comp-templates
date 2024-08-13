double area_of_triangle(PT a, PT b, PT c) {
  return fabs(cross(b - a, c - a) * 0.5);
}
// -1 if strictly inside, 0 if on the polygon, 1 if strictly outside
int is_point_in_triangle(PT a, PT b, PT c, PT p) {
  if (sign(cross(b - a,c - a)) < 0) swap(b, c);
  int c1 = sign(cross(b - a,p - a));
  int c2 = sign(cross(c - b,p - b));
  int c3 = sign(cross(a - c,p - c));
  if (c1<0 || c2<0 || c3 < 0) return 1;
  if (c1 + c2 + c3 != 3) return 0;
  return -1;
}
double perimeter(vector<PT> &p) {
  double ans=0; int n = p.size();
  for (int i = 0; i < n; i++) ans += dist(p[i], p[(i + 1) % n]);
  return ans;
}
double area(vector<PT> &p) {
  double ans = 0; int n = p.size();
  for (int i = 0; i < n; i++) ans += cross(p[i], p[(i + 1) % n]);
  return fabs(ans) * 0.5;
}
// centroid of a (possibly non-convex) polygon, 
// assuming that the coordinates are listed in a clockwise or
// counterclockwise fashion.  Note that the centroid is often known as
// the "center of gravity" or "center of mass".
PT centroid(vector<PT> &p) {
  int n = p.size(); PT c(0, 0);
  double sum = 0;
  for (int i = 0; i < n; i++) sum += cross(p[i], p[(i + 1) % n]);
  double scale = 3.0 * sum;
  for (int i = 0; i < n; i++) {
    int j = (i + 1) % n;
    c = c + (p[i] + p[j]) * cross(p[i], p[j]);
  }
  return c / scale;
}

// 0 if cw, 1 if ccw
bool get_direction(vector<PT> &p) {
  double ans = 0; int n = p.size();
  for (int i = 0; i < n; i++) ans += cross(p[i], p[(i + 1) % n]);
  if (sign(ans) > 0) return 1;
  return 0;
}
// it returns a point such that the sum of distances
// from that point to all points in p  is minimum
// O(n log^2 MX)
PT geometric_median(vector<PT> p) {
  auto tot_dist = [&](PT z) {
    double res = 0;
    for (int i = 0; i < p.size(); i++) res += dist(p[i], z);
    return res;
  };
  auto findY = [&](double x) {
    double yl = -1e5, yr = 1e5;
    for (int i = 0; i < 60; i++) {
      double ym1 = yl + (yr - yl) / 3;
      double ym2 = yr - (yr - yl) / 3;
      double d1 = tot_dist(PT(x, ym1));
      double d2 = tot_dist(PT(x, ym2));
      if (d1 < d2) yr = ym2;
      else yl = ym1;
    }
    return pair<double, double> (yl, tot_dist(PT(x, yl)));
  };
  double xl = -1e5, xr = 1e5;
  for (int i = 0; i < 60; i++) {
    double xm1 = xl + (xr - xl) / 3;
    double xm2 = xr - (xr - xl) / 3;
    double y1, d1, y2, d2;
    auto z = findY(xm1); y1 = z.first; d1 = z.second;
    z = findY(xm2); y2 = z.first; d2 = z.second;
    if (d1 < d2) xr = xm2;
    else xl = xm1;
  }
  return {xl, findY(xl).first };
}
vector<PT> convex_hull(vector<PT> &p) {
  if (p.size() <= 1) return p;
  vector<PT> v = p;
  sort(v.begin(), v.end());
  vector<PT> up, dn;
  for (auto& p : v) {
    while (up.size() > 1 && orientation(up[up.size() - 2], up.back(), p) >= 0) {
      up.pop_back();
    }
    while (dn.size() > 1 && orientation(dn[dn.size() - 2], dn.back(), p) <= 0) {
      dn.pop_back();
    }
    up.push_back(p);
    dn.push_back(p);
  }
  v = dn;
  if (v.size() > 1) v.pop_back();
  reverse(up.begin(), up.end());
  up.pop_back();
  for (auto& p : up) {
    v.push_back(p);
  }
  if (v.size() == 2 && v[0] == v[1]) v.pop_back();
  return v;
}

//checks if convex or not
bool is_convex(vector<PT> &p) {
  bool s[3]; s[0] = s[1] = s[2] = 0;
  int n = p.size();
  for (int i = 0; i < n; i++) {
    int j = (i + 1) % n;
    int k = (j + 1) % n;
    s[sign(cross(p[j] - p[i], p[k] - p[i])) + 1] = 1;
    if (s[0] && s[2]) return 0;
  }
  return 1;
}
// -1 if strictly inside, 0 if on the polygon, 1 if strictly outside
// it must be strictly convex, otherwise make it strictly convex first
int is_point_in_convex(vector<PT> &p, const PT& x) { // O(log n)
  int n = p.size(); assert(n >= 3);
  int a = orientation(p[0], p[1], x), b = orientation(p[0], p[n - 1], x);
  if (a < 0 || b > 0) return 1;
  int l = 1, r = n - 1;
  while (l + 1 < r) {
    int mid = l + r >> 1;
    if (orientation(p[0], p[mid], x) >= 0) l = mid;
    else r = mid;
  }
  int k = orientation(p[l], p[r], x);
  if (k <= 0) return -k;
  if (l == 1 && a == 0) return 0;
  if (r == n - 1 && b == 0) return 0;
  return -1;
}
bool is_point_on_polygon(vector<PT> &p, const PT& z) {
  int n = p.size();
  for (int i = 0; i < n; i++) {
    if (is_point_on_seg(p[i], p[(i + 1) % n], z)) return 1;
  }
  return 0;
}
// returns 1e9 if the point is on the polygon 
int winding_number(vector<PT> &p, const PT& z) { // O(n)
  if (is_point_on_polygon(p, z)) return 1e9;
  int n = p.size(), ans = 0;
  for (int i = 0; i < n; ++i) {
    int j = (i + 1) % n;
    bool below = p[i].y < z.y;
    if (below != (p[j].y < z.y)) {
      auto orient = orientation(z, p[j], p[i]);
      if (orient == 0) return 0;
      if (below == (orient > 0)) ans += below ? 1 : -1;
    }
  }
  return ans;
}
// -1 if strictly inside, 0 if on the polygon, 1 if strictly outside
int is_point_in_polygon(vector<PT> &p, const PT& z) { // O(n)
  int k = winding_number(p, z);
  return k == 1e9 ? 0 : k == 0 ? 1 : -1;
}

// id of the vertex having maximum dot product with z
// polygon must need to be convex
// top - upper right vertex
// for minimum dot product negate z and return -dot(z, p[id])
int extreme_vertex(vector<PT> &p, const PT &z, const int top) { // O(log n)
  int n = p.size();
  if (n == 1) return 0;
  double ans = dot(p[0], z); int id = 0;
  if (dot(p[top], z) > ans) ans = dot(p[top], z), id = top;
  int l = 1, r = top - 1;
  while (l < r) {
    int mid = l + r >> 1;
    if (dot(p[mid + 1], z) >= dot(p[mid], z)) l = mid + 1;
    else r = mid;
  }
  if (dot(p[l], z) > ans) ans = dot(p[l], z), id = l;
  l = top + 1, r = n - 1;
  while (l < r) {
    int mid = l + r >> 1;
    if (dot(p[(mid + 1) % n], z) >= dot(p[mid], z)) l = mid + 1;
    else r = mid;
  }
  l %= n;
  if (dot(p[l], z) > ans) ans = dot(p[l], z), id = l;
  return id;
}
// maximum distance from any point on the perimeter to another point on the perimeter
double diameter(vector<PT> &p) {
  int n = (int)p.size();
  if (n == 1) return 0;
  if (n == 2) return dist(p[0], p[1]);
  double ans = 0;
  int i = 0, j = 1;
  while (i < n) {
    while (cross(p[(i + 1) % n] - p[i], p[(j + 1) % n] - p[j]) >= 0) {
      ans = max(ans, dist2(p[i], p[j]));
      j = (j + 1) % n;
    }
    ans = max(ans, dist2(p[i], p[j]));
    i++;
  }
  return sqrt(ans);
}
// minimum distance between two parallel lines (non necessarily axis parallel)
// such that the polygon can be put between the lines
double width(vector<PT> &p) {
  int n = (int)p.size();
  if (n <= 2) return 0;
  double ans = inf;
  int i = 0, j = 1;
  while (i < n) {
    while (cross(p[(i + 1) % n] - p[i], p[(j + 1) % n] - p[j]) >= 0) j = (j + 1) % n;
    ans = min(ans, dist_from_point_to_line(p[i], p[(i + 1) % n], p[j]));
    i++;
  }
  return ans;
}

// minimum perimeter
double minimum_enclosing_rectangle(vector<PT> &p) {
  int n = p.size();
  if (n <= 2) return perimeter(p);
  int mndot = 0; double tmp = dot(p[1] - p[0], p[0]);
  for (int i = 1; i < n; i++) {
    if (dot(p[1] - p[0], p[i]) <= tmp) {
      tmp = dot(p[1] - p[0], p[i]);
      mndot = i;
    }
  }
  double ans = inf;
  int i = 0, j = 1, mxdot = 1;
  while (i < n) {
    PT cur = p[(i + 1) % n] - p[i];
    while (cross(cur, p[(j + 1) % n] - p[j]) >= 0) j = (j + 1) % n;
    while (dot(p[(mxdot + 1) % n], cur) >= dot(p[mxdot], cur)) mxdot = (mxdot + 1) % n;
    while (dot(p[(mndot + 1) % n], cur) <= dot(p[mndot], cur)) mndot = (mndot + 1) % n;
    ans = min(ans, 2.0 * ((dot(p[mxdot], cur) / cur.norm() - dot(p[mndot], cur) / cur.norm()) + dist_from_point_to_line(p[i], p[(i + 1) % n], p[j])));
    i++;
  }
  return ans;
}
// given n points, find the minimum enclosing circle of the points
// call convex_hull() before this for faster solution
// expected O(n)
circle minimum_enclosing_circle(vector<PT> &p) {
  random_shuffle(p.begin(), p.end());
  int n = p.size();
  circle c(p[0], 0);
  for (int i = 1; i < n; i++) {
    if (sign(dist(c.p, p[i]) - c.r) > 0) {
      c = circle(p[i], 0);
      for (int j = 0; j < i; j++) {
        if (sign(dist(c.p, p[j]) - c.r) > 0) {
          c = circle((p[i] + p[j]) / 2, dist(p[i], p[j]) / 2);
          for (int k = 0; k < j; k++) {
            if (sign(dist(c.p, p[k]) - c.r) > 0) {
              c = circle(p[i], p[j], p[k]);
            }
          }
        }
      }
    }
  }
  return c;
}
// returns a vector with the vertices of a polygon with everything 
// to the left of the line going from a to b cut away.
vector<PT> cut(vector<PT> &p, PT a, PT b) {
  vector<PT> ans;
  int n = (int)p.size();
  for (int i = 0; i < n; i++) {
    double c1 = cross(b - a, p[i] - a);
    double c2 = cross(b - a, p[(i + 1) % n] - a);
    if (sign(c1) >= 0) ans.push_back(p[i]);
    if (sign(c1 * c2) < 0) {
      if (!is_parallel(p[i], p[(i + 1) % n], a, b)) {
        PT tmp; line_line_intersection(p[i], p[(i + 1) % n], a, b, tmp);
        ans.push_back(tmp);
      }
    }
  }
  return ans;
}

// not necessarily convex, boundary is included in the intersection
// returns total intersected length
// it returns the sum of the lengths of the portions of the line that are inside the polygon
double polygon_line_intersection(vector<PT> p, PT a, PT b) {
  int n = p.size();
  p.push_back(p[0]);
  line l = line(a, b);
  double ans = 0.0;
  vector< pair<double, int> > vec;
  for (int i = 0; i < n; i++) {
    int s1 = orientation(a, b, p[i]);
    int s2 = orientation(a, b, p[i + 1]);
    if (s1 == s2) continue;
    line t = line(p[i], p[i + 1]);
    PT inter = (t.v * l.c - l.v * t.c) / cross(l.v, t.v);
    double tmp = dot(inter, l.v);
    int f;
    if (s1 > s2) f = s1 && s2 ? 2 : 1;
    else f = s1 && s2 ? -2 : -1;
    vec.push_back(make_pair((f > 0 ? tmp - eps : tmp + eps), f)); // keep eps very small like 1e-12
  }
  sort(vec.begin(), vec.end());
  for (int i = 0, j = 0; i + 1 < (int)vec.size(); i++){
    j += vec[i].second;
    if (j) ans += vec[i + 1].first - vec[i].first; // if this portion is inside the polygon
    // else ans = 0; // if we want the maximum intersected length which is totally inside the polygon, uncomment this and take the maximum of ans
  }
  ans = ans / sqrt(dot(l.v, l.v));
  p.pop_back();
  return ans;
}
// given a convex polygon p, and a line ab and the top vertex of the polygon
// returns the intersection of the line with the polygon
// it returns the indices of the edges of the polygon that are intersected by the line
// so if it returns i, then the line intersects the edge (p[i], p[(i + 1) % n])
array<int, 2> convex_line_intersection(vector<PT> &p, PT a, PT b, int top) {
  int end_a = extreme_vertex(p, (a - b).perp(), top);
  int end_b = extreme_vertex(p, (b - a).perp(), top);
  auto cmp_l = [&](int i) { return orientation(a, p[i], b); };
  if (cmp_l(end_a) < 0 || cmp_l(end_b) > 0)
    return {-1, -1}; // no intersection
  array<int, 2> res;
  for (int i = 0; i < 2; i++) {
    int lo = end_b, hi = end_a, n = p.size();
    while ((lo + 1) % n != hi) {
      int m = ((lo + hi + (lo < hi ? 0 : n)) / 2) % n;
      (cmp_l(m) == cmp_l(end_b) ? lo : hi) = m;
    }
    res[i] = (lo + !cmp_l(hi)) % n;
    swap(end_a, end_b);
  }
  if (res[0] == res[1]) return {res[0], -1}; // touches the vertex res[0]
  if (!cmp_l(res[0]) && !cmp_l(res[1])) 
    switch ((res[0] - res[1] + (int)p.size() + 1) % p.size()) {
      case 0: return {res[0], res[0]}; // touches the edge (res[0], res[0] + 1)
      case 2: return {res[1], res[1]}; // touches the edge (res[1], res[1] + 1)
    }
  return res; // intersects the edges (res[0], res[0] + 1) and (res[1], res[1] + 1)
}

pair<PT, int> point_poly_tangent(vector<PT> &p, PT Q, int dir, int l, int r) {
  while (r - l > 1) {
    int mid = (l + r) >> 1;
    bool pvs = orientation(Q, p[mid], p[mid - 1]) != -dir;
    bool nxt = orientation(Q, p[mid], p[mid + 1]) != -dir;
    if (pvs && nxt) return {p[mid], mid};
    if (!(pvs || nxt)) {
      auto p1 = point_poly_tangent(p, Q, dir, mid + 1, r);
      auto p2 = point_poly_tangent(p, Q, dir, l, mid - 1);
      return orientation(Q, p1.first, p2.first) == dir ? p1 : p2;
    }
    if (!pvs) {
      if (orientation(Q, p[mid], p[l]) == dir)  r = mid - 1;
      else if (orientation(Q, p[l], p[r]) == dir) r = mid - 1;
      else l = mid + 1;
    }
    if (!nxt) {
      if (orientation(Q, p[mid], p[l]) == dir)  l = mid + 1;
      else if (orientation(Q, p[l], p[r]) == dir) r = mid - 1;
      else l = mid + 1;
    }
  }
  pair<PT, int> ret = {p[l], l};
  for (int i = l + 1; i <= r; i++) ret = orientation(Q, ret.first, p[i]) != dir ? make_pair(p[i], i) : ret;
  return ret;
}
// (ccw, cw) tangents from a point that is outside this convex polygon
// returns indexes of the points
// ccw means the tangent from Q to that point is in the same direction as the polygon ccw direction
pair<int, int> tangents_from_point_to_polygon(vector<PT> &p, PT Q){
  int ccw = point_poly_tangent(p, Q, 1, 0, (int)p.size() - 1).second;
  int cw = point_poly_tangent(p, Q, -1, 0, (int)p.size() - 1).second;
  return make_pair(ccw, cw);
}

// minimum distance from a point to a convex polygon
// it assumes point lie strictly outside the polygon
double dist_from_point_to_polygon(vector<PT> &p, PT z) {
  double ans = inf;
  int n = p.size();
  if (n <= 3) {
    for(int i = 0; i < n; i++) ans = min(ans, dist_from_point_to_seg(p[i], p[(i + 1) % n], z));
    return ans;
  }
  auto [r, l] = tangents_from_point_to_polygon(p, z);
  if(l > r) r += n;
  while (l < r) {
    int mid = (l + r) >> 1;
    double left = dist2(p[mid % n], z), right= dist2(p[(mid + 1) % n], z);
    ans = min({ans, left, right});
    if(left < right) r = mid;
    else l = mid + 1;
  }
  ans = sqrt(ans);
  ans = min(ans, dist_from_point_to_seg(p[l % n], p[(l + 1) % n], z));
  ans = min(ans, dist_from_point_to_seg(p[l % n], p[(l - 1 + n) % n], z));
  return ans;
}
// minimum distance from convex polygon p to line ab
// returns 0 is it intersects with the polygon
// top - upper right vertex
double dist_from_polygon_to_line(vector<PT> &p, PT a, PT b, int top) { //O(log n)
  PT orth = (b - a).perp();
  if (orientation(a, b, p[0]) > 0) orth = (a - b).perp();
  int id = extreme_vertex(p, orth, top);
  if (dot(p[id] - a, orth) > 0) return 0.0; //if orth and a are in the same half of the line, then poly and line intersects
  return dist_from_point_to_line(a, b, p[id]); //does not intersect
}

// calculates the area of the union of n polygons (not necessarily convex). 
// the points within each polygon must be given in CCW order.
// complexity: O(N^2), where N is the total number of points
double rat(PT a, PT b, PT p) {
  return !sign(a.x - b.x) ? (p.y - a.y) / (b.y - a.y) : (p.x - a.x) / (b.x - a.x);
};
double polygon_union(vector<vector<PT>> &p) {
  int n = p.size();
  double ans=0;
  for(int i = 0; i < n; ++i) {
    for (int v = 0; v < (int)p[i].size(); ++v) {
      PT a = p[i][v], b = p[i][(v + 1) % p[i].size()];
      vector<pair<double, int>> segs;
      segs.emplace_back(0,  0), segs.emplace_back(1,  0);
      for(int j = 0; j < n; ++j) {
        if(i != j) {
          for(size_t u = 0; u < p[j].size(); ++u) {
            PT c = p[j][u], d = p[j][(u + 1) % p[j].size()];
            int sc = sign(cross(b - a, c - a)), sd = sign(cross(b - a, d - a));
            if(!sc && !sd) {
              if(sign(dot(b - a, d - c)) > 0 && i > j) {
                segs.emplace_back(rat(a, b, c), 1), segs.emplace_back(rat(a, b, d),  -1);
              }
            } 
            else {
              double sa = cross(d - c, a - c), sb = cross(d - c, b - c);
              if(sc >= 0 && sd < 0) segs.emplace_back(sa / (sa - sb), 1);
              else if(sc < 0 && sd >= 0) segs.emplace_back(sa / (sa - sb),  -1);
            }
          }
        }
      }
      sort(segs.begin(),  segs.end());
      double pre = min(max(segs[0].first, 0.0), 1.0), now, sum = 0;
      int cnt = segs[0].second;
      for(int j = 1; j < segs.size(); ++j) {
        now = min(max(segs[j].first, 0.0), 1.0);
        if (!cnt) sum += now - pre;
        cnt += segs[j].second;
        pre = now;
      }
      ans += cross(a, b) * sum;
    }
  }
  return ans * 0.5;
}

// minimum distance from a convex polygon to another convex polygon
// the polygon doesnot overlap or touch
// tested in https://toph.co/p/the-wall
double dist_from_polygon_to_polygon(vector<PT> &p1, vector<PT> &p2) { // O(n log n)
  double ans = inf;
  for (int i = 0; i < p1.size(); i++) {
    ans = min(ans, dist_from_point_to_polygon(p2, p1[i]));
  }
  for (int i = 0; i < p2.size(); i++) {
    ans = min(ans, dist_from_point_to_polygon(p1, p2[i]));
  }
  return ans;
}

// maximum distance from a convex polygon to another convex polygon
double maximum_dist_from_polygon_to_polygon(vector<PT> &u, vector<PT> &v){ //O(n)
  int n = (int)u.size(), m = (int)v.size();
  double ans = 0;
  if (n < 3 || m < 3) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) ans = max(ans, dist2(u[i], v[j]));
    }
    return sqrt(ans);
  }
  if (u[0].x > v[0].x) swap(n, m), swap(u, v);
  int i = 0, j = 0, step = n + m + 10;
  while (j + 1 < m && v[j].x < v[j + 1].x) j++ ;
  while (step--) {
    if (cross(u[(i + 1)%n] - u[i], v[(j + 1)%m] - v[j]) >= 0) j = (j + 1) % m;
    else i = (i + 1) % n;
    ans = max(ans, dist2(u[i], v[j]));
  }
  return sqrt(ans);
}

// contains all points p such that: cross(b - a, p - a) >= 0
struct HP {
  PT a, b;
  HP() {}
  HP(PT a, PT b) : a(a), b(b) {}
  HP(const HP& rhs) : a(rhs.a), b(rhs.b) {}
  int operator < (const HP& rhs) const {
    PT p = b - a;
    PT q = rhs.b - rhs.a;
    int fp = (p.y < 0 || (p.y == 0 && p.x < 0));
    int fq = (q.y < 0 || (q.y == 0 && q.x < 0));
    if (fp != fq) return fp == 0;
    if (cross(p, q)) return cross(p, q) > 0;
    return cross(p, rhs.b - a) < 0;
  }
  PT line_line_intersection(PT a, PT b, PT c, PT d) {
    b = b - a; d = c - d; c = c - a;
    return a + b * cross(c, d) / cross(b, d);
  }
  PT intersection(const HP &v) {
    return line_line_intersection(a, b, v.a, v.b);
  }
};
int check(HP a, HP b, HP c) {
  return cross(a.b - a.a, b.intersection(c) - a.a) > -eps; //-eps to include polygons of zero area (straight lines, points)
}

// consider half-plane of counter-clockwise side of each line
// if lines are not bounded add infinity rectangle
// returns a convex polygon, a point can occur multiple times though
// complexity: O(n log(n))
vector<PT> half_plane_intersection(vector<HP> h) {
  sort(h.begin(), h.end());
  vector<HP> tmp;
  for (int i = 0; i < h.size(); i++) {
    if (!i || cross(h[i].b - h[i].a, h[i - 1].b - h[i - 1].a)) {
      tmp.push_back(h[i]);
    }
  }
  h = tmp;
  vector<HP> q(h.size() + 10);
  int qh = 0, qe = 0;
  for (int i = 0; i < h.size(); i++) {
    while (qe - qh > 1 && !check(h[i], q[qe - 2], q[qe - 1])) qe--;
    while (qe - qh > 1 && !check(h[i], q[qh], q[qh + 1])) qh++;
    q[qe++] = h[i];
  }
  while (qe - qh > 2 && !check(q[qh], q[qe - 2], q[qe - 1])) qe--;
  while (qe - qh > 2 && !check(q[qe - 1], q[qh], q[qh + 1])) qh++;
  vector<HP> res; 
  for (int i = qh; i < qe; i++) res.push_back(q[i]);
  vector<PT> hull;
  if (res.size() > 2) {
    for (int i = 0; i < res.size(); i++) {
      hull.push_back(res[i].intersection(res[(i + 1) % ((int)res.size())]));
    }
  }
  return hull;
}
// rotate the polygon such that the (bottom, left)-most point is at the first position
void reorder_polygon(vector<PT> &p) {
  int pos = 0;
  for (int i = 1; i < p.size(); i++) {
    if (p[i].y < p[pos].y || (sign(p[i].y - p[pos].y) == 0 && p[i].x < p[pos].x)) pos = i;
  }
  rotate(p.begin(), p.begin() + pos, p.end());
}
// a and b are convex polygons
// returns a convex hull of their minkowski sum
// min(a.size(), b.size()) >= 2
// https://cp-algorithms.com/geometry/minkowski.html
vector<PT> minkowski_sum(vector<PT> a, vector<PT> b) {
  reorder_polygon(a); reorder_polygon(b);
  int n = a.size(), m = b.size();
  int i = 0, j = 0;
  a.push_back(a[0]); a.push_back(a[1]);
  b.push_back(b[0]); b.push_back(b[1]);
  vector<PT> c;
  while (i < n || j < m) {
    c.push_back(a[i] + b[j]);
    double p = cross(a[i + 1] - a[i], b[j + 1] - b[j]);
    if (sign(p) >= 0) ++i;
    if (sign(p) <= 0) ++j;
  }
  return c;
}
// returns the area of the intersection of the circle with center c and radius r
// and the triangle formed by the points c, a, b
double _triangle_circle_intersection(PT c, double r, PT a, PT b) {
  double sd1 = dist2(c, a), sd2 = dist2(c, b);
  if(sd1 > sd2) swap(a, b), swap(sd1, sd2);
  double sd = dist2(a, b);
  double d1 = sqrtl(sd1), d2 = sqrtl(sd2), d = sqrt(sd);
  double x = abs(sd2 - sd - sd1) / (2 * d);
  double h = sqrtl(sd1 - x * x);
  if(r >= d2) return h * d / 2;
  double area = 0;
  if(sd + sd1 < sd2) {
    if(r < d1) area = r * r * (acos(h / d2) - acos(h / d1)) / 2;
    else {
      area = r * r * ( acos(h / d2) - acos(h / r)) / 2;
      double y = sqrtl(r * r - h * h);
      area += h * (y - x) / 2;
    }
  } 
  else {
    if(r < h) area = r * r * (acos(h / d2) + acos(h / d1)) / 2;
    else {
      area += r * r * (acos(h / d2) - acos(h / r)) / 2;
      double y = sqrtl(r * r - h * h);
      area += h * y / 2;
      if(r < d1) {
        area += r * r * (acos(h / d1) - acos(h / r)) / 2;
        area += h * y / 2;
      } 
      else area += h * x / 2;
    }
  }
  return area;
}

// intersection between a simple polygon and a circle
double polygon_circle_intersection(vector<PT> &v, PT p, double r) {
  int n = v.size();
  double ans = 0.00;
  PT org = {0, 0};
  for(int i = 0; i < n; i++) {
    int x = orientation(p, v[i], v[(i + 1) % n]);
    if(x == 0) continue;
    double area = _triangle_circle_intersection(org, r, v[i] - p, v[(i + 1) % n] - p);
    if (x < 0) ans -= area;
    else ans += area;
  }
  return abs(ans);
}
// find a circle of radius r that contains as many points as possible
// O(n^2 log n);
double maximum_circle_cover(vector<PT> p, double r, circle &c) {
  int n = p.size();
  int ans = 0;
  int id = 0; double th = 0;
  for (int i = 0; i < n; ++i) {
    // maximum circle cover when the circle goes through this point
    vector<pair<double, int>> events = {{-PI, +1}, {PI, -1}};
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      double d = dist(p[i], p[j]);
      if (d > r * 2) continue;
      double dir = (p[j] - p[i]).arg();
      double ang = acos(d / 2 / r);
      double st = dir - ang, ed = dir + ang;
      if (st > PI) st -= PI * 2;
      if (st <= -PI) st += PI * 2;
      if (ed > PI) ed -= PI * 2;
      if (ed <= -PI) ed += PI * 2;
      events.push_back({st - eps, +1}); // take care of precisions!
      events.push_back({ed, -1});
      if (st > ed) {
        events.push_back({-PI, +1});
        events.push_back({+PI, -1});
      }
    }
    sort(events.begin(), events.end());
    int cnt = 0;
    for (auto &&e: events) {
      cnt += e.second;
      if (cnt > ans) {
        ans = cnt;
        id = i; th = e.first;
      }
    }
  }
  PT w = PT(p[id].x + r * cos(th), p[id].y + r * sin(th));
  c = circle(w, r); //best_circle
  return ans;
}

// radius of the maximum inscribed circle in a convex polygon
double maximum_inscribed_circle(vector<PT> p) {
  int n = p.size();
  if (n <= 2) return 0;
  double l = 0, r = 20000;
  while (r - l > eps) {
    double mid = (l + r) * 0.5;
    vector<HP> h;
    const int L = 1e9;
    h.push_back(HP(PT(-L, -L), PT(L, -L)));
    h.push_back(HP(PT(L, -L), PT(L, L)));
    h.push_back(HP(PT(L, L), PT(-L, L)));
    h.push_back(HP(PT(-L, L), PT(-L, -L)));
    for (int i = 0; i < n; i++) {
      PT z = (p[(i + 1) % n] - p[i]).perp();
      z = z.truncate(mid);
      PT y = p[i] + z, q = p[(i + 1) % n] + z;
      h.push_back(HP(p[i] + z, p[(i + 1) % n] + z));
    }
    vector<PT> nw = half_plane_intersection(h);
    if (!nw.empty()) l = mid;
    else r = mid;
  }
  return l;
}
// ear decomposition, O(n^3) but faster
vector<vector<PT>> triangulate(vector<PT> p) {
  vector<vector<PT>> v;
  while (p.size() >= 3) {
    for (int i = 0, n = p.size(); i < n; i++) {
      int pre = i == 0 ? n - 1 : i - 1;;
      int nxt = i == n - 1 ? 0 : i + 1;;
      int ori = orientation(p[i], p[pre], p[nxt]);
      if (ori < 0) {
        int ok = 1;
        for (int j = 0; j < n; j++) {
          if (j == i || j == pre || j == nxt)continue;
          if (is_point_in_triangle(p[i], p[pre], p[nxt] , p[j]) < 1) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          v.push_back({p[pre], p[i], p[nxt]});
          p.erase(p.begin() + i);
          break;
        }
      }
    }
  }
  return v;
}

struct star {
  int n;    // number of sides of the star
  double r; // radius of the circumcircle
  star(int _n, double _r) {
    n = _n;
    r = _r;
  }

  double area() {
    double theta = PI / n;
    double s = 2 * r * sin(theta);
    double R = 0.5 * s / tan(theta);
    double a = 0.5 * n * s * R;
    double a2 = 0.25 * s * s / tan(1.5 * theta);
    return a - n * a2;
  }
};

// given a list of lengths of the sides of a polygon in counterclockwise order
// returns the maximum area of a non-degenerate polygon that can be formed using those lengths
double get_maximum_polygon_area_for_given_lengths(vector<double> v) {
  if (v.size() < 3) {
    return 0;
  }
  int m = 0;
  double sum = 0;
  for (int i = 0; i < v.size(); i++) {
    if (v[i] > v[m]) {
      m = i;
    }
    sum += v[i];
  }
  if (sign(v[m] - (sum - v[m])) >= 0) {
    return 0; // no non-degenerate polygon is possible
  }
  // the polygon should be a circular polygon
  // that is all points are on the circumference of a circle
  double l = v[m] / 2, r = 1e6; // fix it correctly
  int it = 60;
  auto ang = [](double x, double r) { // x = length of the chord, r = radius of the circle
    return 2 * asin((x / 2) / r);
  };
  auto calc = [=](double r) {
    double sum = 0;
    for (auto x: v) {
      sum += ang(x, r);
    }
    return sum;
  };
  // compute the radius of the circle
  while (it--) {
    double mid = (l + r) / 2;
    if (calc(mid) <= 2 * PI) {
      r = mid;
    }
    else {
      l = mid;
    }
  }

  if (calc(r) <= 2 * PI - eps) { // the center of the circle is outside the polygon
    auto calc2 = [&](double r) {
      double sum = 0;
      for (int i = 0; i < v.size(); i++) {
        double x = v[i];
        double th = ang(x, r);
        if (i != m) {
          sum += th;
        }
        else {
          sum += 2 * PI - th;
        }
      }
      return sum;
    };
    l = v[m] / 2; r = 1e6;
    it = 60;
    while (it--) {
      double mid = (l + r) / 2;
      if (calc2(mid) > 2 * PI) {
        r = mid;
      }
      else {
        l = mid;
      }
    }
    auto get_area = [=](double r) {
      double ans = 0;
      for (int i = 0; i < v.size(); i++) {
        double x = v[i];
        double area = r * r * sin(ang(x, r)) / 2;
        if (i != m) {
          ans += area;
        }
        else {
          ans -= area;
        }
      }
      return ans;
    };
    return get_area(r);
  }
  else { // the center of the circle is inside the polygon
    auto get_area = [=](double r) {
      double ans = 0;
      for (auto x: v) {
        ans += r * r * sin(ang(x, r)) / 2;
      }
      return ans;
    };
    return get_area(r);
  }
}
