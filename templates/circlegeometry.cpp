struct circle {
  PT p; double r;
  circle() {}
  circle(PT _p, double _r): p(_p), r(_r) {};
  // center (x, y) and radius r
  circle(double x, double y, double _r): p(PT(x, y)), r(_r) {};
  // circumcircle of a triangle
  // the three points must be unique
  circle(PT a, PT b, PT c) {
    b = (a + b) * 0.5;
    c = (a + c) * 0.5;
    line_line_intersection(b, b + rotatecw90(a - b), c, c + rotatecw90(a - c), p);
    r = dist(a, p);
  }
  // inscribed circle of a triangle
  // pass a bool just to differentiate from circumcircle
  circle(PT a, PT b, PT c, bool t) {
    line u, v;
    double m = atan2(b.y - a.y, b.x - a.x), n = atan2(c.y - a.y, c.x - a.x);
    u.a = a;
    u.b = u.a + (PT(cos((n + m)/2.0), sin((n + m)/2.0)));
    v.a = b;
    m = atan2(a.y - b.y, a.x - b.x), n = atan2(c.y - b.y, c.x - b.x);
    v.b = v.a + (PT(cos((n + m)/2.0), sin((n + m)/2.0)));
    line_line_intersection(u.a, u.b, v.a, v.b, p);
    r = dist_from_point_to_seg(a, b, p);
  }
  bool operator == (circle v) { return p == v.p && sign(r - v.r) == 0; }
  double area() { return PI * r * r; }
  double circumference() { return 2.0 * PI * r; }
};
//0 if outside, 1 if on circumference, 2 if inside circle
int circle_point_relation(PT p, double r, PT b) {
  double d = dist(p, b);
  if (sign(d - r) < 0) return 2;
  if (sign(d - r) == 0) return 1;
  return 0;
}
// 0 if outside, 1 if on circumference, 2 if inside circle
int circle_line_relation(PT p, double r, PT a, PT b) {
  double d = dist_from_point_to_line(a, b, p);
  if (sign(d - r) < 0) return 2;
  if (sign(d - r) == 0) return 1;
  return 0;
}
//compute intersection of line through points a and b with
//circle centered at c with radius r > 0
vector<PT> circle_line_intersection(PT c, double r, PT a, PT b) {
  vector<PT> ret;
  b = b - a; a = a - c;
  double A = dot(b, b), B = dot(a, b);
  double C = dot(a, a) - r * r, D = B * B - A * C;
  if (D < -eps) return ret;
  ret.push_back(c + a + b * (-B + sqrt(D + eps)) / A);
  if (D > eps) ret.push_back(c + a + b * (-B - sqrt(D)) / A);
  return ret;
}


//5 - outside and do not intersect
//4 - intersect outside in one point
//3 - intersect in 2 points
//2 - intersect inside in one point
//1 - inside and do not intersect
int circle_circle_relation(PT a, double r, PT b, double R) {
  double d = dist(a, b);
  if (sign(d - r - R) > 0)  return 5;
  if (sign(d - r - R) == 0) return 4;
  double l = fabs(r - R);
  if (sign(d - r - R) < 0 && sign(d - l) > 0) return 3;
  if (sign(d - l) == 0) return 2;
  if (sign(d - l) < 0) return 1;
  assert(0); return -1;
}
vector<PT> circle_circle_intersection(PT a, double r, PT b, double R) {
  if (a == b && sign(r - R) == 0) return {PT(1e18, 1e18)};
  vector<PT> ret;
  double d = sqrt(dist2(a,  b));
  if (d > r + R || d + min(r,  R) < max(r,  R)) return ret;
  double x = (d * d - R * R + r * r) / (2 * d);
  double y = sqrt(r * r - x * x);
  PT v = (b - a) / d;
  ret.push_back(a + v * x  +  rotateccw90(v) * y);
  if (y > 0) ret.push_back(a + v * x - rotateccw90(v) * y);
  return ret;
}
// returns two circle c1, c2 through points a, b and of radius r
// 0 if there is no such circle, 1 if one circle, 2 if two circle
int get_circle(PT a, PT b, double r, circle &c1, circle &c2) {
  vector<PT> v = circle_circle_intersection(a, r, b, r);
  int t = v.size();
  if (!t) return 0;
  c1.p = v[0], c1.r = r;
  if (t == 2) c2.p = v[1], c2.r = r;
  return t;
}
// returns two circle c1, c2 which is tangent to line u,  goes through
// point q and has radius r1; 0 for no circle, 1 if c1 = c2 , 2 if c1 != c2
int get_circle(line u, PT q, double r1, circle &c1, circle &c2) {
  double d = dist_from_point_to_line(u.a, u.b, q);
  if (sign(d - r1 * 2.0) > 0) return 0;
  if (sign(d) == 0) {
    cout << u.v.x << ' ' << u.v.y << '\n';
    c1.p = q + rotateccw90(u.v).truncate(r1);
    c2.p = q + rotatecw90(u.v).truncate(r1);
    c1.r = c2.r = r1;
    return 2;
  }
  line u1 = line(u.a + rotateccw90(u.v).truncate(r1), u.b + rotateccw90(u.v).truncate(r1));
  line u2 = line(u.a + rotatecw90(u.v).truncate(r1), u.b + rotatecw90(u.v).truncate(r1));
  circle cc = circle(q, r1);
  PT p1, p2; vector<PT> v;
  v = circle_line_intersection(q, r1, u1.a, u1.b);
  if (!v.size()) v = circle_line_intersection(q, r1, u2.a, u2.b);
  v.push_back(v[0]);
  p1 = v[0], p2 = v[1];
  c1 = circle(p1, r1);
  if (p1 == p2) {
    c2 = c1;
    return 1;
  }
  c2 = circle(p2, r1);
  return 2;
}

// returns the circle such that for all points w on the circumference of the circle
// dist(w, a) : dist(w, b) = rp : rq
// rp != rq
// https://en.wikipedia.org/wiki/Circles_of_Apollonius
circle get_apollonius_circle(PT p, PT q, double rp, double rq ){
  rq *= rq ;
  rp *= rp ;
  double a = rq - rp ;
  assert(sign(a));
  double g = rq * p.x - rp * q.x ; g /= a ;
  double h = rq * p.y - rp * q.y ; h /= a ;
  double c = rq * p.x * p.x - rp * q.x * q.x + rq * p.y * p.y - rp * q.y * q.y ;
  c /= a ;
  PT o(g, h);
  double r = g * g + h * h - c ;
  r = sqrt(r);
  return circle(o,r);
}
// returns area of intersection between two circles
double circle_circle_area(PT a, double r1, PT b, double r2) {
  double d = (a - b).norm();
  if(r1 + r2 < d + eps) return 0;
  if(r1 + d < r2 + eps) return PI * r1 * r1;
  if(r2 + d < r1 + eps) return PI * r2 * r2;
  double theta_1 = acos((r1 * r1 + d * d - r2 * r2) / (2 * r1 * d)), 
         theta_2 = acos((r2 * r2 + d * d - r1 * r1)/(2 * r2 * d));
  return r1 * r1 * (theta_1 - sin(2 * theta_1)/2.) + r2 * r2 * (theta_2 - sin(2 * theta_2)/2.);
}
// tangent lines from point q to the circle
int tangent_lines_from_point(PT p, double r, PT q, line &u, line &v) {
  int x = sign(dist2(p, q) - r * r);
  if (x < 0) return 0; // point in cricle
  if (x == 0) { // point on circle
    u = line(q, q + rotateccw90(q - p));
    v = u;
    return 1;
  }
  double d = dist(p, q);
  double l = r * r / d;
  double h = sqrt(r * r - l * l);
  u = line(q, p + ((q - p).truncate(l) + (rotateccw90(q - p).truncate(h))));
  v = line(q, p + ((q - p).truncate(l) + (rotatecw90(q - p).truncate(h))));
  return 2;
}
// returns outer tangents line of two circles
// if inner == 1 it returns inner tangent lines
int tangents_lines_from_circle(PT c1, double r1, PT c2, double r2, bool inner, line &u, line &v) {
  if (inner) r2 = -r2;
  PT d = c2 - c1;
  double dr = r1 - r2, d2 = d.norm2(), h2 = d2 - dr * dr;
  if (d2 == 0 || h2 < 0) {
    assert(h2 != 0);
    return 0;
  }
  vector<pair<PT, PT>>out;
  for (int tmp: {- 1, 1}) {
    PT v = (d * dr + rotateccw90(d) * sqrt(h2) * tmp) / d2;
    out.push_back({c1 + v * r1, c2 + v * r2});
  }
  u = line(out[0].first, out[0].second);
  if (out.size() == 2) v = line(out[1].first, out[1].second);
  return 1 + (h2 > 0);
}


// O(n^2 log n)
// https://vjudge.net/problem/UVA-12056
struct CircleUnion {
  int n;
  double x[2020], y[2020], r[2020];
  int covered[2020];
  vector<pair<double, double> > seg, cover;
  double arc, pol;
  inline int sign(double x) {return x < -eps ? -1 : x > eps;}
  inline int sign(double x, double y) {return sign(x - y);}
  inline double SQ(const double x) {return x * x;}
  inline double dist(double x1, double y1, double x2, double y2) {return sqrt(SQ(x1 - x2) + SQ(y1 - y2));}
  inline double angle(double A, double B, double C) {
    double val = (SQ(A) + SQ(B) - SQ(C)) / (2 * A * B);
    if (val < -1) val = -1;
    if (val > +1) val = +1;
    return acos(val);
  }
  CircleUnion() {
    n = 0;
    seg.clear(), cover.clear();
    arc = pol = 0;
  }
  void init() {
    n = 0;
    seg.clear(), cover.clear();
    arc = pol = 0;
  }
  void add(double xx, double yy, double rr) {
    x[n] = xx, y[n] = yy, r[n] = rr, covered[n] = 0, n++;
  }
  void getarea(int i, double lef, double rig) {
    arc += 0.5 * r[i] * r[i] * (rig - lef - sin(rig - lef));
    double x1 = x[i] + r[i] * cos(lef), y1 = y[i] + r[i] * sin(lef);
    double x2 = x[i] + r[i] * cos(rig), y2 = y[i] + r[i] * sin(rig);
    pol += x1 * y2 - x2 * y1;
  }
  double solve() {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        if (!sign(x[i] - x[j]) && !sign(y[i] - y[j]) && !sign(r[i] - r[j])) {
          r[i] = 0.0;
          break;
        }
      }
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i != j && sign(r[j] - r[i]) >= 0 && sign(dist(x[i], y[i], x[j], y[j]) - (r[j] - r[i])) <= 0) {
          covered[i] = 1;
          break;
        }
      }
    }
    for (int i = 0; i < n; i++) {
      if (sign(r[i]) && !covered[i]) {
        seg.clear();
        for (int j = 0; j < n; j++) {
          if (i != j) {
            double d = dist(x[i], y[i], x[j], y[j]);
            if (sign(d - (r[j] + r[i])) >= 0 || sign(d - abs(r[j] - r[i])) <= 0) {
              continue;
            }
            double alpha = atan2(y[j] - y[i], x[j] - x[i]);
            double beta = angle(r[i], d, r[j]);
            pair<double, double> tmp(alpha - beta, alpha + beta);
            if (sign(tmp.first) <= 0 && sign(tmp.second) <= 0) {
              seg.push_back(pair<double, double>(2 * PI + tmp.first, 2 * PI + tmp.second));
            }
            else if (sign(tmp.first) < 0) {
              seg.push_back(pair<double, double>(2 * PI + tmp.first, 2 * PI));
              seg.push_back(pair<double, double>(0, tmp.second));
            }
            else {
              seg.push_back(tmp);
            }
          }
        }
        sort(seg.begin(), seg.end());
        double rig = 0;
        for (vector<pair<double, double> >::iterator iter = seg.begin(); iter != seg.end(); iter++) {
          if (sign(rig - iter->first) >= 0) {
            rig = max(rig, iter->second);
          }
          else {
            getarea(i, rig, iter->first);
            rig = iter->second;
          }
        }
        if (!sign(rig)) {
          arc += r[i] * r[i] * PI;
        }
        else {
          getarea(i, rig, 2 * PI);
        }
      }
    }
    return pol / 2.0 + arc;
  }
} CU; 
