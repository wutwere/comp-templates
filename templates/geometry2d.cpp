template <class T>
struct Point {
  T x, y;

  Point(T x_ = 0, T y_ = 0) : x(x_), y(y_) {}

  Point operator+(Point p) { return Point<T>(x + p.x, y + p.y); }
  Point operator-(Point p) { return Point<T>(x - p.x, y - p.y); }
  Point operator*(T c) { return Point<T>(x * c, y * c); }
  Point operator/(T c) { return Point<T>(x / c, y / c); }

  bool operator==(Point p) { return x == p.x && y == p.y; }
  bool operator<(Point p) { return x < p.x || (x == p.x && y < p.y); } // for sorting

  inline T dot(Point p) { return (T) x * p.x + (T) y * p.y; }
  inline T cross(Point p) { return (T) x * p.y - (T) y * p.x; }
  inline T cross(Point a, Point b) { return (a - *this).cross(b - *this); }
  inline T dist2() { return x * x + y * y; }
};

typedef Point<int64_t> point;

// check if c is between line segment formed by a and b
bool point_between(point a, point b, point c) {
  if (a.cross(c, b)) return false;
  return min(a.x, b.x) <= c.x && c.x <= max(a.x, b.x) && min(a.y, b.y) <= c.y && c.y <= max(a.y, b.y);
}

// check if line segment formed by p[0] and p[1] intersects with line segment formed by p[2] and p[3]
bool intersect(vector<point> &&p) {
  int64_t a = p[0].cross(p[2], p[1]);
  int64_t b = p[0].cross(p[3], p[1]);
  int64_t c = p[2].cross(p[0], p[3]);
  int64_t d = p[2].cross(p[1], p[3]);
  if (a > b) swap(a, b);
  if (c > d) swap(c, d);
  if (!a && !b && !c && !d) {
    return point_between(p[0], p[1], p[2]) || point_between(p[0], p[1], p[3]) || point_between(p[2], p[3], p[0]) || point_between(p[2], p[3], p[1]);
  } else {
    return a <= 0 && 0 <= b && c <= 0 && 0 <= d;
  }
}

// check point inside polygon
int is_in_polygon(vector<point> &poly, point x) {
  const point far = point(1e9 + 5, 0);
  int check = 0;
  bool boundary = 0;
  for (int j = 0; j < n; j++) {
    if (point_between(poly[j], poly[(j + 1) % n], x)) {
      boundary = 1;
      break;
    }
    if (intersect({x, far, poly[j], poly[(j + 1) % n]})) check++;
  }
  if (boundary) return 0; // on boundary
  else if (check & 1) return 1; // inside
  else return -1; // outside
}

// find 2 * area of a polygon (points must be adjacent)
int64_t polygon_area(vector<point> &p) {
  int64_t ans = p.back().cross(p[0]);
  for (int i = 0; i < n; i++) {
    ans += p[i].cross(p[i + 1]);
  }
  return abs(ans);
}

// find convex hull of points (returned adjacent)
vector<point> convex_hull(vector<point> p) {
  sort(p.begin(), p.end());
  vector<point> top, bot;
  for (int i = 0; i < n; i++) {
    while (top.size() > 1 && top.back().cross(end(top)[-2], p[i]) < 0)
      top.pop_back();
    top.push_back(p[i]);
    while (bot.size() > 1 && bot.back().cross(end(bot)[-2], p[n - 1 - i]) < 0)
      bot.pop_back();
    bot.push_back(p[n - 1 - i]);
  }
  top.insert(top.end(), bot.begin() + 1, bot.end() - 1);
  return top;
}

typedef Point<double> P;
#define arg(p, q) atan2(p.cross(q), p.dot(q))
double circlePoly(P c, double r, vector<P> ps) {
  auto tri = [&](P p, P q) {
    auto r2 = r * r / 2;
    P d = q - p;
    auto a = d.dot(p)/d.dist2(), b = (p.dist2()-r*r)/d.dist2();
    auto det = a * a - b;
    if (det <= 0) return arg(p, q) * r2;
    auto s = max(0., -a-sqrt(det)), t = min(1., -a+sqrt(det));
    if (t < 0 || 1 <= s) return arg(p, q) * r2;
    P u = p + d * s, v = p + d * t;
    return arg(p,u) * r2 + u.cross(v)/2 + arg(v,q) * r2;
  };
  auto sum = 0.0;
  for (int i = 0; i < ps.size(); i++)
    sum += tri(ps[i] - c, ps[(i + 1) % ps.size()] - c);
  return sum;
}

