double simpson(double (*f)(double), double a, double b) {
  double c = (a + b) / 2;
  return (f(a) + 4 * f(c) + f(b)) * (b - a) / 6;
}

double rec(double (*f)(double), double a, double b, double eps, double S) {
  double c = (a + b) / 2;
  double S1 = simpson(f, a, c);
  double S2 = simpson(f, c, b), T = S1 + S2;
  if (abs(T - S) <= 15 * eps || b - a < 1e-10)
    return T + (T - S) / 15;
  return rec(f, a, c, eps / 2, S1) + rec(f, c, b, eps / 2, S2);
}

double integrate(double (*f)(double), double a, double b, double eps = 1e-8) {
  return rec(f, a, b, eps, simpson(f, a, b));
}

double y;
double f(double x) {
  return x * x + y * y <= 1;
}
double g(double yy) {
  y = yy;
  return integrate(f, -1, 1);
}

integrate(g, -1, 1); // returns pi
