#!/usr/bin/env python3
"""Extract the fast-neutron data atomsim needs from an ENDF-6 evaluation.

Usage: extract_endf.py <evaluation.endf> <output.xs>

Writes, for incident energies 25 keV - 20 MeV:
  XS  MF3 cross-sections (elastic, inelastic levels and continuum, (n,2n),
      fission, capture) with their Q-values
  NU  MF1 MT452 total nu-bar(E)
  ANG MF4 MT2 elastic angular distributions in the centre-of-mass frame,
      as 32 equiprobable bins of mu per incident energy
Source: ENDF/B-VIII.0 (Brown et al., Nucl. Data Sheets 148, 1, 2018).
"""
import math
import sys

EMIN, EMAX = 2.5e4, 2.0e7  # eV
BINS = 32


def num(field):
    """ENDF floats omit the 'e': '1.234567+5' or '-2.5-3'."""
    s = field.strip()
    if not s:
        return 0.0
    for i in range(len(s) - 1, 0, -1):
        if s[i] in "+-" and s[i - 1] not in "eE":
            return float(s[:i] + "e" + s[i:])
    return float(s)


def sections(path):
    out = {}
    with open(path) as f:
        for line in f:
            mf, mt = int(line[70:72] or 0), int(line[72:75] or 0)
            if mf and mt:
                out.setdefault((mf, mt), []).append(line)
    return out


class Reader:
    def __init__(self, lines):
        self.lines, self.i = lines, 0

    def cont(self):
        line = self.lines[self.i]
        self.i += 1
        return [num(line[k:k + 11]) for k in range(0, 66, 11)]

    def values(self, n):
        vals = []
        while len(vals) < n:
            vals += self.cont()
        return vals[:n]

    def tab1(self):
        c = self.cont()
        nr, np_ = int(c[4]), int(c[5])
        self.values(2 * nr)  # interpolation laws: lin-lin is used throughout
        xy = self.values(2 * np_)
        return c, list(zip(xy[0::2], xy[1::2]))

    def tab2(self):
        c = self.cont()
        self.values(2 * int(c[4]))
        return c

    def list_(self):
        c = self.cont()
        return c, self.values(int(c[4]))


def legendre_table(lmax, xs):
    """P_l(x) for l = 0..lmax at each x, by the three-term recursion."""
    table = [[1.0] * len(xs), list(xs)]
    for l in range(2, lmax + 1):
        table.append([((2 * l - 1) * x * p1 - (l - 1) * p0) / l
                      for x, p0, p1 in zip(xs, table[l - 2], table[l - 1])])
    return table


def equiprobable(mu, pdf):
    """Edges of BINS equal-probability bins of a tabulated pdf(mu)."""
    pdf = [max(p, 0.0) for p in pdf]
    cdf = [0.0]
    for k in range(1, len(mu)):
        cdf.append(cdf[-1] + 0.5 * (pdf[k] + pdf[k - 1]) * (mu[k] - mu[k - 1]))
    total = cdf[-1]
    edges, k = [-1.0], 0
    for b in range(1, BINS):
        target = total * b / BINS
        while cdf[k + 1] < target:
            k += 1
        frac = (target - cdf[k]) / (cdf[k + 1] - cdf[k]) if cdf[k + 1] > cdf[k] else 0.0
        edges.append(mu[k] + frac * (mu[k + 1] - mu[k]))
    return edges + [1.0]


def angular(lines):
    r = Reader(lines)
    head = r.cont()
    ltt = int(head[3])
    r.cont()  # LI, LCT (2 = centre of mass)
    grid = [-1 + 2 * k / 400 for k in range(401)]
    poly = legendre_table(100, grid)
    table = []
    if ltt in (1, 3):
        ne = int(r.tab2()[5])
        for _ in range(ne):
            c, coeffs = r.list_()
            pdf = [0.5] * len(grid)
            for l, a in enumerate(coeffs, start=1):
                w = (2 * l + 1) / 2 * a
                pdf = [p + w * pl for p, pl in zip(pdf, poly[l])]
            table.append((c[1], equiprobable(grid, pdf)))
    if ltt in (2, 3):
        ne = int(r.tab2()[5])
        for _ in range(ne):
            c, pts = r.tab1()
            table.append((c[1], equiprobable([p[0] for p in pts], [p[1] for p in pts])))
    return table


def clip(pts):
    """Points inside [EMIN, EMAX], keeping one neighbour on each side."""
    keep = [k for k, (e, _) in enumerate(pts) if EMIN <= e <= EMAX]
    if not keep:
        return []
    lo, hi = max(keep[0] - 1, 0), min(keep[-1] + 1, len(pts) - 1)
    return pts[lo:hi + 1]


def main(src, dst):
    sec = sections(src)
    awr = Reader(sec[(3, 1)]).cont()[1]
    out = ["# ENDF/B-VIII.0 fast-neutron data extracted by tools/extract_endf.py from "
           + src.split("/")[-1], "AWR %.6f" % awr]
    wanted = [2, 16, 18, 102] + list(range(51, 91)) + [91]
    for mt in wanted:
        if (3, mt) not in sec:
            continue
        r = Reader(sec[(3, mt)])
        r.cont()
        c, pts = r.tab1()
        pts = clip(pts)
        if not pts or max(y for _, y in pts) <= 0:
            continue
        out.append("XS %d %.6e %d" % (mt, c[1], len(pts)))
        out += ["%.6e %.6e" % p for p in pts]

    r = Reader(sec[(1, 452)])
    head = r.cont()
    if int(head[3]) == 2:
        _, pts = r.tab1()
    else:  # polynomial in E
        _, coeffs = r.list_()
        pts = [(e, sum(a * e ** k for k, a in enumerate(coeffs)))
               for e in (EMIN * (EMAX / EMIN) ** (k / 50) for k in range(51))]
    pts = clip(pts)
    out.append("NU %d" % len(pts))
    out += ["%.6e %.6f" % p for p in pts]

    ang = [row for row in angular(sec[(4, 2)]) if EMIN / 2 <= row[0] <= EMAX * 1.5]
    out.append("ANG %d %d" % (len(ang), BINS + 1))
    out += ["%.6e " % e + " ".join("%.5f" % x for x in edges) for e, edges in ang]
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
