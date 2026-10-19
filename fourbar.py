#!/usr/bin/env python3
"""Fourbar: offline four-bar linkage analyser driven by a plain-text .lnk file.

Geometry: the fixed pivots are O2 = (0, 0) and O4 = (ground, 0). The crank
O2A turns by theta, and the coupler AB and rocker O4B close the loop. The
rocker angle phi is measured from +x at O4. Coupler points are given as a
distance along AB and an angle (degrees) from the AB direction, positive
counter-clockwise.

Exit codes: 0 success, 1 file or lint error, 2 gate failed or no assembly.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

LINKS = ("ground", "crank", "coupler", "rocker")
DEFAULT_STEP = 1.0
DEFAULT_MIN_TRANSMISSION = 40.0
DEAD_TOLERANCE_DEG = 1.0
RESERVED_POINTS = {"A", "B", "O2", "O4"}
NAME_RE = re.compile(r"[A-Za-z][A-Za-z0-9_]*")
CLASS_TEXT = {
    "crank-rocker": "crank rotates fully, rocker oscillates",
    "rocker-crank": "rocker rotates fully, crank oscillates",
    "double-crank": "both ground-pivoted links rotate fully",
    "double-rocker": "neither ground-pivoted link rotates fully",
    "triple-rocker": "non-Grashof: no link rotates fully",
}


@dataclass(frozen=True)
class Point:
    name: str
    dist: float
    angle: float


@dataclass
class Linkage:
    name: str = "unnamed"
    ground: float = 0.0
    crank: float = 0.0
    coupler: float = 0.0
    rocker: float = 0.0
    min_transmission: float = DEFAULT_MIN_TRANSMISSION
    units: str = ""
    points: list[Point] = field(default_factory=list)


@dataclass(frozen=True)
class Pose:
    A: tuple[float, float]
    B: tuple[float, float]
    phi: float
    mu: float
    points: dict[str, tuple[float, float]]


# ---------------------------------------------------------------- parsing

def _number(token: str) -> float | None:
    try:
        value = float(token)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def parse_text(text: str, source: str = "<input>") -> tuple[Linkage | None, list[str]]:
    """Parse .lnk text. Returns (linkage or None, list of 'source:line: message')."""
    link = Linkage()
    errors: list[tuple[int, str]] = []
    seen: set[str] = set()
    point_names: set[str] = set()

    for lineno, raw in enumerate(text.splitlines(), start=1):
        parts = raw.split("#", 1)[0].split()
        if not parts:
            continue
        key, args = parts[0].lower(), parts[1:]

        if key == "point":
            if len(args) != 3:
                errors.append((lineno, "point needs NAME DISTANCE ANGLE"))
                continue
            name, dist_text, angle_text = args
            dist, angle = _number(dist_text), _number(angle_text)
            if not NAME_RE.fullmatch(name):
                errors.append((lineno, f"invalid point name {name!r}"))
            elif name in RESERVED_POINTS or name in point_names:
                errors.append((lineno, f"point name {name!r} is reserved or already used"))
            elif dist is None or dist < 0:
                errors.append((lineno, "point distance must be a number >= 0"))
            elif angle is None:
                errors.append((lineno, f"point angle is not a finite number: {angle_text!r}"))
            else:
                link.points.append(Point(name, dist, angle))
                point_names.add(name)
            continue

        if key not in ("name", "units", "min_transmission") + LINKS:
            errors.append((lineno, f"unknown directive {key!r}"))
            continue
        if key in seen:
            errors.append((lineno, f"duplicate '{key}'"))
            continue
        seen.add(key)

        if key == "name":
            if args:
                link.name = " ".join(args)
            else:
                errors.append((lineno, "'name' needs a value"))
        elif key == "units":
            if len(args) == 1:
                link.units = args[0]
            else:
                errors.append((lineno, "'units' takes one word"))
        else:
            if len(args) != 1:
                errors.append((lineno, f"'{key}' takes one number"))
                continue
            value = _number(args[0])
            if value is None:
                errors.append((lineno, f"'{key}' is not a finite number: {args[0]!r}"))
            elif key == "min_transmission":
                if 0 < value <= 90:
                    link.min_transmission = value
                else:
                    errors.append((lineno, "'min_transmission' must be in (0, 90]"))
            elif value <= 0:
                errors.append((lineno, f"'{key}' must be greater than zero"))
            else:
                setattr(link, key, value)

    for key in LINKS:
        if key not in seen:
            errors.append((0, f"missing required '{key}'"))

    messages = [
        f"{source}:{line}: {msg}" if line else f"{source}: {msg}"
        for line, msg in sorted(errors, key=lambda item: (item[0] == 0, item[0]))
    ]
    return (link if not messages else None), messages


def load(path: str) -> tuple[Linkage | None, list[str]]:
    try:
        text = Path(path).read_text(encoding="utf-8")
    except (OSError, ValueError) as exc:
        return None, [f"{path}: cannot read file ({exc})"]
    return parse_text(text, source=path)


# ---------------------------------------------------------------- geometry

def solve(link: Linkage, theta_deg: float, branch: str = "open") -> Pose | None:
    """Return the assembly at a crank angle, or None if the links cannot close."""
    a, b, c, d = link.crank, link.coupler, link.rocker, link.ground
    t = math.radians(theta_deg)
    ax, ay = a * math.cos(t), a * math.sin(t)
    px, py = ax - d, ay
    reach = math.hypot(px, py)
    if reach == 0:
        return None
    # Rocker circle about O4 meets the coupler circle about A: u.P = k.
    k = (c * c + reach * reach - b * b) / (2 * c)
    ratio = k / reach
    if abs(ratio) > 1 + 1e-9:
        return None
    ratio = max(-1.0, min(1.0, ratio))
    psi = math.atan2(py, px)
    offset = math.acos(ratio)
    phi = psi - offset if branch == "crossed" else psi + offset
    bx, by = d + c * math.cos(phi), c * math.sin(phi)

    # Transmission angle at B, between BA and BO4.
    cos_mu = ((ax - bx) * (d - bx) + (ay - by) * (-by)) / (b * c)
    mu = math.degrees(math.acos(max(-1.0, min(1.0, cos_mu))))

    ex, ey = (bx - ax) / b, (by - ay) / b
    nx, ny = -ey, ex
    points = {}
    for p in link.points:
        alpha = math.radians(p.angle)
        ca, sa = math.cos(alpha), math.sin(alpha)
        points[p.name] = (
            ax + p.dist * (ca * ex + sa * nx),
            ay + p.dist * (ca * ey + sa * ny),
        )
    return Pose(A=(ax, ay), B=(bx, by), phi=math.degrees(phi), mu=mu, points=points)


def sweep(link: Linkage, step: float = DEFAULT_STEP, branch: str = "open"):
    """Sample the crank over one turn. Returns [(theta_deg, Pose or None)]."""
    n = max(4, round(360.0 / step))
    out = []
    for i in range(n):
        theta = i * 360.0 / n
        out.append((theta, solve(link, theta, branch)))
    return out


def classify(link: Linkage) -> dict:
    values = {name: getattr(link, name) for name in LINKS}
    s = min(values.values())
    l = max(values.values())
    total = sum(values.values())
    pq = total - s - l
    grashof = s + l <= pq + 1e-9
    shortest = next(name for name in LINKS if values[name] == s)
    if not grashof:
        kind = "triple-rocker"
    elif shortest == "ground":
        kind = "double-crank"
    elif shortest == "crank":
        kind = "crank-rocker"
    elif shortest == "rocker":
        kind = "rocker-crank"
    else:
        kind = "double-rocker"
    return {
        "class": kind,
        "grashof": grashof,
        "shortest": shortest,
        "change_point": abs(s + l - pq) <= 1e-9,
        "sum_sl": round(s + l, 6),
        "sum_pq": round(pq, 6),
    }


def feasible_runs(ok: list[bool]) -> list[list[int]]:
    """Contiguous runs of feasible sample indices, treating the crank as circular."""
    n = len(ok)
    if not any(ok):
        return []
    if all(ok):
        return [list(range(n))]
    start = ok.index(False)
    runs: list[list[int]] = []
    current: list[int] = []
    for j in range(1, n + 1):
        i = (start + j) % n
        if ok[i]:
            current.append(i)
        elif current:
            runs.append(current)
            current = []
    return runs


def unwrap(angles: list[float]) -> list[float]:
    out = [angles[0]]
    for a in angles[1:]:
        delta = (a - out[-1] + 180.0) % 360.0 - 180.0
        out.append(out[-1] + delta)
    return out


def _toggle_gap(mu: float) -> float:
    return min(mu, 180.0 - mu)


def _rocker_motion(samples, runs, n) -> dict:
    best = {"full_rotation": False, "swing_deg": None, "theta_range": None}
    for run in runs:
        phis = [samples[i][1].phi for i in run]
        if len(run) == n:
            loop = unwrap(phis + phis[:1])
            if abs(loop[-1] - loop[0]) > 180:
                return {"full_rotation": True, "swing_deg": None, "theta_range": None}
            span = loop[:-1]
        else:
            span = unwrap(phis)
        swing = max(span) - min(span)
        if best["swing_deg"] is None or swing > best["swing_deg"]:
            best = {
                "full_rotation": False,
                "swing_deg": round(swing, 4),
                "theta_range": [round(samples[run[0]][0], 4), round(samples[run[-1]][0], 4)],
            }
    return best


def _dead_points(samples, runs, n) -> list[dict]:
    dead = []
    for run in runs:
        length = len(run)
        full = length == n
        for p in range(length):
            if not full and p in (0, length - 1):
                continue
            i = run[p]
            gap = _toggle_gap(samples[i][1].mu)
            neighbours = (run[(p - 1) % length], run[(p + 1) % length])
            if gap <= DEAD_TOLERANCE_DEG and all(
                gap <= _toggle_gap(samples[j][1].mu) for j in neighbours
            ):
                dead.append({"theta": round(samples[i][0], 4), "gap_deg": round(gap, 4)})
    return dead


def analyse(link: Linkage, step: float = DEFAULT_STEP, branch: str = "open") -> dict:
    samples = sweep(link, step, branch)
    n = len(samples)
    ok = [pose is not None for _, pose in samples]
    cls = classify(link)
    result = {
        "name": link.name,
        "units": link.units,
        "branch": branch,
        "step": step,
        "lengths": {k: getattr(link, k) for k in LINKS},
        "grashof": cls["grashof"],
        "shortest": cls["shortest"],
        "class": cls["class"],
        "change_point": cls["change_point"],
        "sum_sl": cls["sum_sl"],
        "sum_pq": cls["sum_pq"],
        "samples": n,
        "feasible_samples": sum(ok),
        "full_crank_rotation": all(ok),
        "crank_runs": [],
        "rocker": {"full_rotation": False, "swing_deg": None, "theta_range": None},
        "min_transmission": None,
        "max_transmission": None,
        "limit": link.min_transmission,
        "dead_points": [],
        "points": {},
    }
    if not any(ok):
        return result

    runs = feasible_runs(ok)
    result["crank_runs"] = [
        [round(samples[run[0]][0], 4), round(samples[run[-1]][0], 4)] for run in runs
    ]
    result["rocker"] = _rocker_motion(samples, runs, n)
    result["dead_points"] = _dead_points(samples, runs, n)

    feasible = [(theta, pose) for theta, pose in samples if pose is not None]
    lo_theta, lo = min(feasible, key=lambda item: item[1].mu)
    hi_theta, hi = max(feasible, key=lambda item: item[1].mu)
    result["min_transmission"] = {
        "value": round(lo.mu, 4),
        "theta": round(lo_theta, 4),
        "ok": lo.mu >= link.min_transmission,
    }
    result["max_transmission"] = {"value": round(hi.mu, 4), "theta": round(hi_theta, 4)}

    for p in link.points:
        xs = [pose.points[p.name][0] for _, pose in feasible]
        ys = [pose.points[p.name][1] for _, pose in feasible]
        result["points"][p.name] = {
            "x": [round(min(xs), 4), round(max(xs), 4)],
            "y": [round(min(ys), 4), round(max(ys), 4)],
        }
    return result


def gate(result: dict, limit: float) -> tuple[int, str]:
    if result["feasible_samples"] == 0:
        return 2, "FAIL: no crank position closes the linkage"
    lowest = result["min_transmission"]["value"]
    if lowest < limit:
        return 2, f"FAIL: minimum transmission angle {lowest:.2f} deg is below {limit:.2f} deg"
    return 0, f"PASS: minimum transmission angle {lowest:.2f} deg meets {limit:.2f} deg"


# ---------------------------------------------------------------- output

def format_report(r: dict) -> str:
    unit = f" ({r['units']})" if r["units"] else ""
    lengths = r["lengths"]
    lines = [
        f"Four-bar: {r['name']}{unit}",
        f"  ground {lengths['ground']:.2f}  crank {lengths['crank']:.2f}"
        f"  coupler {lengths['coupler']:.2f}  rocker {lengths['rocker']:.2f}",
        f"Grashof: {'yes' if r['grashof'] else 'no'}"
        f"  (shortest + longest = {r['sum_sl']:.2f}, other two = {r['sum_pq']:.2f})",
    ]
    change = "; change-point (can lock at a toggle)" if r["change_point"] else ""
    lines.append(f"Class: {r['class']} - {CLASS_TEXT[r['class']]}{change}")
    lines.append(
        f"Branch: {r['branch']}   sweep: {r['samples']} samples, "
        f"{r['feasible_samples']} feasible"
    )
    if r["feasible_samples"] == 0:
        lines.append("No crank position closes the linkage.")
        return "\n".join(lines)

    if r["full_crank_rotation"]:
        lines.append("Crank range: full rotation (0 to 360 deg)")
    else:
        spans = "; ".join(f"{a:.1f} to {b:.1f} deg" for a, b in r["crank_runs"])
        lines.append(f"Crank range: {spans}")

    rocker = r["rocker"]
    if rocker["full_rotation"]:
        lines.append("Rocker: rotates fully through 360 deg")
    elif rocker["swing_deg"] is not None:
        a, b = rocker["theta_range"]
        lines.append(
            f"Rocker swing: {rocker['swing_deg']:.2f} deg (crank {a:.1f} to {b:.1f} deg)"
        )

    lo, hi = r["min_transmission"], r["max_transmission"]
    verdict = "OK" if lo["ok"] else "FAIL"
    lines.append(
        f"Transmission angle: min {lo['value']:.2f} deg at crank {lo['theta']:.1f}"
        f" (limit {r['limit']:.1f}, {verdict}); max {hi['value']:.2f} deg"
        f" at crank {hi['theta']:.1f}"
    )
    if r["dead_points"]:
        items = ", ".join(
            f"crank {d['theta']:.1f} deg (gap {d['gap_deg']:.2f})" for d in r["dead_points"]
        )
        lines.append(f"Near-toggle (dead) positions: {items}")
    else:
        lines.append("Near-toggle (dead) positions: none")
    for name, box in r["points"].items():
        lines.append(
            f"Point {name}: x {box['x'][0]:.2f}..{box['x'][1]:.2f}, "
            f"y {box['y'][0]:.2f}..{box['y'][1]:.2f}"
        )
    return "\n".join(lines)


def format_pose(link: Linkage, theta: float, branch: str) -> str | None:
    pose = solve(link, theta, branch)
    if pose is None:
        return None
    ax, ay = pose.A
    bx, by = pose.B
    lines = [
        f"Crank {theta:.2f} deg ({branch} branch)",
        f"  O2 (0.000, 0.000)  A ({ax:.3f}, {ay:.3f})  B ({bx:.3f}, {by:.3f})"
        f"  O4 ({link.ground:.3f}, 0.000)",
        f"  rocker angle {pose.phi:.2f} deg   transmission angle {pose.mu:.2f} deg",
    ]
    for name, (x, y) in pose.points.items():
        lines.append(f"  {name} ({x:.3f}, {y:.3f})")
    return "\n".join(lines)


def sweep_csv(link: Linkage, step: float, branch: str) -> str:
    buf = io.StringIO()
    writer = csv.writer(buf)
    names = [p.name for p in link.points]
    writer.writerow(
        ["theta_deg", "feasible", "rocker_deg", "transmission_deg", "B_x", "B_y"]
        + [f"{name}_{axis}" for name in names for axis in "xy"]
    )
    for theta, pose in sweep(link, step, branch):
        if pose is None:
            writer.writerow([f"{theta:.4f}", 0] + [""] * (4 + 2 * len(names)))
            continue
        row = [
            f"{theta:.4f}", 1, f"{pose.phi:.6f}", f"{pose.mu:.6f}",
            f"{pose.B[0]:.6f}", f"{pose.B[1]:.6f}",
        ]
        for name in names:
            x, y = pose.points[name]
            row += [f"{x:.6f}", f"{y:.6f}"]
        writer.writerow(row)
    return buf.getvalue()


def ascii_trace(link: Linkage, point: str, step: float, branch: str,
                width: int = 60, height: int = 18) -> str | None:
    pts = [pose.points[point] for _, pose in sweep(link, step, branch) if pose is not None]
    if not pts:
        return None
    x0, x1 = min(x for x, _ in pts), max(x for x, _ in pts)
    y0, y1 = min(y for _, y in pts), max(y for _, y in pts)
    if x1 - x0 < 1e-9:
        x0, x1 = x0 - 1, x1 + 1
    if y1 - y0 < 1e-9:
        y0, y1 = y0 - 1, y1 + 1
    grid = [[" "] * width for _ in range(height)]
    for x, y in pts:
        col = round((x - x0) / (x1 - x0) * (width - 1))
        row = round((y1 - y) / (y1 - y0) * (height - 1))
        grid[row][col] = "*"
    border = "+" + "-" * width + "+"
    body = ["|" + "".join(row) + "|" for row in grid]
    header = (
        f"Trace of {point} over feasible crank angles "
        f"(x {x0:.2f}..{x1:.2f}, y {y0:.2f}..{y1:.2f})"
    )
    return "\n".join([header, border, *body, border])


def svg_drawing(link: Linkage, theta: float, branch: str = "open",
                step: float = DEFAULT_STEP, trace: str | None = None,
                width_px: int = 600) -> str:
    pose = solve(link, theta, branch)
    if pose is None:
        raise ValueError(f"no assembly at crank {theta:g} deg ({branch} branch)")
    feasible = [p for _, p in sweep(link, step, branch) if p is not None]

    xs, ys = [0.0, link.ground], [0.0, 0.0]
    for p in [*feasible, pose]:
        for x, y in (p.A, p.B, *p.points.values()):
            xs.append(x)
            ys.append(y)
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    margin = 0.08 * max(x1 - x0, y1 - y0, 1.0)
    x0, x1, y0, y1 = x0 - margin, x1 + margin, y0 - margin, y1 + margin
    w, h = x1 - x0, y1 - y0
    height_px = round(width_px * h / w)
    stroke = 0.012 * max(w, h)
    radius = 0.025 * max(w, h)
    font = 0.04 * max(w, h)

    def seg(a, b, color, dash=""):
        extra = f' stroke-dasharray="{dash}"' if dash else ""
        return (f'<line x1="{a[0]:.4f}" y1="{a[1]:.4f}" x2="{b[0]:.4f}" y2="{b[1]:.4f}"'
                f' stroke="{color}" stroke-width="{stroke:.4f}"{extra}/>')

    def dot(p, color):
        return f'<circle cx="{p[0]:.4f}" cy="{p[1]:.4f}" r="{radius:.4f}" fill="{color}"/>'

    o2, o4 = (0.0, 0.0), (link.ground, 0.0)
    body = [seg(o2, o4, "#999999", "4 3"), seg(o2, pose.A, "#1f6feb"),
            seg(pose.A, pose.B, "#222222"), seg(o4, pose.B, "#cf222e")]
    if trace is not None:
        trail = " ".join(f"{p.points[trace][0]:.4f},{p.points[trace][1]:.4f}" for p in feasible)
        body.append(f'<polyline points="{trail}" fill="none" stroke="#1a7f37"'
                    f' stroke-width="{stroke / 2:.4f}" stroke-dasharray="{stroke * 2:.4f} {stroke:.4f}"/>')
    body += [dot(o2, "#333333"), dot(o4, "#333333"), dot(pose.A, "#1f6feb"), dot(pose.B, "#cf222e")]
    body += [dot(p, "#1a7f37") for p in pose.points.values()]

    labels = [("O2", o2), ("O4", o4), ("A", pose.A), ("B", pose.B)]
    labels += list(pose.points.items())
    text = [f'<text x="{x + radius:.4f}" y="{-y - radius:.4f}" font-size="{font:.4f}"'
            f' font-family="sans-serif">{name}</text>' for name, (x, y) in labels]

    return "\n".join([
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width_px}" height="{height_px}"'
        f' viewBox="{x0:.4f} {-y1:.4f} {w:.4f} {h:.4f}">',
        f'<rect x="{x0:.4f}" y="{-y1:.4f}" width="{w:.4f}" height="{h:.4f}" fill="white"/>',
        '<g transform="scale(1,-1)">',
        *body,
        "</g>",
        *text,
        "</svg>",
        "",
    ])


# ---------------------------------------------------------------- CLI

def _step(text: str) -> float:
    value = float(text)
    if not 0.05 <= value <= 90:
        raise argparse.ArgumentTypeError("step must be between 0.05 and 90 degrees")
    return value


def cmd_check(link: Linkage, args) -> int:
    cls = classify(link)
    print(f"OK {args.file}: {link.name} is a {cls['class']} ({CLASS_TEXT[cls['class']]})")
    return 0


def cmd_analyse(link: Linkage, args) -> int:
    result = analyse(link, args.step, args.branch)
    print(json.dumps(result, indent=2) if args.json else format_report(result))
    return 0


def cmd_at(link: Linkage, args) -> int:
    text = format_pose(link, args.angle, args.branch)
    if text is None:
        print(f"no assembly at crank {args.angle:g} deg ({args.branch} branch): "
              "links cannot close", file=sys.stderr)
        return 2
    print(text)
    return 0


def cmd_sweep(link: Linkage, args) -> int:
    print(sweep_csv(link, args.step, args.branch), end="")
    return 0


def cmd_trace(link: Linkage, args) -> int:
    if args.point not in {p.name for p in link.points}:
        names = ", ".join(p.name for p in link.points) or "none defined"
        print(f"unknown point {args.point!r} (available: {names})", file=sys.stderr)
        return 1
    text = ascii_trace(link, args.point, args.step, args.branch, args.width, args.height)
    if text is None:
        print("no feasible crank positions to trace", file=sys.stderr)
        return 2
    print(text)
    return 0


def cmd_svg(link: Linkage, args) -> int:
    if args.trace is not None and args.trace not in {p.name for p in link.points}:
        print(f"unknown point {args.trace!r}", file=sys.stderr)
        return 1
    try:
        svg = svg_drawing(link, args.angle, args.branch, args.step, args.trace)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if args.output:
        Path(args.output).write_text(svg, encoding="utf-8")
        print(f"wrote {args.output}")
    else:
        print(svg, end="")
    return 0


def cmd_gate(link: Linkage, args) -> int:
    limit = args.min_transmission if args.min_transmission is not None else link.min_transmission
    result = analyse(link, args.step, args.branch)
    code, message = gate(result, limit)
    print(message)
    return code


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fourbar",
        description="Offline four-bar linkage analyser driven by a plain-text .lnk file.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def add(name, func, help_text, sweep_opts=True):
        p = sub.add_parser(name, help=help_text, description=help_text)
        p.add_argument("file", help="linkage .lnk file")
        if sweep_opts:
            p.add_argument("--step", type=_step, default=DEFAULT_STEP,
                           help="crank sweep step in degrees (default 1)")
            p.add_argument("--branch", choices=("open", "crossed"), default="open",
                           help="assembly branch (default open)")
        p.set_defaults(func=func)
        return p

    add("check", cmd_check, "lint the file and report the mechanism class", sweep_opts=False)
    p = add("analyse", cmd_analyse, "full analysis: class, ranges, transmission, dead points")
    p.add_argument("--json", action="store_true", help="print the analysis as JSON")
    p = add("at", cmd_at, "joint positions at one crank angle")
    p.add_argument("angle", type=float, help="crank angle in degrees")
    add("sweep", cmd_sweep, "CSV of every sampled crank angle")
    p = add("trace", cmd_trace, "ASCII plot of a coupler point's path")
    p.add_argument("point", help="coupler point name")
    p.add_argument("--width", type=int, default=60)
    p.add_argument("--height", type=int, default=18)
    p = add("svg", cmd_svg, "SVG drawing of the mechanism at one crank angle")
    p.add_argument("angle", type=float, help="crank angle in degrees")
    p.add_argument("--trace", help="draw the path of this coupler point")
    p.add_argument("-o", "--output", help="write SVG to this file instead of stdout")
    p = add("gate", cmd_gate, "CI gate: fail if the minimum transmission angle is too low")
    p.add_argument("--min-transmission", type=float, default=None,
                   help="override the file's min_transmission limit (degrees)")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    link, errors = load(args.file)
    if link is None:
        for message in errors:
            print(message, file=sys.stderr)
        print(f"{len(errors)} error(s) in {args.file}", file=sys.stderr)
        return 1
    return args.func(link, args)


if __name__ == "__main__":
    sys.exit(main())
