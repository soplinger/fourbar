"""Tests for fourbar. Known answers come from closed-form geometry, not from the solver output."""

import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import fourbar  # noqa: E402
from fourbar import Linkage, Point, analyse, classify, gate, main, parse_text, solve  # noqa: E402

DEMO = str(ROOT / "examples" / "crank-rocker.lnk")
BROKEN = str(ROOT / "examples" / "broken.lnk")


def make(ground, crank, coupler, rocker, points=()):
    return Linkage(ground=ground, crank=crank, coupler=coupler, rocker=rocker,
                   points=list(points))


def demo_link():
    link, errors = fourbar.load(DEMO)
    assert errors == []
    return link


# ----------------------------------------------------------------- parsing

def test_demo_file_parses_into_expected_lengths():
    link = demo_link()
    assert (link.ground, link.crank, link.coupler, link.rocker) == (4.0, 1.0, 3.0, 2.5)
    assert link.name == "wiper-linkage"
    assert link.units == "mm"
    assert link.min_transmission == 40
    assert [p.name for p in link.points] == ["T", "P"]


def test_error_messages_carry_line_numbers():
    text = "ground 4\ncrank -1\ncoupler 3\nrocker 2\n"
    link, errors = parse_text(text, "x.lnk")
    assert link is None
    assert errors == ["x.lnk:2: 'crank' must be greater than zero"]


def test_missing_required_link_is_reported_without_a_line():
    link, errors = parse_text("ground 4\ncrank 1\n", "x.lnk")
    assert link is None
    assert "x.lnk: missing required 'coupler'" in errors
    assert "x.lnk: missing required 'rocker'" in errors


def test_broken_file_reports_every_problem():
    link, errors = fourbar.load(BROKEN)
    assert link is None
    joined = "\n".join(errors)
    assert ":4: 'crank' must be greater than zero" in joined
    assert ":6: unknown directive 'rockr'" in joined
    assert ":7: 'min_transmission' must be in (0, 90]" in joined
    assert ":8: point name 'A' is reserved or already used" in joined
    assert ":9: point angle is not a finite number: 'abc'" in joined
    assert ":10: duplicate 'ground'" in joined


# ----------------------------------------------------------- classification

@pytest.mark.parametrize(
    "lengths, expected",
    [
        ((4, 1, 3, 2.5), "crank-rocker"),
        ((1, 3, 4, 3.5), "double-crank"),
        ((4, 3, 1, 3), "double-rocker"),
        ((4, 3, 3, 1), "rocker-crank"),
        ((4, 1, 2, 2), "triple-rocker"),
    ],
)
def test_grashof_classes(lengths, expected):
    ground, crank, coupler, rocker = lengths
    assert classify(make(ground, crank, coupler, rocker))["class"] == expected


def test_change_point_is_flagged_when_sums_are_equal():
    info = classify(make(4, 1, 3, 2))
    assert info["change_point"] is True
    assert info["grashof"] is True
    assert classify(demo_link())["change_point"] is False


# ----------------------------------------------------------------- geometry

def test_transmission_angle_matches_law_of_cosines_at_crank_zero():
    # At theta = 0 the coupler-rocker triangle has AO4 = ground - crank = 3.
    link = demo_link()
    pose = solve(link, 0.0)
    expected = math.degrees(math.acos((3 ** 2 + 2.5 ** 2 - 3 ** 2) / (2 * 3 * 2.5)))
    # Written out for clarity: (9 + 6.25 - 9) / 15 = 6.25 / 15
    assert pose.mu == pytest.approx(expected, abs=1e-6)
    assert pose.mu == pytest.approx(math.degrees(math.acos(6.25 / 15)), abs=1e-6)


def test_joint_distances_match_link_lengths_at_many_angles():
    link = demo_link()
    o4 = (link.ground, 0.0)
    for theta in range(0, 360, 7):
        pose = solve(link, float(theta))
        assert pose is not None
        ax, ay = pose.A
        bx, by = pose.B
        assert math.hypot(ax, ay) == pytest.approx(link.crank, abs=1e-9)
        assert math.hypot(bx - ax, by - ay) == pytest.approx(link.coupler, abs=1e-9)
        assert math.hypot(bx - o4[0], by - o4[1]) == pytest.approx(link.rocker, abs=1e-9)


def test_crossed_branch_is_a_second_valid_assembly():
    link = demo_link()
    open_pose = solve(link, 30.0, "open")
    crossed_pose = solve(link, 30.0, "crossed")
    assert open_pose.B != pytest.approx(crossed_pose.B)
    assert open_pose.mu == pytest.approx(crossed_pose.mu, abs=1e-9)


def test_infeasible_crank_returns_none():
    # Ground 4, crank 1, coupler 2, rocker 2: at theta = 180 the crank reaches 5,
    # which is more than coupler + rocker can span from O4.
    link = make(4, 1, 2, 2)
    assert solve(link, 0.0) is not None
    assert solve(link, 180.0) is None


def test_coupler_point_at_zero_distance_sits_on_crank_tip():
    link = make(4, 1, 3, 2.5, [Point("Z", 0.0, 0.0)])
    pose = solve(link, 50.0)
    assert pose.points["Z"] == pytest.approx(pose.A, abs=1e-9)


def test_coupler_point_distance_from_crank_tip_is_its_offset():
    link = make(4, 1, 3, 2.5, [Point("P", 1.5, 37.0)])
    pose = solve(link, 110.0)
    px, py = pose.points["P"]
    assert math.hypot(px - pose.A[0], py - pose.A[1]) == pytest.approx(1.5, abs=1e-9)


# ----------------------------------------------------------------- analysis

def test_demo_has_full_crank_rotation_and_known_transmission_extremes():
    result = analyse(demo_link())
    assert result["feasible_samples"] == 360
    assert result["full_crank_rotation"] is True
    assert result["rocker"]["full_rotation"] is False
    assert 0 < result["rocker"]["swing_deg"] < 360
    lo = result["min_transmission"]
    hi = result["max_transmission"]
    # Transmission angle depends only on AO4, which runs from 3 (theta 0) to 5 (theta 180).
    assert lo["value"] == pytest.approx(math.degrees(math.acos(6.25 / 15)), abs=0.01)
    assert lo["theta"] == pytest.approx(0.0)
    assert hi["value"] == pytest.approx(math.degrees(math.acos((15.25 - 25) / 15)), abs=0.01)
    assert hi["theta"] == pytest.approx(180.0)


def test_non_grashof_mechanism_cannot_turn_the_crank_fully():
    result = analyse(make(4, 1, 2, 2))
    assert result["class"] == "triple-rocker"
    assert result["full_crank_rotation"] is False
    assert 0 < result["feasible_samples"] < 360


def test_near_toggle_position_is_found_at_straight_crank():
    # Coupler 3 plus rocker 2 equals the crank-to-ground reach of 5 at theta = 180.
    link = make(4, 1, 3, 2)
    result = analyse(link)
    thetas = [d["theta"] for d in result["dead_points"]]
    assert any(abs(t - 180.0) < 1e-6 for t in thetas)


def test_gate_passes_and_fails_around_the_known_minimum():
    result = analyse(demo_link())
    assert gate(result, 60.0)[0] == 0
    assert gate(result, 80.0)[0] == 2


def test_sweep_csv_has_header_and_one_row_per_sample():
    lines = fourbar.sweep_csv(demo_link(), 1.0, "open").splitlines()
    assert lines[0].startswith("theta_deg,feasible,rocker_deg,transmission_deg")
    assert lines[0].endswith("T_x,T_y,P_x,P_y")
    assert len(lines) == 361


def test_svg_contains_mechanism_and_trace():
    svg = fourbar.svg_drawing(demo_link(), 0.0, trace="T")
    assert svg.startswith("<svg")
    assert "<polyline" in svg
    assert svg.rstrip().endswith("</svg>")
    assert "<polyline" not in fourbar.svg_drawing(demo_link(), 0.0)


def test_ascii_trace_plots_points():
    text = fourbar.ascii_trace(demo_link(), "P", 1.0, "open")
    assert text.startswith("Trace of P")
    assert "*" in text


# ---------------------------------------------------------------------- CLI

def test_cli_check_accepts_demo(capsys):
    assert main(["check", DEMO]) == 0
    assert "crank-rocker" in capsys.readouterr().out


def test_cli_check_rejects_broken_file(capsys):
    assert main(["check", BROKEN]) == 1
    assert "error(s)" in capsys.readouterr().err


def test_cli_analyse_json_is_valid(capsys):
    assert main(["analyse", DEMO, "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["class"] == "crank-rocker"


def test_cli_at_prints_joints(capsys):
    assert main(["at", DEMO, "0"]) == 0
    out = capsys.readouterr().out
    assert "Crank 0.00 deg" in out
    assert "O4 (4.000, 0.000)" in out


def test_cli_gate_exit_codes():
    assert main(["gate", DEMO, "--min-transmission", "60"]) == 0
    assert main(["gate", DEMO, "--min-transmission", "80"]) == 2


def test_cli_rejects_bad_step_with_usage_error():
    with pytest.raises(SystemExit):
        main(["analyse", DEMO, "--step", "0"])


def test_cli_trace_rejects_unknown_point(capsys):
    assert main(["trace", DEMO, "Q"]) == 1
    assert "unknown point" in capsys.readouterr().err
