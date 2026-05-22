"""
make_interpolators.py
=====================
Command-line tool for generating and testing interpolator JSON files.

Usage examples
--------------
# Generate the default uniform-grid interpolators.json:
  python make_interpolators.py generate

# Generate the adaptive PCHIP interpolator_adaptive.json:
  python make_interpolators.py generate-adaptive

# Test an existing interpolators.json against direct integration:
  python make_interpolators.py test

# Test an existing interpolator_adaptive.json:
  python make_interpolators.py test-adaptive

# All options have --help:
  python make_interpolators.py generate --help
"""

import argparse
import sys


def cmd_generate(args):
    from interpolators import generate_interpolators_file

    generate_interpolators_file(
        output_file=args.output,
        n0_values=range(args.n0_min, args.n0_max + 1),
        l_max=args.l_max,
        qr_min=args.qr_min,
        qr_max=args.qr_max,
        num_qr=args.num_qr,
        num_workers=args.workers,
    )


def cmd_generate_adaptive(args):
    from interpolators import generate_adaptive_interpolators_file

    generate_adaptive_interpolators_file(
        output_file=args.output,
        n0_values=range(args.n0_min, args.n0_max + 1),
        l_max=args.l_max,
        qr_min=args.qr_min,
        qr_max=args.qr_max,
        rel_tol=args.rel_tol,
        abs_tol=args.abs_tol,
        initial_points=args.initial_points,
        max_points=args.max_points,
        num_workers=args.workers,
        auto_qr_max=(not args.disable_auto_qr_max),
        auto_threshold=args.auto_threshold,
        low_q_switch=args.low_q_switch,
        low_q_anchor_points=args.low_q_anchor_points,
        auto_low_q_switch=(not args.disable_auto_low_q_switch),
        low_q_change_threshold=args.low_q_change_threshold,
        phase_tol=args.phase_tol,
        use_q_scaling=(not args.disable_q_scaling),
        phase_amp_rel_floor=args.phase_amp_rel_floor,
        min_interval_rel=args.min_interval_rel,
        min_interval_abs=args.min_interval_abs,
        max_new_points_per_iter=args.max_new_points_per_iter,
        cutoff_file=args.cutoff_file,
        auto_use_cutoffs=(not args.disable_auto_cutoffs),
        cutoff_field=args.cutoff_field,
        cutoff_threshold=args.cutoff_threshold,
        cutoff_pad_factor=args.cutoff_pad_factor,
        cutoff_min_qr_max=args.cutoff_min_qr_max,
        zero_component_threshold=args.zero_component_threshold,
        use_tail_fit=args.use_tail_fit,
        tail_fit_start_q=args.tail_fit_start_q,
        fit_num_cycles=args.fit_num_cycles,
        fit_scan_qr_max=args.fit_scan_qr_max,
        fit_scan_points=args.fit_scan_points,
        use_point_cache=(not args.disable_point_cache),
        point_cache_dir=args.point_cache_dir,
    )


def cmd_test(args):
    from interpolators import test_precomputed_interpolators

    result = test_precomputed_interpolators(
        threshold_percent=args.threshold,
        input_file=args.input,
        n0_values=range(args.n0_min, args.n0_max + 1),
        l_max=args.l_max,
        max_interpolators=args.max_interpolators,
    )
    if result["flagged_count"] > 0:
        sys.exit(1)


def cmd_test_adaptive(args):
    from interpolators import test_precomputed_adaptive_interpolators

    result = test_precomputed_adaptive_interpolators(
        threshold_percent=args.threshold,
        input_file=args.input,
        n0_values=range(args.n0_min, args.n0_max + 1),
        l_max=args.l_max,
        max_interpolators=args.max_interpolators,
        plot_output_dir=args.plot_dir,
        plot_points=args.plot_points,
        plot_qr_max=args.plot_qr_max,
        num_workers=args.workers,
        use_point_cache=not args.disable_point_cache,
        point_cache_dir=args.point_cache_dir,
    )
    if result["flagged_count"] > 0:
        sys.exit(1)


def cmd_analyze_cutoffs(args):
    from interpolators import generate_cutoff_recommendations_file

    thresholds = [float(v.strip()) for v in args.thresholds.split(",") if v.strip()]
    generate_cutoff_recommendations_file(
        output_file=args.output,
        n0_values=range(args.n0_min, args.n0_max + 1),
        l_max=args.l_max,
        qr_min=args.qr_min,
        qr_max=args.qr_max,
        n_probe=args.n_probe,
        thresholds=tuple(thresholds),
        num_workers=args.workers,
        integration_dps=args.integration_dps,
    )


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def _add_common_n0_l(p):
    p.add_argument("--n0-min", type=int, default=1, metavar="N", help="Minimum n0 value (default: 1)")
    p.add_argument("--n0-max", type=int, default=5, metavar="N", help="Maximum n0 value (default: 5)")
    p.add_argument("--l-max", type=int, default=10, metavar="L", help="Maximum l value (default: 10)")


def _add_qr_range(p):
    p.add_argument("--qr-min", type=float, default=0.0, help="Minimum qr value (default: 0.0)")
    p.add_argument("--qr-max", type=float, default=250.0, help="Maximum qr value (default: 250.0)")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="make_interpolators.py",
        description="Generate and test interpolator JSON files for I_integral.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ---- generate ----
    p_gen = sub.add_parser("generate", help="Build uniform-grid interpolators.json")
    p_gen.add_argument("-o", "--output", default="interpolators.json", help="Output file (default: interpolators.json)")
    _add_common_n0_l(p_gen)
    _add_qr_range(p_gen)
    p_gen.add_argument("--num-qr", type=int, default=4500, help="Number of uniform qr points (default: 4500)")
    p_gen.add_argument("--workers", type=int, default=None, metavar="N", help="Worker processes (default: all CPU cores)")
    p_gen.set_defaults(func=cmd_generate)

    # ---- generate-adaptive ----
    p_adp = sub.add_parser("generate-adaptive", help="Build adaptive PCHIP interpolator_adaptive.json")
    p_adp.add_argument("-o", "--output", default="interpolator_adaptive.json", help="Output file (default: interpolator_adaptive.json)")
    _add_common_n0_l(p_adp)
    _add_qr_range(p_adp)
    p_adp.add_argument("--rel-tol", type=float, default=1e-3, help="Refinement tolerance (default: 1e-3)")
    p_adp.add_argument("--abs-tol", type=float, default=1e-10, help="Absolute error tolerance floor (default: 1e-10)")
    p_adp.add_argument("--initial-points", type=int, default=24, help="Seed grid points (default: 24)")
    p_adp.add_argument("--max-points", type=int, default=8000, help="Maximum grid points per interpolator (default: 8000)")
    p_adp.add_argument("--workers", type=int, default=None, metavar="N", help="Worker processes (default: all CPU cores)")
    p_adp.add_argument("--auto-threshold", type=float, default=1e-10, help="Significance threshold relative to peak amplitude (default: 1e-10)")
    p_adp.add_argument("--low-q-switch", type=float, default=0.1, help="Use dedicated low-qr patch up to this qr (default: 0.1)")
    p_adp.add_argument("--low-q-anchor-points", type=int, default=32, help="Anchor points for low-qr patch (default: 32)")
    p_adp.add_argument("--low-q-change-threshold", type=float, default=1e-4, help="Auto low-qr switch trigger as fraction of peak amplitude (default: 1e-4)")
    p_adp.add_argument("--phase-tol", type=float, default=0.35, help="Maximum midpoint phase error in radians before refinement (default: 0.35)")
    p_adp.add_argument("--phase-amp-rel-floor", type=float, default=1e-3, help="Apply phase refinement only above this fraction of peak amplitude (default: 1e-3)")
    p_adp.add_argument("--min-interval-rel", type=float, default=1e-5, help="Do not refine intervals smaller than this fraction of effective qr_max (default: 1e-5)")
    p_adp.add_argument("--min-interval-abs", type=float, default=1e-6, help="Absolute minimum interval width for refinement (default: 1e-6)")
    p_adp.add_argument("--max-new-points-per-iter", type=int, default=128, help="Maximum new midpoint samples added per refinement iteration (default: 128)")
    p_adp.add_argument("--disable-auto-qr-max", action="store_true", help="Disable automatic effective qr_max detection")
    p_adp.add_argument("--disable-auto-low-q-switch", action="store_true", help="Disable automatic low-qr switch detection and use --low-q-switch directly")
    p_adp.add_argument("--disable-q-scaling", action="store_true", help="Disable q-scaled interpolation mode")
    p_adp.add_argument("--cutoff-file", type=str, default="interpolator_cutoffs.json", help="Per-channel cutoff JSON file (default: interpolator_cutoffs.json)")
    p_adp.add_argument("--cutoff-field", type=str, default="cutoffs_relaxed_5pct", help="Cutoff field to use (default: cutoffs_relaxed_5pct)")
    p_adp.add_argument("--cutoff-threshold", type=str, default="0.2", help="Threshold key from cutoff field (default: 0.2)")
    p_adp.add_argument("--cutoff-pad-factor", type=float, default=1.02, help="Multiply channel cutoff by this safety factor (default: 1.02)")
    p_adp.add_argument("--cutoff-min-qr-max", type=float, default=1.0, help="Minimum per-channel qr_max floor (default: 1.0)")
    p_adp.add_argument("--disable-auto-cutoffs", action="store_true", help="Disable automatic use of per-channel cutoff file")
    p_adp.add_argument("--zero-component-threshold", type=float, default=1e-2, help="RMS fraction below which a component is treated as exactly zero (default: 1e-2)")
    p_adp.add_argument("--use-tail-fit", dest="use_tail_fit", action="store_true", default=False, help="Enable damped oscillatory tail fit (default: disabled)")
    p_adp.add_argument("--no-use-tail-fit", dest="use_tail_fit", action="store_false", help="Disable damped oscillatory tail fit")
    p_adp.add_argument("--tail-fit-start-q", type=float, default=25.0, help="Start qr for damped tail fit (default: 25.0)")
    p_adp.add_argument("--fit-num-cycles", type=int, default=0, help="If >0, train fit using only first N minima/maxima cycles")
    p_adp.add_argument("--fit-scan-qr-max", type=float, default=120.0, help="Upper qr for initial cycle detection scan (default: 120)")
    p_adp.add_argument("--fit-scan-points", type=int, default=240, help="Points in initial cycle detection scan (default: 240)")
    p_adp.add_argument("--disable-point-cache", action="store_true", help="Disable persistent on-disk cache of evaluated I_integral points")
    p_adp.add_argument("--point-cache-dir", type=str, default="integral_point_cache", help="Directory for persistent per-channel I_integral point cache")
    p_adp.set_defaults(func=cmd_generate_adaptive)

    # ---- test ----
    p_test = sub.add_parser("test", help="Test uniform-grid interpolators against direct integration")
    p_test.add_argument("-i", "--input", default="interpolators.json", help="JSON file to test (default: interpolators.json)")
    _add_common_n0_l(p_test)
    p_test.add_argument("--threshold", type=float, default=10.0, help="Flag threshold in percent (default: 10.0)")
    p_test.add_argument("--max-interpolators", type=int, default=25, help="Cap on combinations tested (default: 25)")
    p_test.set_defaults(func=cmd_test)

    # ---- test-adaptive ----
    p_tadp = sub.add_parser("test-adaptive", help="Test adaptive PCHIP interpolators against direct integration")
    p_tadp.add_argument("-i", "--input", default="interpolator_adaptive.json", help="JSON file to test (default: interpolator_adaptive.json)")
    _add_common_n0_l(p_tadp)
    p_tadp.add_argument("--threshold", type=float, default=1.0, help="Flag threshold in percent (default: 1.0)")
    p_tadp.add_argument("--max-interpolators", type=int, default=25, help="Cap on combinations tested (default: 25)")
    p_tadp.add_argument("--plot-dir", type=str, default="adaptive_test_plots_main", help="Output directory for per-channel adaptive performance plots")
    p_tadp.add_argument("--plot-points", type=int, default=150, help="Samples per channel for adaptive test plots (default: 120)")
    p_tadp.add_argument("--plot-qr-max", type=float, default=None, help="Override test/plot maximum qr for all channels")
    p_tadp.add_argument("--workers", type=int, default=None, metavar="N", help="Worker processes (default: all CPU cores)")
    p_tadp.add_argument("--disable-point-cache", action="store_true", default=False, help="Disable on-disk integral point cache for test evaluations")
    p_tadp.add_argument("--point-cache-dir", type=str, default="integral_point_cache", help="Directory for per-channel integral point cache (default: integral_point_cache)")
    p_tadp.set_defaults(func=cmd_test_adaptive)

    # ---- analyze-cutoffs ----
    p_cut = sub.add_parser("analyze-cutoffs", help="Analyze direct-probe per-channel qr cutoff recommendations")
    p_cut.add_argument("-o", "--output", default="interpolator_cutoffs.json", help="Output file (default: interpolator_cutoffs.json)")
    _add_common_n0_l(p_cut)
    _add_qr_range(p_cut)
    p_cut.add_argument("--n-probe", type=int, default=120, help="Probe points per log/linear branch (default: 120)")
    p_cut.add_argument("--thresholds", type=str, default="1e-2,1e-3,1e-4,1e-5,1e-6", help="Comma-separated relative envelope thresholds")
    p_cut.add_argument("--integration-dps", type=int, default=35, help="mpmath decimal precision for direct probes (default: 35)")
    p_cut.add_argument("--workers", type=int, default=None, metavar="N", help="Worker processes (default: all CPU cores)")
    p_cut.set_defaults(func=cmd_analyze_cutoffs)

    return parser


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)
