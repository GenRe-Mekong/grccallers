#!/usr/bin/env python3
"""
Top-level dispatcher for the grccallers CLI suite.

Usage:
    grccallers <tool> [options]
    grccallers --help
    grccallers k13 --help
"""
import argparse
import sys

from . import k13 as _k13


class _CustomFormatter(argparse.RawDescriptionHelpFormatter,
                       argparse.ArgumentDefaultsHelpFormatter):
    pass


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="grccallers",
        description=(
            "GRC mutation callers — variant-first callers for malaria "
            "drug-resistance genes.\n\n"
            "Run 'grccallers <tool> --help' for tool-specific options."
        ),
        formatter_class=_CustomFormatter,
    )
    sub = parser.add_subparsers(dest="command", metavar="TOOL", title="available tools")

    k13_parser = sub.add_parser(
        "k13",
        help="Kelch13 propeller-domain mutation caller.",
        description=_k13.DESCRIPTION,
        epilog=_k13.EPILOG,
        formatter_class=_CustomFormatter,
    )
    _k13.add_args(k13_parser)
    k13_parser.set_defaults(_run=_k13.run)

    return parser


def main(argv=None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args._run(args)


if __name__ == "__main__":
    main()
