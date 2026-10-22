import argparse

from functions import ExtractSCENICPlusTables


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scplus_object")
    parser.add_argument("output_dir")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    ExtractSCENICPlusTables(args.scplus_object, args.output_dir, verbose=not args.quiet)
