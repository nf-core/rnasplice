#!/usr/bin/env python
# Adapted from nf-core/rnaseq check_samplesheet.py script.
# Please see nf-core/rnaseq for author and license.

"""Provide a command line tool to validate and transform tabular contrastsheets."""


import argparse
import csv
import logging
import re
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    def __init__(
        self,
        contrast_col="contrast",
        treatment_col="treatment",
        control_col="control",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            contrast_col (str): The name of the column that contains the contrast name (default "contrast").
            treatment_col (str): The name of the column that contains the treatment condition (default "treatment").
            second_col (str): The name of the column that contains the control condition (default "control").
        """
        super().__init__(**kwargs)
        self._contrast_col = contrast_col
        self._treatment_col = treatment_col
        self._control_col = control_col
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_contrast(row)
        self._validate_treatment(row)
        self._validate_control(row)
        self.modified.append(row)

    def _validate_contrast(self, row):
        """Assert that the contrast name exists and convert spaces to underscores."""
        if len(row[self._contrast_col]) <= 0:
            raise AssertionError("Contrast name is required.")
        # Sanitize contrasts slightly.
        row[self._contrast_col] = row[self._contrast_col].replace(" ", "_")

    def _validate_treatment(self, row):
        """Assert that the treatment name exists and has a syntactically valid name."""
        if len(row[self._treatment_col]) <= 0:
            raise AssertionError("Treatment name is required.")
        # Sanity check condition.
        self._validate_condition_value(row[self._treatment_col])

    def _validate_control(self, row):
        """Assert that the control name exists and has a syntactically valid name."""
        if len(row[self._control_col]) <= 0:
            raise AssertionError("Control name is required.")
        # Sanity check condition.
        self._validate_condition_value(row[self._control_col])

    def _validate_condition_value(self, condition):
        regex = "^(([A-Za-z]|[.][._A-Za-z])[._A-Za-z0-9]*)|[.]$"
        assert bool(re.search(regex, condition)), (
            f"The condition column has an invalid name: {condition}\n"
            f"A syntactically valid name consists of letters, numbers and the dot or underline characters and starts with a letter or the dot not followed by a number."
        )


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


# def sniff_format(handle):
#     """
#     Detect the tabular format.

#     Args:
#         handle (text file): A handle to a `text file`_ object. The read position is
#         expected to be at the beginning (index 0).

#     Returns:
#         csv.Dialect: The detected tabular format.

#     .. _text file:
#         https://docs.python.org/3/glossary.html#term-text-file

#     """
#     peek = read_head(handle)
#     handle.seek(0)
#     sniffer = csv.Sniffer()
#     if not sniffer.has_header(peek):
#         logger.critical("The given contrast sheet does not appear to contain a header.")
#         sys.exit(1)
#     dialect = sniffer.sniff(peek)
#     return dialect


def check_contrastsheet(file_in, file_out):
    """
    Check that the tabular contrastsheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row.

    Args:
        file_in (pathlib.Path): The given tabular contrastsheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed contrastsheet should
            be created; always in CSV format.

    Example:
        This function checks that the contrastsheet follows the following structure:

        contrast,treatment,control
        A-B,A,B
        A-C,A,C
        B-C,B,C

    """
    required_columns = {"contrast", "treatment", "control"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle)
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The contrast sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)

    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular contrastsheet.",
        epilog="Example: python check_contrastsheet.py contrastsheet.csv contrastsheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input contrastsheet in CSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output contrastsheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_contrastsheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
