#!/usr/bin/env python

import argparse
import csv
import logging
import sys
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    def __init__(
        self,
        condition_col="condition",
        treatment_col="treatment",
        control_col="control",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self._condition_col = condition_col
        self._treatment_col = treatment_col
        self._control_col = control_col

    def validate(self, row, condition_list, samplesheet_headers):
        self._check_control_col(row, condition_list)
        self._check_treatment_col(row, condition_list)

    def create_condition_list(self, reader1):
        samplesheet_headers = reader1.fieldnames
        condition_list = list()
        for row in reader1:
            condition_column = row[self._condition_col]
            condition_list.append(condition_column)
            condition_list = list(set(condition_list))
        return condition_list, samplesheet_headers

    def _check_control_col(self, row, condition_list):
        if row[self._control_col] not in condition_list:
            raise AssertionError(
                f"{row[self._control_col]} in contrastsheet does not match with any value in the condition column in samplesheet \n"
                f"Please ensure that the values of control column are consistent with the values in the condition column (In samplesheet)"
            )

    def _check_treatment_col(self, row, condition_list):
        if row[self._treatment_col] not in condition_list:
            raise AssertionError(
                f"{row[self._treatment_col]} in contrastsheet does not match with any value in the condition column in samplesheet. \n"
                f"Please ensure that the values of treatment column are consistent with the values in the condition column (In samplesheet)"
            )

def check_sheets(samplesheet, contrastsheet):
    with samplesheet.open(newline="") as in_samplesheet, contrastsheet.open(newline="") as in_contrastsheet:
        reader1 = csv.DictReader(in_samplesheet)
        reader2 = csv.DictReader(in_contrastsheet)
        checker = RowChecker()
        # check per row
        condition_list, samplesheet_headers = checker.create_condition_list(reader1)
        for i, row in enumerate(reader2):
            try:
                checker.validate(row, condition_list, samplesheet_headers)
            except AssertionError as error:
                logger.critical(f"{str(error)} on line {i + 2}")
                sys.exit(1)
        # check block col
        try:
            checker.check_block_unique(reader2)
        except AssertionError as error:
            logger.critical(f"{str(error)}.")
            sys.exit(1)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate whether the samplesheet and contrastsheet match.",
        epilog="Example: python samplesheet_contrastsheet_checker.py samplesheet.csv contrastsheet.csv",
    )
    parser.add_argument(
        "samplesheet",
        metavar="SAMPLESHEET",
        type=Path,
        help="Tabular input samplesheet in CSV format.",
    )
    parser.add_argument(
        "contrastsheet",
        metavar="CONTRASTSHEET",
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
    if not args.samplesheet.is_file():
        logger.error(f"The given input file {args.samplesheet} was not found!")
        sys.exit(2)
    if not args.contrastsheet.is_file():
        logger.error(f"The given input file {args.contrastsheet} was not found!")
        sys.exit(2)
    check_sheets(args.samplesheet, args.contrastsheet)


if __name__ == "__main__":
    sys.exit(main())
