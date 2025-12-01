import argparse
import csv
import json
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Convert CSV file to JSONL format.")
    parser.add_argument("--input_csv", required=True, help="Path to the input CSV file.")
    parser.add_argument("--output_jsonl", required=True, help="Path to the output JSONL file.")
    parser.add_argument(
        "--encoding",
        default="utf-8",
        help="File encoding for input and output (default: utf-8)",
    )
    return parser.parse_args()


def csv_to_jsonl(input_csv: str, output_jsonl: str, encoding: str = "utf-8") -> None:
    input_path = Path(input_csv)
    output_path = Path(output_jsonl)

    if not input_path.is_file():
        raise FileNotFoundError(f"Input CSV file not found: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding=encoding, newline="") as csv_file, output_path.open(
        "w", encoding=encoding
    ) as jsonl_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            # 直接将每一行字典写为一行 JSON
            json.dump(row, jsonl_file, ensure_ascii=False)
            jsonl_file.write("\n")


if __name__ == "__main__":
    args = parse_args()
    csv_to_jsonl(args.input_csv, args.output_jsonl, args.encoding)
