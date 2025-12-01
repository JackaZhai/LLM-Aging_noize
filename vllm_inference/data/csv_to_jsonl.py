import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Read one or multiple Excel/CSV files and convert to JSONL. "
            "Supports your 5 xls tables (基本信息, 生活方式, 体格检查, 血液尿液指标, 既往病史)."
        )
    )

    parser.add_argument(
        "--input_files",
        nargs="+",
        required=True,
        help=(
            "List of input files (xls/xlsx/csv). "
            "If multiple, they will be merged by id column."
        ),
    )
    parser.add_argument(
        "--output_jsonl",
        required=True,
        help="Path to the output JSONL file.",
    )
    parser.add_argument(
        "--id_column",
        default="id",
        help="Name of the common ID column used to merge tables (default: id).",
    )
    parser.add_argument(
        "--encoding",
        default="utf-8",
        help="Encoding for CSV files (Excel will ignore this, default: utf-8)",
    )
    return parser.parse_args()


def load_table(path: Path, encoding: str = "utf-8") -> pd.DataFrame:
    """根据扩展名自动读取 xls/xlsx/csv 为 DataFrame。"""

    suffix = path.suffix.lower()
    if suffix in [".xls", ".xlsx"]:
        return pd.read_excel(path)
    if suffix == ".csv":
        return pd.read_csv(path, encoding=encoding)
    raise ValueError(f"Unsupported file type: {suffix} for {path}")


def merge_tables(paths, id_column: str, encoding: str = "utf-8") -> pd.DataFrame:
    """按 id_column 依次左连接合并多张表。"""

    dfs = []
    for p in paths:
        if not p.is_file():
            raise FileNotFoundError(f"Input file not found: {p}")
        df = load_table(p, encoding=encoding)
        if id_column not in df.columns:
            raise KeyError(f"Column '{id_column}' not found in file: {p}")
        dfs.append(df)

    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on=id_column, how="left")
    return merged


def df_to_jsonl(df: pd.DataFrame, output_jsonl: str) -> None:
    output_path = Path(output_jsonl)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # DataFrame 中的 NaN 转成 None，避免写出 "NaN"
    records = df.where(pd.notnull(df), None).to_dict(orient="records")

    with output_path.open("w", encoding="utf-8") as f:
        for row in records:
            json.dump(row, f, ensure_ascii=False)
            f.write("\n")


def main():
    args = parse_args()
    input_paths = [Path(p) for p in args.input_files]

    df_merged = merge_tables(input_paths, id_column=args.id_column, encoding=args.encoding)
    df_to_jsonl(df_merged, args.output_jsonl)


if __name__ == "__main__":
    main()
