import pandas as pd


def append_results(results, csv_path):
    write_header = not csv_path.exists()
    pd.DataFrame(results).to_csv(
        csv_path,
        mode="w" if write_header else "a",
        header=write_header,
        index=False,
    )
