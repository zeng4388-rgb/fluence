from pathlib import Path

import numpy as np
from astropy.table import Table


try:
    BASE_DIR = Path(__file__).resolve().parent
except NameError:
    BASE_DIR = Path.cwd()
CATALOG_DIR = BASE_DIR / "catalog" / "catalog"
OUTPUT_CAT = BASE_DIR / "chime_frb_cat.txt"
OUTPUT_SVY = BASE_DIR / "chime_tel_svy.txt"


def as_text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore")
    return str(value)


def format_float(value, default=np.nan):
    if value is None:
        value = default
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def process_table(table, survey_name, file_handle):
    required = ["flux", "bc_width", "dm_fitb", "dm_exc_ne2001"]
    mask = np.ones(len(table), dtype=bool)
    for column in required:
        values = np.asarray(table[column], dtype=float)
        mask &= np.isfinite(values)

    filtered = table[mask]
    for row in filtered:
        name = as_text(row["tns_name"])
        gl = format_float(row["gl"])
        gb = format_float(row["gb"])
        s_val = format_float(row["flux"])
        s_err = format_float(row["flux_err"]) if "flux_err" in table.colnames else 0.0

        width_ms = format_float(row["bc_width"]) * 1e3
        fluence = format_float(row["fluence"])
        fluence_err = format_float(row["fluence_err"]) if "fluence_err" in table.colnames else 0.0

        dm = format_float(row["dm_fitb"])
        dm_ne2001 = format_float(row["dm_exc_ne2001"])
        dm_ymw16 = format_float(row["dm_exc_ymw16"])
        snr = format_float(row["bonsai_snr"])

        file_handle.write(
            f"{name:<12} {gl:<10.3f} {gb:<10.3f} {s_val:<8.3f} {s_err:<8.3f} {s_err:<8.3f} "
            f"{width_ms:<8.3f} {0.0:<8.3f} {0.0:<8.3f} {fluence:<8.3f} {fluence_err:<8.3f} {fluence_err:<8.3f} "
            f"{dm:<10.3f} {dm_ne2001:<12.3f} {dm_ymw16:<12.3f} {survey_name:<10} {snr:<8.3f} "
            f"{1.4:<8.1f} {50.0:<8.1f} {400.0:<10.1f} {2:<8} {9.0:<8.1f} {'CHIME':<15}\n"
        )


def main():
    tables = [
        (Table.read(CATALOG_DIR / "chimefrbcat1.fits"), "CHIME_1"),
        (Table.read(CATALOG_DIR / "chimefrbcat2.fits"), "CHIME_2"),
    ]

    with OUTPUT_CAT.open("w", encoding="utf-8", newline="") as fout:
        fout.write(
            f"{'FRBname':<12} {'gl':<10} {'gb':<10} {'S':<8} {'Seu':<8} {'Sel':<8} {'W':<8} {'Weu':<8} "
            f"{'Wel':<8} {'F':<8} {'Feu':<8} {'Fel':<8} {'DM':<10} {'DM_NE2001':<12} {'DM_YMW16':<12} "
            f"{'SURVEY':<10} {'SN':<8} {'Gain':<8} {'Tsys':<8} {'BW':<10} {'Npol':<8} {'SN0':<8} {'Ref':<15}\n"
        )

        for table, survey_name in tables:
            process_table(table, survey_name, fout)

    with OUTPUT_SVY.open("w", encoding="utf-8", newline="") as fsvy:
        fsvy.write("SURVEY     FOV        TIME       Gain    Tsys     Npol     SN0    BW\n")
        fsvy.write("CHIME_1    200        8760       1.4     50       2        9.0     400\n")
        fsvy.write("CHIME_2    200        17520      1.4     50       2        9.0     400\n")

    print("Conversion complete!")


if __name__ == "__main__":
    main()
