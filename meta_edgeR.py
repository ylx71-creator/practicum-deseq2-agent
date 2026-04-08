import re
import pandas as pd
from pathlib import Path

base = Path(__file__).resolve().parent

# 读取文件
lines = (base/"data/GSE164641_series_matrix.txt").read_text(
    encoding="utf-8", errors="ignore"
).splitlines()

# 提取 sample_id
gsm = re.findall(r"GSM\d+", next(l for l in lines if l.startswith("!Sample_geo_accession")))

# 提取 High / Average
cat_line = next(l for l in lines if "risk category" in l.lower())
cats = re.findall(r"risk category: ([A-Za-z]+)", cat_line)

# 对齐 + 保存
meta = pd.DataFrame({"sample_id": gsm[:len(cats)], "condition": cats})
meta.to_csv(base/"output/metadata_for_edgeR.tsv", sep="\t", index=False)

print(meta.head())

