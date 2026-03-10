from pathlib import Path
from fastapi import FastAPI
from pydantic import BaseModel

from step1_pipeline import run_pipeline  # 你的 pipeline 文件名如果不同，这里改一下

app = FastAPI(title="DESeq2 Pipeline API")

class RunRequest(BaseModel):
    threshold: float = 20.0
    padj_cutoff: float = 0.05
    lfc_cutoff: float = 1.0

@app.post("/run_deseq2")
def run_deseq2(req: RunRequest):
    base_dir = Path(__file__).resolve().parent
    summary = run_pipeline(
        base_dir=base_dir,
        threshold=req.threshold,
        padj_cutoff=req.padj_cutoff,
        lfc_cutoff=req.lfc_cutoff,
    )
    return summary