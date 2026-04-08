from fastapi import FastAPI
from pydantic import BaseModel

from step1_pipeline import run_pipeline

app = FastAPI(title="DESeq2 Pipeline API")


class RunRequest(BaseModel):
    padj_cutoff: float = 0.05
    lfc_cutoff: float = 1.0


@app.get("/")
def home():
    return {"message": "DESeq2 Pipeline API is running.", "docs": "/docs"}


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/run_deseq2")
def run_deseq2(req: RunRequest):
    summary = run_pipeline(
        padj_cutoff=req.padj_cutoff,
        lfc_cutoff=req.lfc_cutoff,
    )
    return summary