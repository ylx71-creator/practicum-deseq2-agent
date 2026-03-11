# practicum-deseq2-agent

RNA-seq differential expression analysis pipeline using **DESeq2**, exposed through a **FastAPI API service** and designed for future integration into an **agentic AI workflow**.

This project builds a modular pipeline for omics analysis and wraps the workflow into an API that can be called by automated agents.

---

# Project Overview

This repository implements an RNA-seq analysis pipeline that:

- preprocesses gene expression data  
- performs differential expression analysis using **PyDESeq2**  
- converts gene IDs to gene symbols  
- exposes the pipeline as a **REST API**

This project represents the **first stage of an agentic AI system for automated omics analysis.**

---

# Architecture

Current workflow:

User / Agent  
↓  
FastAPI API (api_app.py)  
↓  
DESeq2 Analysis Pipeline (step1_pipeline.py)  
↓  
Differential Expression Results  

The pipeline is exposed through an API endpoint so that it can be integrated into automated workflows.

---

# Project Structure


---

# API Endpoint

Run Differential Expression Analysis

POST /run_deseq2

Request body example:

{
  "threshold": 20,
  "padj_cutoff": 0.05,
  "lfc_cutoff": 1
}

Parameters:

threshold — clinical score threshold used to group samples  
padj_cutoff — adjusted p-value cutoff  
lfc_cutoff — log2 fold change cutoff  

---

# Running the API

Start the FastAPI server:

uvicorn api_app:app --reload

Open the interactive API documentation:

http://127.0.0.1:8000/docs

You can then run the pipeline through the API interface.

---



# Technologies Used

Python  
FastAPI  
PyDESeq2  
Pandas  
NumPy  

---

# Future Work

 

---

