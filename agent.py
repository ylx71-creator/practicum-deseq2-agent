import requests


API_URL = "http://127.0.0.1:8000/run_deseq2"


def parse_user_request(user_request: str) -> dict:
    """
    Very simple rule-based parser.
    Later this can be replaced by an LLM.
    """
    user_request = user_request.lower()

    params = {
        "threshold": 20.0,
        "padj_cutoff": 0.05,
        "lfc_cutoff": 1.0,
    }

    if "25" in user_request:
        params["threshold"] = 25.0
    if "0.01" in user_request:
        params["padj_cutoff"] = 0.01
    if "2.0" in user_request or "2 " in user_request:
        params["lfc_cutoff"] = 2.0

    return params


def should_run_deseq2(user_request: str) -> bool:
    keywords = ["deseq2", "differential", "gene", "rna", "expression"]
    text = user_request.lower()
    return any(k in text for k in keywords)


def call_deseq2_api(params: dict) -> dict:
    response = requests.post(API_URL, json=params)
    response.raise_for_status()
    return response.json()


def summarize_result(result: dict) -> str:
    n_sig = result.get("n_sig", "unknown")
    n_genes = result.get("n_genes", "unknown")
    n_used = result.get("n_used", "unknown")

    return (
        f"The DESeq2 pipeline finished successfully. "
        f"It analyzed {n_used} samples and {n_genes} genes, "
        f"and identified {n_sig} significant genes."
    )


def main():
    user_request = input("Enter your request: ").strip()

    if not should_run_deseq2(user_request):
        print("This request does not appear to require the DESeq2 pipeline.")
        return

    print("\n[Agent] Request recognized as a DESeq2 task.")
    params = parse_user_request(user_request)
    print("[Agent] Parsed parameters:", params)

    try:
        result = call_deseq2_api(params)
        print("\n[Agent] API returned:")
        print(result)

        summary = summarize_result(result)
        print("\n[Agent] Summary:")
        print(summary)

    except requests.exceptions.RequestException as e:
        print("\n[Agent] Failed to call API:")
        print(str(e))


if __name__ == "__main__":
    main()