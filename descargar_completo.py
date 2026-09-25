import os
from dotenv import load_dotenv
from huggingface_hub import snapshot_download

load_dotenv()

# Download challenge files into a local './data' folder
dataset_dir = snapshot_download(
    repo_id="SageBio/mva-hackathon-2026-data",
    repo_type="dataset",
    local_dir="./data",
    ignore_patterns=["README.md", ".gitattributes"],
    token=os.environ["HF_TOKEN"]
)

print(f"Dataset downloaded to: {dataset_dir}")
