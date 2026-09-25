import os
from dotenv import load_dotenv
from huggingface_hub import hf_hub_download

load_dotenv()

# Download a specific file into a local './data' folder
file_path = hf_hub_download(
    repo_id="SageBio/mva-hackathon-2026-data",
    filename="WGS_EX2312012_HGWCNDSX7.vcf.gz",  # Replace with target file
    repo_type="dataset",
    local_dir="./data",
    token=os.environ["HF_TOKEN"]
)

print(f"File downloaded to: {file_path}")
