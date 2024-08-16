import os
import gdown

def download_file(file_id, output_path):
    url = f"https://drive.google.com/uc?id={file_id}"
    gdown.download(url, output_path, quiet=False)

def main():
    # File IDs and their corresponding output paths
    files = [
        {"id": "18wWAgTyXgmgcS4mX6Pqk3D4WtVbhiGc9", "path": "data/mse3-reload-hvg-level2-preint.h5ad"},
        {"id": "1gnkwNtzAqEK9-H3FXgPqgPxWmdq1Uugp", "path": "example_data/Bashore_postQC_noCITE_noR.h5ad"}
    ]

    for file in files:
        # Create the directory if it doesn't exist
        os.makedirs(os.path.dirname(file["path"]), exist_ok=True)
        
        print(f"Downloading file to {file['path']}...")
        download_file(file["id"], file["path"])
        print(f"Download complete for {file['path']}")

if __name__ == "__main__":
    main()
