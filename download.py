import os
import gdown

def download_file(file_id, output_path):
    url = f"https://drive.google.com/uc?id={file_id}"
    gdown.download(url, output_path, quiet=False)

def main():
    # Create the output directory if it doesn't exist
    os.makedirs("output", exist_ok=True)

    # File IDs and their corresponding output paths
    files = [
        {"id": "1HFowENDJbQJTM9YmAbBxhvXwmBTdNCe5", "path": "example_data/Hu_subset.h5ad"},
        {"id": "1jJPHL2VEc4nW5z2pfgb8vJqC3lyt29iI", "path": "data/Big-Atlas-level12-log1p-hvg.h5ad"},
        {"id": "1iUcYvDPHZQnEKHmaMjxgesDmSrcKuBRI", "path": "data/Big-Atlas-level122-log1p-hvg.h5ad"}
    ]

    for file in files:
        # Create the directory if it doesn't exist
        os.makedirs(os.path.dirname(file["path"]), exist_ok=True)
        
        print(f"Downloading file to {file['path']}...")
        download_file(file["id"], file["path"])
        print(f"Download complete for {file['path']}")

if __name__ == "__main__":
    main()
