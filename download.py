import requests

def download_file_from_google_drive(file_id, destination):
    URL = "https://drive.google.com/uc?export=download"

    session = requests.Session()

    response = session.get(URL, params={'id': file_id}, stream=True)
    token = get_confirm_token(response)

    if token:
        params = {'id': file_id, 'confirm': token}
        response = session.get(URL, params=params, stream=True)

    save_response_content(response, destination)

def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value
    return None

def save_response_content(response, destination):
    CHUNK_SIZE = 32768

    with open(destination, "wb") as f:
        for chunk in response.iter_content(CHUNK_SIZE):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)

# Specify the IDs of the files you want to download
file_ids = {
    '18wWAgTyXgmgcS4mX6Pqk3D4WtVbhiGc9': "mse3-reload-hvg-level2-preint.h5ad",  
    '1gnkwNtzAqEK9-H3FXgPqgPxWmdq1Uugp': "Bashore_postQC_noCITE_noR.h5ad",  
}

# Specify the directory where you want to download the files
download_dir1 = 'data'  
download_dir2 = 'example_data' 

# Download the files
for file_id, file_name in file_ids.items():
    if file_name == "Bashore_postQC_noCITE_noR.h5ad":
        destination = f"{download_dir2}/{file_name}"
    elif file_name == "mse3-reload-hvg-level2-preint.h5ad":
        destination = f"{download_dir1}/{file_name}"
    download_file_from_google_drive(file_id, destination)
    print(f'Downloaded {file_name} to {destination}')