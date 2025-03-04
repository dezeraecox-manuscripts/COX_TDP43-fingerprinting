import os
import zipfile
import shutil
import urllib.request as request
from contextlib import closing

output_folder = f'data/'

def download_file(url, output_folder):
    with closing(request.urlopen(url)) as r:
        with open(f'{output_folder}{url.split("/")[-1]}', 'wb') as f:
            shutil.copyfileobj(r, f)

    if url.split('.')[-1] == 'zip':
        filename = url.split("/")[-1]
        with zipfile.ZipFile(f'{output_folder}{filename}', 'r') as zip_ref:
            zip_ref.extractall(f'{output_folder}{filename.split(".")[-2]}/')
    
if __name__ == "__main__":

    datasets = {
        'microscopy-tissue': 'https://zenodo.org/records/14960397/files/microscopy-tissue.zip',
        'microscopy-RRMs': 'https://zenodo.org/records/14960397/files/microscopy-RRMs.zip',
        'proteomics': 'https://zenodo.org/records/14960397/files/proteomics.zip',
        'orthogonals': 'https://zenodo.org/records/14960397/files/orthogonals.zip',
    }

    for dataset, url in datasets.items():
        download_file(f'{url}', output_folder=f'{output_folder}{dataset}/') 
        
        filename = url.split('/')[-1]
        os.rename(f'{output_folder}{filename}', f'{output_folder}{dataset}.{filename}')
