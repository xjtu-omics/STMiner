from setuptools import setup, find_packages
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='STMiner',
    version='1.0.9',
    author='Peisen Sun',
    url='https://github.com/xjtu-omics/STMiner',
    license='GPL-3.0 license',
    description='Python package for spatial transcriptomics data analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author_email='sunpeisen@stu.xjtu.edu.cn',
    packages=find_packages(),
    platforms=['Linux', 'Mac', 'Windows'],
    keywords=['STMiner',
                      'bioinformatics',
                      'GMM',
                      'Hellinger distance',
                      'Spatial transcriptomics'],
    install_requires=['anndata==0.10.7',
                              'bioservices==1.11.2',
                              'matplotlib==3.7.2',
                              'networkx==3.1',
                              'numba==0.57.0',
                              'numpy==1.24.3',
                              'pandas==2.0.3',
                              'Pillow==10.3.0',
                              'plotly==5.9.0',
                              'POT==0.9.1',
                              'scanpy==1.10.1',
                              'scikit_learn==1.3.0',
                              'scipy==1.11.1',
                              'seaborn==0.13.2',
                              'setuptools==68.0.0',
                              'tifffile==2021.4.8',
                              'tqdm==4.65.0',
                              'umap_learn==0.5.3']
)
