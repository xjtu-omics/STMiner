from pathlib import Path

from setuptools import find_packages, setup


this_directory = Path(__file__).resolve().parent
readme_path = this_directory / 'README.md'
if not readme_path.exists():
    readme_path = this_directory / 'readme.md'
long_description = readme_path.read_text(encoding='utf-8')

setup(
    name='stminer',
    version='1.1.4',
    author='Peisen Sun',
    url='https://github.com/xjtu-omics/STMiner',
    license='GPL-3.0-only',
    description='Python package for spatial transcriptomics data analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author_email='sunpeisen@stu.xjtu.edu.cn',
    packages=find_packages(),
    python_requires='>=3.10',
    platforms=['Linux', 'Mac', 'Windows'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
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
                      'tifffile==2021.4.8',
                      'tqdm==4.65.0',
                      'umap_learn==0.5.3',
                      'scikit-misc']
)
