from setuptools import setup, find_packages

setup(
    name='STMiner',
    version='1.0.0',
    author='PSSUN',
    url='https://github.com/PSSUN/STMiner',
    license='MIT License',
    description='Python package for spatial transcriptomics data analysis',
    author_email='sunpeisen@stu.xjtu.edu.cn',
    packages=find_packages(),
    platforms=['Linux', 'Mac', 'Windows'],
    keywords=['STMiner',
              'gaussian mixture models',
              'GMM',
              'hellinger distance',
              'Spatial transcriptomics']
)
