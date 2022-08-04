from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='commot',
    version='0.0.2',
    packages=find_packages(exclude=['tests*']),
    include_package_data=True,
    description='Cell-cell communications in space of transcriptomics data via collective optimal transport.',
    author='Zixuan Cang',
    author_email='cangzx@gmail.com',
    classifiers=[
    	"Programming Language :: Python :: 3.8",
    	"Programming Language :: Python :: 3.9",
    	"License :: OSI Approved :: MIT License",
    	"Operating System :: POSIX :: Linux",
    	"Operating System :: MacOS",
    	"Topic :: Scientific/Engineering",
    	"Topic :: Scientific/Engineering :: Bio-Informatics"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires = [
        "anndata>=0.7.6",
        "anndata2ri==1.0.6",
        "cython>=0.29.24",
        "karateclub>=1.2.2",
        "leidenalg>=0.8.4",
        "networkx>=2.5.1",
        "numpy>=1.19.5",
        "pandas>=1.4.1",
        "plotly>=5.3.1",
        "pot>=0.8.0",
        "pysal>=2.6.0",
        "python-igraph>=0.9.9",
        "python-louvain>=0.15",
        "rpy2==3.4.2",
        "scanpy==1.8.2",
        "scikit-learn==1.0.2"],
)
