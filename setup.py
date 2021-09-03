from setuptools import setup, find_packages

setup(
    name='commot',
    version='0.1',
    packages=find_packages(exclude=['tests*']),
    description='Cell-cell communications in space of transcriptomics data via collective optimal transport.',
    author='Zixuan Cang',
    author_email='cangzx@gmail.com'
)
