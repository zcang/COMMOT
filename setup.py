from setuptools import setup, find_packages

setup(
    name='commot',
    version='0.0.2',
    packages=find_packages(exclude=['tests*']),
    include_package_data=True,
    description='Cell-cell communications in space of transcriptomics data via collective optimal transport.',
    author='Zixuan Cang',
    author_email='cangzx@gmail.com'
)
