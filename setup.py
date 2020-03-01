from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='domino_hagai',
    version="0.1",
    author="Hagai Levi",
    author_email="hagai.levi.007@gmail.com",
    long_description=long_description,
    long_description_content_type="text/markdown",
    description='DOMINO: Discovery of Modules In Networks using Omic',
    url='https://github.com/hag007/DOMINO',
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    packages = find_packages(),
    install_requires=[
        'networkx==2.4',
        'matplotlib==3.1.3',
        'numpy==1.18.1',
        'scipy==1.4.1',
        'pandas==1.0.1',
        'pcst-fast==1.0.7',
        'statsmodels==0.11.0'],
    entry_points = {
        "console_scripts": [
            "domino=src.runner:main_domino",
            "slicer=src.runner:main_slicer",
        ]
    }



)