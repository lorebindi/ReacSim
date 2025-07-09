from setuptools import setup, find_packages

setup(
    name="ReacSim",
    version="0.1.0",
    description="Gillespie simulator of chemical reactions system",
    author="Lorenzo Bindi, Elia Leonardi, Francesco Lucchesi",
    author_email="tuo@email.com",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "python-libsbml",
        "libroadrunner",
        "scipy"
    ],
    python_requires=">=3.7",
)