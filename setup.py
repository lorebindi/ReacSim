from setuptools import setup, find_packages

setup(
    name="ReacSim",
    version="0.1.0",
    description="Gillespie simulator of chemical reaction systems",
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
    entry_points={
        "console_scripts": [
            "reacsim=ReacSim.main:main"
        ]
    }
)