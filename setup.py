from setuptools import setup, find_packages

setup(
    name="ReacSim",  # Cambia con il nome del tuo progetto
    version="0.1.0",
    description="simulatore Gillespie di reazioni chimiche",
    author="Lorenzo Bindi, Elia Leonardi, Francesco Lucchesi",
    author_email="tuo@email.com",
    packages=find_packages(),  # Trova tutti i package nella cartella
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "python-libsbml",
        "libroadrunner",
        "scipy"
        # Aggiungi qui tutte le librerie di cui hai bisogno
        # es: "roadrunner", "chempy", "pint", ecc.
    ],
    python_requires=">=3.7",
)