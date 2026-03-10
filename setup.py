from setuptools import setup


setup(
    name="bigwig-compressor",
    version="0.1.0",
    description="Attenuate high-expression genes in RNA bigWig tracks while preserving bigWig output.",
    py_modules=["attenuate_rna_bigwig"],
    python_requires=">=3.10",
    install_requires=[
        "matplotlib==3.10.8",
        "numpy",
        "pyBigWig==0.3.25",
    ],
    entry_points={
        "console_scripts": [
            "attenuate-rna-bigwig=attenuate_rna_bigwig:main",
        ]
    },
)
