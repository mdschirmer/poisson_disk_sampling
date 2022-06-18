import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pds",
    version="0.1",
    author="Markus D. Schirmer",
    author_email='software@markus-schirmer.com',
    description='Poisson disk sampling approach for brain parcellations.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/mdschirmer/poisson_disk_sampling',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
)
