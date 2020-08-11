import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ffbio-dnanto",
    version="0.0.1",
    author="Daniel Antonio Negrón",
    author_email="dnegron2@gmu.edu",
    description="flat-file sequence/database utils",
    long_description="scripts that process flat-file biological sequence databases",
    long_description_content_type="text/markdown",
    url="https://github.com/dnanto/ffbio",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
)
