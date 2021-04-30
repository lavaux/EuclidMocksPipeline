import setuptools


with open("README.md", mode="r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="euclid-observation-systematics",
    version="0.0.1",
    author="Pierluigi Monacoi <pierluigi.monaco@inaf.it>, Tiago Castro, Guilhem Lavaux <guilhem.lavaux@iap.fr>",
    author_email="pierluigi.monaco@inaf.it",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    package_dir={"euclid_obssys": "src/euclid_obssys"},
    scripts=[
        "scripts/euclidGenerateInput.py",
        "scripts/createSDHOD_Catalog.py
    ],
    package_data={"euclid_obssys":["templates/input.py"]},
    packages=setuptools.find_packages(where="src"),
    python_requires="<3",
)

