import setuptools


with open("README.md", mode="r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="euclid-observation-systematics",
    version="0.0.1",
    author="Pierluigi Monacoi <pierluigi.monaco@inaf.it>, Guilhem Lavaux <guilhem.lavaux@iap.fr>",
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
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires="<3",
)

