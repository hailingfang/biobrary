import setuptools

setuptools.setup(
    name="biobrary",
    version="0.0.1",
    author="Benjamin Fang",
    author_email="benjaminfang.ol@outlook.com",
    description="library for bioinformatics.",
    url="https://github.com/benjaminfang/biobrary",
    packages=setuptools.find_packages(
        where="./src"),
    python_requires=">=3.7",
    classifers=["Development Status :: 1 - Planning",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
    ],
)
