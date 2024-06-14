from setuptools import find_packages, setup

with open("README.md", "r") as f:
    long_description = f.read()
    
setup(
    name="transmark",
    version="0.1.0",
    description="A tool to find protein coding ORFs",
    packages=find_packages(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Markusjsommer/TD2",
    author="markus",
    author_email = "markusjsommer@gmail.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.11",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
    ],
    install_requires=["psauron>=1.0.2",
                      "typing-extensions>=4.9.0",
                      "tqdm>=4.66.1",
                      "scipy>=1.10.1",
                      "biopython>=1.83",
                      "numpy>=1.24.4",
                      "pandas>=2.0.3"],
    extras_require={
        "dev": ["pytest>=7.0", "twine>=4.0.2", "pytest-cov>=4.0", "wheel"],
    },
    python_requires=">=3.8",
    entry_points={
        'console_scripts': [
            'transmark = app.transmark:main',
        ],
    },
    #entry_points={
    #        'console_scripts': [
    #        'transmark = app.transmark:main',
    #    ],
    #},
    #include_package_data=True,
    #package_data={'': ['/data/model_state_dict.pt']},
    
)

