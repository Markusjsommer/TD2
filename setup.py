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
                      "pandas>=2.0.3",
                      "wheel"],
    extras_require={
        "dev": ["pytest>=7.0", "twine>=4.0.2", "pytest-cov>=4.0", "wheel"],
    },
    python_requires=">=3.8",
    entry_points={
        'console_scripts': [
            'transmark.LongOrfs = transmark.LongOrfs:main',
            'transmark.Predict = transmark.Predict:main',
        ],
    },

    include_package_data=True,
    package_data={'': ['tables/table_1.json',
                       'tables/table_2.json',
                       'tables/table_3.json',
                       'tables/table_4.json',
                       'tables/table_5.json',
                       'tables/table_6.json',
                       'tables/table_9.json',
                       'tables/table_10.json',
                       'tables/table_11.json',
                       'tables/table_12.json',
                       'tables/table_13.json',
                       'tables/table_14.json',
                       'tables/table_15.json',
                       'tables/table_16.json',
                       'tables/table_21.json',
                       'tables/table_22.json',
                       'tables/table_23.json',
                       'tables/table_24.json',
                       'tables/table_25.json',
                       'tables/table_26.json',
                       'tables/table_27.json',
                       'tables/table_28.json',
                       'tables/table_29.json',
                       'tables/table_30.json',
                       'tables/table_31.json',
                       'tables/table_33.json']},
    
)

