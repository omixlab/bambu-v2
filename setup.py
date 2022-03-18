from setuptools import setup, find_packages

setup(
    name="bambu-qsar",
    version='0.0.13',
    packages=find_packages(),
    author="Isadora Leitzke Guidotti, Frederico Schmitt Kremer",
    author_email="leitzke.gi@gmail.com, fred.s.kremer@gmail.com",
    description="bambu (bioassays model builder) is CLI tool to build QSAR models based on PubChem BioAssays datasets",
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    keywords="bioinformatics machine-learning data science drug discovery QSAR",
    entry_points = {'console_scripts':[
        'bambu-download   = bambu.download:main',
        'bambu-preprocess = bambu.preprocess:main',
        'bambu-train      = bambu.train:main',
        'bambu-predict    = bambu.predict:main',
        'bambu-validate   = bambu.validate:main',
        ]},
    install_requires = [
        requirement.strip('\n') for requirement in open("requirements.txt")
    ]
)
