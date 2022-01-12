from setuptools import setup, find_packages

setup(
    name="bioassays-model-builder",
    version='0.0.5',
    packages=find_packages(),
    author="Isadora Leitzke Guidotti",
    author_email="leitzke.gi@gmail.com",
    description="bambu",
    keywords="bioinformatics machine-learning data science drug discovery QSAR",
    entry_points = {'console_scripts':[
        'bambu-download   = bambu.download:main',
        'bambu-preprocess = bambu.preprocess:main',
        'bambu-train      = bambu.train:main',
        'bambu-predict    = bambu.predict:main',
        ]},
    install_requires = [
        requirement.strip('\n') for requirement in open("requirements.txt")
    ]
)
