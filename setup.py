from setuptools import setup, find_packages

setup(
    name="bambu-qsar",
    version='0.0.5',
    packages=find_packages(),
    author="Isadora Leitzke Guidotti",
    author_email="leitzke.gi@gmail.com",
    description="bambu",
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
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
