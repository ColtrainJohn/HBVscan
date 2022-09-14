from setuptools import setup, find_packages
from subprocess import run
import pathlib

root = str(pathlib.Path(__file__).parent.resolve())
run(['chmod', '+x', f'{root}/hbv/data/blastn'])

setup(
    name='hbv',
    version='0.0.2',
    description='HBV sequence classifier',
    author='Semen Selezov',
    author_email='1selezov1@gmail.com',
    packages=['hbv'],
    package_data={'hbv' : ['hbv*', 'blastn']},
    include_package_data=True,
    install_requires=[
       'biopython>=1.79',
       'pandas>=1.1.3'
    ]
)
