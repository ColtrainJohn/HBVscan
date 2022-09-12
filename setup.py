from setuptools import setup, find_packages

setup(
    name='HBVisor',
    version='0.1',
    description='HBV sequence classifier',
    url='http://github.com/',
    author='Semen Selezov',
    author_email='selezov@gmail.com',
    license='Ok',
    packages=find_packages(),
    zip_safe=False,
    install_requires=[
       'biopython>=1.79',
       'pandas>=1.1.3'
       ],
)
