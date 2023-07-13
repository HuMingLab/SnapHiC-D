from setuptools import setup, find_packages

with open("requirements.txt", "r") as f:
    requirements = f.read().strip().split("\n")

setup(
    name = 'SnapHiC-D',
    version = '0.1.0',
    author = 'Lindsay Lee, Hongyu Yu',
    author_email = 'hongyuyu@unc.edu',
    license = 'GNU',
    description = 'A computational pipeline to identify differential chromatin contacts from single cell Hi-C data',
    long_description = open('README.md').read(),
    long_description_content_type = "text/markdown",
    url = 'https://github.com/lindsayhrlee/SnapHiC-D',
    py_modules = ['run_SnapHiC_D', 'SnapHiC_D'],
    packages = find_packages(),
    install_requires = [requirements],
    python_requires = '>=3.6.8',
    entry_points = '''
        [console_scripts]
        SnapHiC-D=run_SnapHiC_D:main
    '''
)