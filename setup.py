from setuptools import setup, find_packages
from os.path import dirname, join
import io

with open('README.md', encoding='utf-8') as readme_file:
    readme = readme_file.read()

def get_version(relpath):
    """Read version info from a file without importing it."""
    for line in io.open(join(dirname(__file__), relpath), encoding="utf-8"):
        if "__version__" in line:
            # Expect version.py to define a dictionary with key "metacooc"
            version_dict = eval(line.split("=")[1])
            return version_dict["metacooc"]

setup(
    name='metacooc',
    version=get_version("metacooc/version.py"),
    description='Co-occurrence analysis of microorganisms in metagenomes',
    long_description=readme,
    long_description_content_type='text/markdown',
    url="https://github.com/yourusername/metacooc",  # update with your repository URL
    author='Your Name',  # update with your name
    license='GPL3+',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3',
    ],
    keywords="metagenomics co-occurrence bioinformatics",
    packages=find_packages(),
    data_files=[(".", ["README.md", "LICENCE.txt"])],
    include_package_data=True,
    install_requires=[
        'pandas>=1.0',
        'numpy>=1.15',
        'scipy>=1.0',
        'matplotlib>=3.0'
    ],
    entry_points={
        'console_scripts': [
            'metacooc = metacooc.__main__:main'
        ]
    },
)
