from setuptools import setup, find_packages


###################################################################
NAME = 'msstitch'
PACKAGES = find_packages(where='src')
KEYWORDS = ['mass spectrometry', 'proteomics', 'processing']
CLASSIFIERS = [
    'Development Status :: 5 - Production/Stable',
    'Intended Audience :: Science/Research',
    'Natural Language :: English',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: Implementation :: CPython',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]
INSTALL_REQUIRES = ['numpy', 'lxml', 'biopython']
METADATA = {
    'version': '3.6',
    'title': 'msstitch',
    'description': 'MS proteomics post processing utilities',
    'uri': 'https://github.com/lehtiolab/msstitch',
    'author': 'Jorrit Boekel',
    'email': 'jorrit.boekel@scilifelab.se',
    'license': 'MIT',
    'copyright': 'Copyright (c) 2013 Jorrit Boekel',
}

CLI = {'console_scripts': ['msstitch=app.msstitch:main']}

###################################################################

from os import path
with open(path.join(path.abspath(path.dirname(__file__)), 'README.md'), encoding='utf-8') as fp:
    long_description = fp.read()


if __name__ == '__main__':
    setup(
        name=NAME,
        description=METADATA['description'],
        license=METADATA['license'],
        url=METADATA['uri'],
        version=METADATA['version'],
        author=METADATA['author'],
        author_email=METADATA['email'],
        maintainer=METADATA['author'],
        maintainer_email=METADATA['email'],
        keywords=KEYWORDS,
        packages=PACKAGES,
        package_dir={'': 'src'},
        long_description=long_description,
        long_description_content_type='text/markdown',
        classifiers=CLASSIFIERS,
        install_requires=INSTALL_REQUIRES,
        entry_points=CLI,
    )
