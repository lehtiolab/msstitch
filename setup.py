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
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: Implementation :: CPython',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]
INSTALL_REQUIRES = ['numpy', 'lxml', 'biopython', 'pyyaml']
METADATA = {
    'version': '2.10',
    'title': 'msstitch',
    'description': 'MS proteomics post processing utilities',
    'uri': 'https://github.com/glormph/msstitch',
    'author': 'Jorrit Boekel',
    'email': 'jorrit.boekel@scilifelab.se',
    'license': 'MIT',
    'copyright': 'Copyright (c) 2013 Jorrit Boekel',
}

CLI = {'console_scripts': ['msspercolator=app.pycolator:main',
                           'msslookup=app.mslookup:main',
                           'msspsmtable=app.mzidtsv:main',
                           'msspeptable=app.peptable:main',
                           'mssprottable=app.prottable:main',
                           ]}
###################################################################

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst', format='md')
except(IOError, ImportError):
    long_description = open('README.md').read()


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
        classifiers=CLASSIFIERS,
        install_requires=INSTALL_REQUIRES,
        entry_points=CLI,
    )
