#!/usr/bin/env python


from distutils.core import setup
import os


long_description='Please read the documentation on Github.'
if os.path.exists('README.md'):
    long_description = open('README.md').read()


setup(name='pyDigHolo',
    version='0.1',
    description="Python wrapping for Joel's Carpenter digHolo C++ module: https://github.com/joelacarpenter/digHolo/",
    author='Sebastien M. Popoff',
    author_email='sebastien.popoff@espci.psl.eu',
    url='https://www.wavefrontshaping.net',
    license = 'MIT',
    packages = ['pyDigHolo'],
    long_description=long_description,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: MIT License",
    ],
    install_requires=[
          'numpy',
      ],
    )