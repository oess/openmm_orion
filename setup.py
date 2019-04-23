# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

# !/usr/bin/env python
import re
import ast


from setuptools import setup, find_packages

try:  # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError:  # for pip <= 9.0.3
    from pip.req import parse_requirements


def get_reqs(reqs):
    return [str(ir.req) for ir in reqs]


try:
    install_reqs = get_reqs(parse_requirements("requirements_dev.txt"))
except TypeError:
    try:  # For pip >= 10
        from pip._internal.download import PipSession
    except ImportError:  # For pip <= 9.0.3
        from pip.download import PipSession
    install_reqs = get_reqs(parse_requirements("requirements_dev.txt", session=PipSession()))


def get_version():
    _version_re = re.compile(r'__version__\s+=\s+(.*)')
    with open('MDOrion/__init__.py', 'rb') as f:
        version = str(ast.literal_eval(_version_re.search(f.read().decode('utf-8')).group(1)))
        return version


setup(
    name="OpenEye MD Floes",
    version=get_version(),
    packages=find_packages(exclude=['tests*']),
    include_package_data=True,
    author="Gaetano Calabro, Christopher Bayly",
    author_email="gcalabro@eyesopen.com",
    description='Orion cubes to perform MD and MD analysis',
    install_requires=install_reqs,
    license='Other/Proprietary License',
    url='https://github.com/oess/openmm_orion',
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
    ],
    entry_points='''
        [console_scripts]
        mdocli=MDOcli.command_line:main
    ''',
    zip_safe=False
)
