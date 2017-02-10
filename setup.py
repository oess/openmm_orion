#!/usr/bin/env python
import re, ast, os
from os.path import relpath, join
from pip.req import parse_requirements
from setuptools import setup, find_packages


def get_reqs(reqs):
    return [str(ir.req) for ir in reqs]

try:
    install_reqs = get_reqs(parse_requirements("requirements.txt"))
except TypeError:
    from pip.download import PipSession
    install_reqs = get_reqs(
        parse_requirements("requirements.txt", session=PipSession())
    )


def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files

def get_version():
    _version_re = re.compile(r'__version__\s+=\s+(.*)')
    with open('OpenMMCubes/__init__.py', 'rb') as f:
        version = str(ast.literal_eval(_version_re.search(f.read().decode('utf-8')).group(1)))
        return version

setup(
    name="OpenMMCubes",
    version='0.1.2',
    #version=get_version(),
    packages=find_packages('examples', exclude=['tests*']),
    include_package_data=True,
    package_data={ 'examples' : find_package_data('examples/data', 'examples') },
    author="Christopher Bayly, Nathan M. Lim, John Chodera",
    author_email="bayly@eyesopen.com",
    description='Prepare complex for MD with OpenMM',
    install_requires=install_reqs,
    license='Other/Proprietary License',
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
    ]
)
