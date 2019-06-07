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

import os

import copy

import shutil

from invoke import task, run

from json import loads, dump

from glob import iglob

from importlib.machinery import SourceFileLoader

import importlib

import MDOrion

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FLOES_DIR = os.path.join(PACKAGE_DIR, "floes")

@task
def flake8(ctx):
    run("flake8 MDOrion")


@task
def test_cubes(ctx):
    """
    run cube tests
    """

    clean(ctx)
    run("py.test -s -v -m 'local'")


@task
def test_floes(ctx, test="all"):
    """
    run tests
    """
    clean(ctx)
    if test == "all":
        run("py.test -s -v ./tests ")
    else:
        run("""py.test -s -v -m "{}" ./tests """.format(test))


@task
def test_orion(ctx, profile="", test="orion"):
    """
    run tests
    """
    clean(ctx)

    if profile is "":
        if "ORION_PROFILE" in os.environ:
            profile = os.getenv("ORION_PROFILE")
        else:
            profile = 'default'

    print("Using Orion Profile: {}".format(profile))

    if test == "all":
        run("ORION_PROFILE={} py.test -s -v --orion ./tests".format(profile))
    else:
        run("""ORION_PROFILE={} py.test -s -v --orion --no-cleanup -m "{}" ./tests""".format(profile, test))


@task
def setversion(ctx, new_version):
    """
    Set the package version
    """

    # clean(ctx)

    fn = os.path.join("./MDOrion", "__init__.py")

    with open(fn, "r") as f:
        lines = f.readlines()

    lines = ["__version__ = '{}'\n".format(new_version) if '__version__' in line else line for line in lines]

    with open(fn, "w") as f:
        f.writelines(lines)

    spec = loads(open('manifest.json', 'r').read())

    importlib.reload(MDOrion)

    spec['version'] = MDOrion.__version__
    dump(spec, open('manifest.json', 'w'), sort_keys=True, indent=4)


@task
def release(ctx):
    """
    Create a package for the distribution. All the floes where
    the release variable is set to True are included in the package
    """

    # clean(ctx)

    # Un-wanted Package list
    pkg_un = ['OpenEye-Artemis', 'OpenEye-floe-pkg-tools']

    with open(os.path.join(PACKAGE_DIR, "requirements_dev.txt"), "r") as f:
        requirements_lines = f.readlines()

    original_requirements = copy.deepcopy(requirements_lines)

    for idx in range(0, len(requirements_lines)):
        for pkg in pkg_un:
            if pkg in requirements_lines[idx]:
                requirements_lines[idx] = "# " + requirements_lines[idx]

    with open("requirements_dev.txt", "w") as f:
        f.writelines(requirements_lines)

    run("python setup.py sdist")

    with open("requirements_dev.txt", "w") as f:
        f.writelines(original_requirements)


@task
def docs(ctx):
    clean_docs(ctx)
    curdir = os.getcwd()
    run('cube_doc MDOrion docs/source')
    run('floe_doc "MDOrion Floes" floes docs/source')
    os.chdir('docs')
    run("make html")
    os.chdir(curdir)

# @task
# def make_package_docs(ctx):
#     docs(ctx)
#     run("make_package_docs docs/build/html/cubes")


@task
def clean(ctx):
    """
    Clean up doc and package builds
    """
    clean_pyc(ctx)
    clean_docs(ctx)
    clean_pycache(ctx)
    shutil.rmtree("dist", ignore_errors=True)
    shutil.rmtree("build", ignore_errors=True)
    egg_path = "{}.egg-info".format("OpenEye_MD_Floes".replace("-", "_"))
    if os.path.isfile(egg_path):
        os.remove(egg_path)
    elif os.path.isdir(egg_path):
        shutil.rmtree(egg_path)
    shutil.rmtree(".pytest_cache", ignore_errors=True)


@task
def clean_pyc(ctx):
    """
    cleans out .pyc files
    """
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".pyc"):
                filename = os.path.join(root, file)
                if os.path.exists(filename):
                    os.unlink(filename)


@task
def clean_pycache(ctx):
    """
    cleans out __pycache__ dirs
    """
    for dirpath, dirs, files in os.walk(os.getcwd()):
        for dir in dirs:
            if dir == '__pycache__':
                del_dir = os.path.join(dirpath, dir)
                shutil.rmtree(del_dir)


@task
def clean_docs(ctx):
    doc_dir = "docs/build/html"
    _clean_out_dir(doc_dir)

    if os.path.isdir("docs/build/doctrees"):
        shutil.rmtree("docs/build/doctrees")


def _clean_out_dir(dir_path):
    if os.path.isdir(dir_path):
        for the_file in os.listdir(dir_path):
            file_path = os.path.join(dir_path, the_file)
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(e)
