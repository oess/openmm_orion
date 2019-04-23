# The OpenEye Sphinx Theme

This is the OpenEye theme for documentation. Set up this repo as a submodule within the docs directory of your project.

Navigate into the docs directory of your project and run the following command:

    git submodule add ssh://git@git.eyesopen.com/common/doc-theme.git doc-theme


For example, in [OENotebook](https://git.eyesopen.com/projects/OTHER/repos/oenotebook) this is the structure of the docs directory.

```
docs/
├── Makefile
├── _build
├── _static
├── conf.py
├── doc-theme
├── index.rst
├── make.bat
└── release_notes.rst
```

#### The following changes need to be made to `conf.py` within the docs directory:

* Add the repo directory to the `html_theme_path` variable


    html_theme_path = ['./doc-theme']

* Add the _static subdirectory to `html_static_path`


    html_static_path = ['_static', './doc-theme/_static']


* This theme is an extension of the sphinx rtd theme


    html_theme = 'sphinx_rtd_theme'
