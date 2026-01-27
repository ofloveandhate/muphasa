# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'muphasa'
copyright = 'The authors of muphasa'
author = 'The authors of muphasa'


import sys
import os

conf_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(conf_dir, "..", ".."))  # docs/source -> repo
src_python = os.path.join(repo_root, "src", "python")

sys.path.insert(0, src_python)

import importlib

try:
	mph = importlib.import_module("mph")
	print("mph imported from:", mph.__file__)
except:
	raise Exception('failed to import `mph`.  Try `pip install -e .` before running sphinx (`make html`), so that the Cython part of the library is built, before generating the documentation')

# # sys.path.insert(0,os.path.abspath('../../src'))
# sys.path.insert(0, os.path.abspath('../src/python/mph'))


# # Get the directory containing conf.py
# docs_dir = os.path.dirname(os.path.abspath(__file__))
# # Go up one level to the project root, then to src/python/mph
# project_root = os.path.dirname(docs_dir)
# src_dir = os.path.join(project_root, '..', 'src', 'python', 'mph')

# # Add the source directory to the path
# sys.path.insert(0, src_dir)

# print(f"Added to sys.path: {src_dir}")
# print(f"Source directory exists: {os.path.exists(src_dir)}")

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
	'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
    'sphinx.ext.autosectionlabel',
    'sphinxcontrib.bibtex']

autodoc_default_options = {
    "members": True,
    "undoc-members": True,      # include members without docstrings
    "inherited-members": True,  # optional
    "show-inheritance": True,   # optional
    "imported-members": True,
}

bibtex_bibfiles = ['./muphasa.bib']

# see http://www.sphinx-doc.org/en/stable/ext/autodoc.html#confval-autoclass_content
autoclass_content = 'both'

source_suffix = {'.rst': 'restructuredtext'}

templates_path = ['_templates']
exclude_patterns = []


# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
try:
    import git  #package gitpython
    repo = git.Repo(search_parent_directories=True) 
    last_commit = str(repo.head.commit) 
    version = last_commit[:7]
    release = last_commit # full version
except:
    last_commit = 'gitpython not installed'
    version = last_commit # The short X.Y version.
    release = version # The full version, including alpha/beta/rc tags.


# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
