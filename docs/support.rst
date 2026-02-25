Contributing
============

Feel free to submit merge requests to have your code included in this repository! Your code must comply with pep8. You can check this by running the following inside the Pepper directory:

.. code-block:: shell

   python3 -m pip install .[dev] --user
   python3 -m flake8



Writing Documentation
-----------------------

The docs are hosted on `cms-pepper.docs.cern.ch <https://cms-pepper.docs.cern.ch/>`__ and consists of two main parts:
   
1. The User Guides and Examples
2. The API Reference

The documentation is built with sphinx and stored in the ``docs`` directory of the repository. Most of the documentation is written in reStructuredText (reST) format, which is a lightweight markup language. You can find more information about reST syntax in the `reST documentation <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`__.
The API reference is autobuilt using sphinx-autodoc. To modify the API reference, you need to add docstrings to your code. These docstrings will be extracted by sphinx-autodoc and included in the API reference automatically.
The docs are automatically compiled and deployed to the documentation server on every push to the default branch.

Compiling the Documentation Locally
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To compile the documentation locally, you need to have a virtual environment installed with the ``.[docs] --user`` options.
You can install the required dependencies using pip:

.. code-block:: shell
 
   python3 -m pip install .[docs] --user


Once you have the dependencies installed, activate your virtual environment, and compile the documentation by running the 
following command inside the ``docs`` directory:

.. code-block:: shell

   make html



This will generate the HTML documentation in the ``docs/_build/html`` directory. 

.. hint:: 

   You can use the ``make clean`` command to remove the generated files and start the build process from scratch.
   
You can open the site using the python ``http.server`` module:

.. code-block:: shell

   python -m http.server --directory docs/_build/html


.. hint:: 

   In case you're developing on a remote server, you can still use port forwarding to access the documentation. For example, you can use the following command to forward the port:

   .. code-block:: shell

      ssh -L 8000:localhost:8888 user@remote-server

   Then, you can open the site using the python ``http.server`` module:

   .. code-block:: shell

      python -m http.server 8888 --directory docs/_build/html

   Then open http://localhost:8000 in your web browser.