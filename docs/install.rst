.. _installation:

Installation
============

It is recommended to use an isolated virtual environment with Pepper. An example environment setup for NAF or LXPLUS 
is included in the repository under `example/environment.sh`. To install Pepper as an editable package, follow these steps:

.. code-block:: bash

    git clone <repository url> pepper
    cd pepper
    source example/environment.sh
    # We recomend installing inside a virtual environment to avoid version conflicts with your other projects
    python3 -m venv .venv
    source .venv/bin/activate
    python3 -m pip install --upgrade --upgrade-strategy eager --editable .

This will update all dependencies to the latest version. 
Now ``pepper`` can be imported as any other python package from any location. 
Because of the ``--editable`` option, if you edit files inside your cloned ``pepper`` 
directory, the changes will be in effect already the next time you ``import pepper`` 
- i.e. there is no need to rebuild or reinstall the package after making changes.

It is recommended that you add a command sourcing your virtual environment to 
your environment setup script, i.e. by adding

.. code-block:: bash

    source /path/to/pepper/.venv/bin/activate

In ``example/environment.sh``, there is a placeholder line for this purpose.