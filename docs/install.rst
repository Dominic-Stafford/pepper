.. _installation:

Installation
============

It is recommended to use an isolated virtual environment with Pepper. An example environment setup for NAF or LXPLUS 
is included in the repository under `example/environment.sh`. To install Pepper as an editable package and in 
a virtual environment, use the ``install.sh`` script:

.. code-block:: bash

    git clone <repository url> pepper
    cd pepper
    chmod +x install.sh
    ./install.sh
    # Activate your newly created environment
    source environment.sh

This will update all dependencies to the latest version. 
Now ``pepper`` can be imported as any other python package from any location. 
Because of the ``--editable`` option, if you edit files inside your cloned ``pepper`` 
directory, the changes will be in effect already the next time you ``import pepper`` 
- i.e. there is no need to rebuild or reinstall the package after making changes.

Everytime you restart the shell, you need to activate the environment again. 
This just means running the `source environment.sh` command again. 
You can test your installation by running the included example script:

.. code-block:: bash

    cd example
    python3 -m pepper.runproc example_processor.py test_config.json
