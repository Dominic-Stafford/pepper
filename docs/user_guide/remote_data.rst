Accessing Remote Data
=====================
The data and MC NanoAODs are stored on various CERN Tier 2 computing clusters worldwide and Pepper is able to access data sets from remote servers using XRootD. You should specify ``"file_mode": "local+xrootd"`` and ``"xrootddomain": "xrootd-cms.infn.it"`` in your config to enable it. 
There are four requirements to get it working (also in conjunction with HTCondor):
 - xrootd must be installed. You can check with ``python3 -m pip show xrootd``
 - The CMS Grid environment needs to be sourced (also inside the HTCondor job). This is done by the [example environment script](example/environment.sh).
 - The environment variable ``X509_USER_PROXY`` needs to be set to a file path accessible by Condor (/tmp/ is not). As above this is also needed inside the Condor job and is cone by the script.
 - A VOMS proxy needs to be created at the path pointed to by ``X509_USER_PROXY``. To do this please once run: ``voms-proxy-init --voms cms --out $X509_USER_PROXY``.

If you get the error 'sslv3 alert certificate expired', please run the voms-proxy-init command again.

File Transfer Request
----------------------
Accessing files on a remote server using XRootD is much slower than accessing local files. If you want to speed up the 
data processing it is recommended to request files to be transferred to a local storage location. Pepper provides a script to
request a file transfer for all data and MC files defined in a single config file. It requires access to ``Rucio`` which 
is obtained through the following commands:

.. code-block:: bash

    source /cvmfs/cms.cern.ch/rucio/setup-py3.sh
    export RUCIO_ACCOUNT=<YOUR_CERN_USERNAME>

Afterwards, the script is invoked like this:

.. code-block:: bash

    python3 scripts/rucio_create_rules.py <path/to/config.json> <site>

where site could for example could be ``T2_DE_DESY`` for a transfer to DESY.

.. note::
   
   The script will create Rucio requests for all data and MC files defined in the config file to be transferred to the specified site.

After the script has been run, requests for all file transfers will be created. These need to be approved by a site administrator before the transfer starts.
Typically, this process will take a couple of working-days depending on the site. It is also worth noting that the requested rules will have a lifetime 
of maximally six months. After this time, the files will be deleted from the local storage unless a new request is made.

.. note::

    The user needs to re-request file transfers after six months if they are still needed.

