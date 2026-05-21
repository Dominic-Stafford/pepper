Understanding the ``Config`` Object
===================================
The ``Config`` object is not merely a dictionary of settings corresponding to the analysis configuration JSON. 
It is a specialized class that encapsulates configuration data and provides methods for accessing, validating these settings.
The flexible design also allows users to extend the ``Config`` fields beyond those defined in the default configuration schema
and beyond what is directly representable in JSON. This means users can add custom fields or complex data structures to the 
configuration without being constrained by the limitations of JSON serialization. Similar to how one would extend processor 
classes to implement custom analysis logic, users can subclass the ``Config`` class to introduce new configuration options
or modify existing ones to better suit their analysis needs. It follows the same hierarchy with a base ``Config`` class and
and a ``ConfigBasicPhysics`` subclass that implements standard HEP analysis configuration options.

--------------------------
The ``Config`` Base Class
--------------------------

The ``Config`` class implements the standard Python ``MutableMapping`` to provide dictionary-like access to configuration settings.
The standard python-dict behaviour is extended by the getter which allows for accessing single fields or lists of fields:

.. code-block:: python3

    >>> config["lumi"]  # Access single field like Python dicts
    59.74  # some value in fb^-1
    >>> config[["lumi", "year"]]  # Access multiple fields
    [59.74, "2018"]

In the constructor method of the ``Config`` class, the following groups are defined:

- Required keys are defined. If these keys are not present, the method ``check_required_args`` will raise a ``ConfigError``.
- The special variables are defined. The standard special-variables (usually prefixed with ``$``) are explained in more detail :ref:`here <special-variables>`. 
  In short, they are substituted with another key at runtime which greatly reduces overhead when subclassing configs. 
  You can read more about config inheritance :ref:`here <config-inheritance-example>`.
- Finally, special behaviour for config keys are defined. 
  These are functions applied to the JSON-serializable string when a specific key is queried. 
  For instance, we could define a behaviour that takes a path to a file and reads the contents 
  into memory instead of having to manually do that every time the key is queried.

The list of standard required keys, behaviours and special variables can be seen in the reference describing the :class:`pepper.config.Config` class.

.. hint:: 

    The ``Config`` class automatically caches the result of configuration lookups such that subsequent accesses to a key will return the cached value instead of re-evaluating the configuration. 
    This is particularly useful for expensive or time-consuming configuration behaviours that would otherwise need to be re-evaluated on each access.

-----------------------------------
The ``ConfigBasicPhysics`` Subclass
-----------------------------------

The ``ConfigBasicPhysics`` subclass extends the base ``Config`` class by adding HEP analysis-specific configuration options.
It introduces additional required keys, behaviours, and special variables relevant to typical HEP analyses. For example,
it adds behaviours that automatically builds Coffea Corrector objects when querying correction-related keys. This makes
it incrediblty easy to apply standard corrections without having to manually instantiate these objects in the analysis code.
The list of additional required keys, behaviours and special variables can be seen in the reference describing the :class:`pepper.config_basic.ConfigBasicPhysics` class.

.. note::

    Users can further extend both the ``Config`` and ``ConfigBasicPhysics`` classes by subclassing them to introduce custom configuration options or modify existing ones to better suit their analysis needs.\
    This is particularly useful for analyses with specialized requirements not covered by the default configuration schema.


-------------------
Config Inheritance
-------------------

The config ``json`` files support inheritance, which allows users to create new configuration files that inherit settings from existing ones. 
This is achieved through the use of the ``import`` special variable, which specifies the path to the parent configuration file. 
When a configuration file is loaded, any settings defined in the child file will override those in the parent file, while any settings not defined in the child will be inherited from the parent. 
This feature promotes reusability and modularity in configuration management, allowing users to create a base configuration with common settings and then extend it for specific analyses without having to duplicate all the settings.

This can be particularly useful when working on a large analysis with multiple data-taking eras for instance. 
You can create a base configuration file with common settings for all eras and then create separate configuration files for each era that inherit from the base configuration and only specify the settings that differ between eras.
Please refer to the example page :ref:`config-inheritance-example` for an example of how to use config inheritance in practice.
