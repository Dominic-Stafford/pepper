Examples
========

This section collects worked examples that build on top of the
:doc:`/user_guide/index`. Each example focuses on a specific task and walks
through it end to end with concrete configuration files and code, rather than
describing the available options exhaustively.

If you are looking for the complete list of available configuration options,
see the :ref:`configuration-reference`. If you are new to Pepper, start with
:doc:`/user_guide/getting_started`.

.. toctree::
   :hidden:

   Config Inheritance <config_inheritance>
   Extending the Config Class <extending_config>
   Using the Bundled Scripts <bundled_scripts>
   Applying Scale Factors <scale_factors>
   CMS Analyses Built with Pepper <cms_analyses>

.. grid:: 2 3 3 4
   :gutter: 2

   .. grid-item-card::
      :link: config_inheritance
      :link-type: doc
      :text-align: center

      Config Inheritance
      ^^^^^^^^^^^^^^^^^^
      Reuse a base configuration across multiple analyses or data-taking eras using the ``import`` key.

   .. grid-item-card::
      :link: extending_config
      :link-type: doc
      :text-align: center
 
      Extending ``Config``
      ^^^^^^^^^^^^^^^^^^^^
      Subclass :class:`pepper.Config` to add custom fields, behaviours, or special variables that JSON cannot express.
 

   .. grid-item-card::
      :link: bundled_scripts
      :link-type: doc
      :text-align: center
 
      Using the Bundled Scripts
      ^^^^^^^^^^^^^^^^^^^^^^^^^
      Walk through the scripts shipped in :repo:`scripts/` in the order
      you would run them for a new analysis.

   .. grid-item-card::
      :link: scale_factors
      :link-type: doc
      :text-align: center
 
      Applying Scale Factors
      ^^^^^^^^^^^^^^^^^^^^^^
      Walk through the process of configuring and applying scale factors to your analysis.

   .. grid-item-card::
      :link: cms_analyses
      :link-type: doc
      :text-align: center
 
      CMS Analyses Built with Pepper
      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      A gallery of CMS analyses that use Pepper, with links to their analysis-specific repositories.

.. Planned Examples
.. ----------------

.. The following examples are planned but not yet written. If you would like to
.. contribute one, or if you have a use case you would like to see documented,
.. please get in touch on :mattermost:`Mattermost <pepper-users>`.

.. - **CMS analyses built with Pepper.** A gallery of public CMS analyses that
..   use Pepper, with links to their analysis-specific repositories.