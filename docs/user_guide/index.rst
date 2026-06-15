User Guide
==========

.. toctree::
   :hidden:

   Get Started <getting_started>
   Concepts <conceptual_overview>
   Configuration <config>
   The Config Object <config_object>
   Pepper Scripts <scripts>
   Run Your Analysis <running_a_processor>
   Accessing Remote Data <remote_data>
   Scale Factors <scale_factors>
   Histograms <histograms>
   HTCondor <htcondor>
   Memory Optimization <optimization>
   Migrating from Awkward 1 <migrating>
   MC Normalization <mc_normalization>

Welcome to the Pepper user guides. Here, you can find various introductions to some of the features available in Pepper. 
In case you have suggestions for what could also be included here, please reach out on :mattermost:`Mattermost <pepper-users>`!

.. grid:: 2 3 3 4
   :gutter: 2

   .. grid-item-card::
      :link: getting_started
      :link-type: doc
      :text-align: center

      Get Started
      ^^^^^^^^^^^
      Get started here!


   .. grid-item-card::
      :link: conceptual_overview
      :link-type: doc
      :text-align: center

      Concepts
      ^^^^^^^^
      Understand the core concepts.

   .. grid-item-card::
      :link: config
      :link-type: doc
      :text-align: center

      Configuration
      ^^^^^^^^^^^^^
      Understand the options when configuring an analysis.

   .. grid-item-card::
      :link: config_object
      :link-type: doc
      :text-align: center

      The Config Object
      ^^^^^^^^^^^^^^^^^^
      Understand the flexibility of the ``Config`` object and how to extend it.

.. grid:: 2 3 3 4
   :gutter: 2

   .. grid-item-card::
      :link: scripts
      :link-type: doc
      :text-align: center

      Pepper Scripts
      ^^^^^^^^^^^^^^^
      Perform common tasks related to running your analysis.

   .. grid-item-card::
      :link: running_a_processor
      :link-type: doc
      :text-align: center

      Run Your Analysis
      ^^^^^^^^^^^^^^^^^^^^^
      How to run your analysis - locally and on the cluster.

   .. grid-item-card::
      :link: remote_data
      :link-type: doc
      :text-align: center

      Remote Data
      ^^^^^^^^^^^^^^
      Learn how to access and process data from remote sources.

   .. grid-item-card::
      :link: scale_factors
      :link-type: doc
      :text-align: center

      Scale Factors
      ^^^^^^^^^^^^^
      Learn about scale factors and how to apply them in your analysis.

   .. grid-item-card::
      :link: optimization
      :link-type: doc
      :text-align: center

      Optimization
      ^^^^^^^^^^^^
      Optimize memory consumption of the analysis.

   .. grid-item-card::
      :link: migrating
      :link-type: doc
      :text-align: center

      Migration Guide
      ^^^^^^^^^^^^^^^
      Migrate existing processor code from Awkward 1 / Coffea 0.7.


.. grid:: 2 3 3 4
   :gutter: 2

   .. grid-item-card::
      :link: histograms
      :link-type: doc
      :text-align: center

      Histogramming
      ^^^^^^^^^^^^^
      How to configure and produce histograms with Pepper.

   .. grid-item-card::
      :link: htcondor
      :link-type: doc
      :text-align: center

      HTCondor
      ^^^^^^^^^^^^^
      How to scale up your analysis using HTCondor.

   .. grid-item-card::
      :link: mc_normalization
      :link-type: doc
      :text-align: center

      MC Normalization
      ^^^^^^^^^^^^^^^^
      How to configure MC normalization in Pepper.

Planned User Guide
---------------------

The following examples are planned but not yet written. If you would like to
contribute one, or if you have a use case you would like to see documented,
please get in touch on :mattermost:`Mattermost <pepper-users>`.

- **Systematic variations.** A walkthrough of how to configure and apply systematic variations in Pepper, including the
  ``compute_systematics`` config key and how it interacts with the processor code.