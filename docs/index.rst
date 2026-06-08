.. Pepper documentation master file, created by
   sphinx-quickstart on Wed Dec  3 10:40:35 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Pepper's documentation!
==================================


A Python framework for HEP analysis on CMS NanoAOD datasets, built on `coffea <https://coffea-hep.readthedocs.io/>`_ 
and `Awkward Array <https://awkward-array.org/>`_.

Pepper handles the routine machinery of a CMS analysis -- object selection, scale factors, 
systematic uncertainties, histogramming, batch submission 
-- so that analysts can focus on the parts that are specific to their measurement or search.


.. grid:: 2 2 4 4
   :gutter: 2

   .. grid-item-card::
      :link: user_guide/getting_started
      :link-type: doc
      :text-align: center

      Get Started
      ^^^^^^^^^^^
      Install Pepper and run your first analysis.

   .. grid-item-card::
      :link: user_guide/conceptual_overview
      :link-type: doc
      :text-align: center

      Concepts
      ^^^^^^^^
      Understand how Pepper extends coffea.

   .. grid-item-card::
      :link: examples/index
      :link-type: doc
      :text-align: center

      Examples
      ^^^^^^^^
      Worked examples covering common tasks.

   .. grid-item-card::
      :link: api/index
      :link-type: doc
      :text-align: center

      API Reference
      ^^^^^^^^^^^^^
      Full module-level reference.

This documentation website is currently work-in-progress. If you encounter any issues or have suggestions for improvement, 
please feel free to reach out to Mads on `e-mail <mailto:mads.baattrup@desy.de?subject=Feedback\ on\ Pepper\ Documentation>`__ or Mattermost.


.. admonition:: CMS-endorsed framework
   :class: note

   Pepper is listed among the 
   `supported frameworks of the CMS collaboration <https://cms-analysis.docs.cern.ch/guidelines/frameworks/frameworks/>`_
   for use in physics analyses.

Citing Pepper
-------------

If you use Pepper in your research, please consider citing it in your publications. Here is a BibTeX entry you can use:

.. code-block:: bibtex

   @software{pepper2025,
      author       = {Pepper Maintainers},
      title        = {Pepper - ParticlE Physics ProcEssoR},
      year         = 2025,
      url          = {https://gitlab.cern.ch/cms-analysis/general/pepper},
   }


.. toctree::
   :hidden:
   :maxdepth: 2

   install
   user_guide/index
   examples/index
   api/index
   support

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
