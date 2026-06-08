.. _cms-analyses:

CMS Analyses Built with Pepper
==============================

The list below collects public and internal CMS analyses that use Pepper as their analysis framework. 
It serves both as a reference for prospective users 
-- to see whether an analysis similar to yours has already been done in Pepper -- 
and as a pointer to worked, real-world code beyond the small examples shipped with the framework itself.

If your analysis uses Pepper and is not listed here, please reach out on :mattermost:`Mattermost <pepper-users>` so we can add it. 
Internal CADI entries are welcome alongside published results.

Legend
------
 
- **Status** -- the maturity of the analysis within CMS:
  
  - *Published* refers to a final paper or one submitted to a journal.
  - *Under review* refers to an analysis past the start of collaboration-wide review (CWR).
  - *Approved* refers to an analysis endorsed by its analysis working group (AWG).
  - *In progress* refers to an earlier-stage analysis with a CADI or AN entry but not yet AWG-endorsed.

- **Repo** -- where available, a link to the analysis-specific code repository. 
  An em-dash indicates that no public repository exists.


Analyses
--------
 
.. list-table::
   :header-rows: 1
   :widths: 18 52 15 10
 
   * - CADI / AN
     - Description
     - Status
     - Repo
   * - B2G (early)
     - Search for top-philic resonances in the boosted all-hadronic final state
     - In progress
     - --
   * - EXO-22-014
     - Search for dark matter produced in association with a single top
       quark or a top-quark pair using the full Run 2 data (dilepton channel)
     - Published
     - --
   * - EXO-24-020
     - Search for long-lived tau sleptons with a NN-based displaced tau tagger using the full Run 2 data
     - Published
     - --
   * - HIG-22-013
     - Search for :math:`A/H \to t\bar{t}` (dilepton channel)
     - Published
     - :repo-gitlab:`https://gitlab.cern.ch/ljeppe/pepper/-/tree/ahtt_run2ul_final`
   * - SMP / AN-2024/221
     - Differential cross section of electroweak :math:`Z+2` jets production with
       the full Run 2 data (dielectron, dimuon, and combined channels)
     - Approved
     - --
   * - TOP-22-008
     - Search for the SM :math:`tWZ` process in multi-lepton final states
     - Published
     - --
   * - TOP-22-012
     - Early Run 3 inclusive :math:`t\bar{t}` cross section (dilepton and
       lepton+jets channels)
     - Published
     - --
   * - TOP-23-002 / AN-2021/217
     - Inclusive and differential measurements of the :math:`t\bar{t}\gamma`
       cross section and the :math:`t\bar{t}\gamma/t\bar{t}` cross section
       ratio at :math:`\sqrt{s} = 13` TeV
     - Published
     - --
   * - TOP-23-004 / AN-2023/021
     - Simultaneous measurement of single and pair production of top quarks
       in association with a :math:`Z` boson
     - Published
     - --
   * - TOP-24-007
     - Observation of a structure at the :math:`t\bar{t}` threshold
     - Published
     - :repo-gitlab:`https://gitlab.cern.ch/ljeppe/pepper/-/tree/ahtt_run2ul_final`
   * - TOP-24-009
     - :math:`tWZ` cross section combining Run 2 and Run 3 data
     - Published
     - --
   * - TOP-25-004 / AN-2024/007
     - :math:`tW\gamma` measurement in Run 2
     - Approved
     - --
   * - TOP-25-012 / AN-2024/182
     - :math:`t\bar{t}` cross section at 13 TeV using events with at least
       one lepton
     - Approved
     - --
   * - TOP / AN-2025/088
     - :math:`t\bar{t}Z` cross section at 13 and 13.6 TeV with hadronic :math:`Z` decays
     - In progress
     - --
   * - TOP (early)
     - :math:`m_{\mathrm{top}}` using energy correlators (boosted
       lepton+jets channel)
     - In progress
     - --
   * - TOP (early)
     - :math:`t\bar{t}` differential cross section including the boosted
       region (dilepton channel)
     - In progress
     - --


Adding an analysis
------------------
 
To add an analysis to the gallery, open a merge request against the
Pepper documentation with a new row in the table above. The four fields
to provide are:
 
- **CADI / AN**: the CMS CADI identifier (e.g. ``TOP-23-001``). 
  For analyses without a CADI yet, the AN number (e.g. ``AN-2024/123``) is also acceptable, 
  optionally prefixed with the target PAG. For early-stage analyses without either, use the PAG name followed by ``(early)``.
- **Description**: a single sentence naming the process measured or searched for, 
  the final state, and the centre-of-mass energy where relevant. 
  Use :math:`\dots` for particle names.
- **Status**: ``Published``, ``Under review``, ``Approved``, or
  ``In progress``, as defined in the legend above.
- **Repo**: the public analysis repository if one exists, otherwise an em-dash. 
  Format as ``:repo-gitlab:`<URL>``` for CERN GitLab or ``:repo-github:`<URL>``` for GitHub.
 
Keeping descriptions to a single sentence is intentional: 
the gallery is a directory, not a summary. 
Readers interested in a specific analysis should follow the repo link.
