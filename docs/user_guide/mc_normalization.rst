.. _mc_normalization:

Monte Carlo Normalization
================

Monte Carlo (MC) samples are produced with arbitrary numbers of events. 
To compare them to data, simulated yields must be scaled so that the predicted 
yield in any region equals the cross section of the process times the integrated 
luminosity of the data. This page explains how Pepper performs that scaling, 
the two ways of providing the normalization, and when each is appropriate.

The configuration keys touched here -- ``luminosity``, ``crosssections``,
``mc_lumifactors`` -- are documented in the :ref:`configuration-reference`.


What is a lumifactor?
---------------------

For each MC dataset describing some process, Pepper computes a single scalar called the
*lumifactor*:

.. math::

   f = \frac{\mathcal{L} \cdot \sigma}{\sum w_{\mathrm{gen}}}

where :math:`\mathcal{L}` is the integrated luminosity of the data sample 
(the ``luminosity`` config key), 
:math:`\sigma` is the cross section of the MC process (from ``crosssections`` config key), 
and :math:`\sum w_{\mathrm{gen}}` is the sum of generator weights over
every event in the dataset *before* any selection.

At the end of the run, every histogram and cutflow entry for an MC dataset is 
multiplied by :math:`f` for that dataset. 
The result is a yield in units of "expected events at 
luminosity :math:`\mathcal{L}`," which is directly comparable to the data yield.

This is the only normalization Pepper applies on top of the per-event
weight chain (pileup, scale factors, etc.). 
The lumifactor is a constant per dataset, computed once and applied at the end.


Two ways to provide :math:`\sum w_{\mathrm{gen}}`
-------------------------------------------------

The cross-sections and luminosity are configuration inputs. 
The :math:`\sum w_{\mathrm{gen}}` denominator has to come from the MC files themselves 
-- there is no way around opening every file in a dataset and reading the sum of generator weights. 
Pepper supports two strategies for this:

**Runtime (posterior) computation (recommended)**: set ``"mc_lumifactors": false`` in the config. 
Pepper computes :math:`\sum w_{\mathrm{gen}}` on the fly, 
summing generator weights from every chunk as it is processed, 
and applies the lumifactor when the output is saved.

**Pre-computed**: run :repo:`scripts/compute_mc_lumifactors.py` once,
which opens every MC file, accumulates the sum of generator weights per dataset, 
and writes the result to a JSON file. Point the ``mc_lumifactors`` config key at that file. 
The processor reads it at startup and uses the stored values.

The two strategies produce identical results when every event in all
files are successfully processed. 


Why runtime computation is the recommended default
---------------------------------------------------

Production runs commonly involve thousands of input files distributed across hundreds of Condor jobs. 
A small fraction of those jobs will fail for transient reasons 
-- node crashes, storage timeouts, corrupt remote replicas -- even after retries. 
Pepper supports two policies for handling unrecoverable MC failures
(:ref:`htcondor`): exit the whole run (``exit_on_failed_jobs: "all"``),
or tolerate MC failures and continue (``exit_on_failed_jobs: "data"``). 
The second policy is attractive to avoid having to re-run the entire production because of one bad MC file.

With pre-computed lumifactors, ``"data"`` mode is dangerous. 
The denominator :math:`\sum w_{\mathrm{gen}}` in the lumifactor JSON was computed over the *full* MC sample, 
but only a subset of that sample ended up in the histograms. 
The output is normalized to too manyevents and the predicted yield is biased low. 
Pepper will warn when this combination is configured, but the result is still wrong.

With runtime computation, :math:`\sum w_{\mathrm{gen}}` is summed from
exactly the same events that fill the histograms. 
If a chunk fails, both its contribution to the histograms *and* its contribution to the
denominator are missing, and the ratio is unbiased. 
The normalization remains statistically correct (with larger uncertainties from the
smaller effective sample) regardless of which jobs failed.

The recommended configuration for a typical production run is therefore:

.. code-block:: json

   {
       "mc_lumifactors": false,
       "exit_on_failed_jobs": "data"
   }

This combination tolerates MC job failures without producing biased normalizations, 
and is the simplest workflow: no pre-processing script to run, no JSON file to keep in sync with the dataset list.


When to use ``compute_mc_lumifactors.py`` instead
-------------------------------------------------

Runtime computation is not always sufficient or desireable. The following cases require pre-computed lumifactors:

**Skimmed or pre-filtered NanoAOD.** Some NanoAOD productions are filtered (i.e. skimmed) before being written 
-- only events passing a loose preselection appear in the files. 
The events in the file are then not a complete sample of the underlying generator events, 
and the :math:`\sum w_{\mathrm{gen}}` summed at runtime is too small. 
The correct denominator is the sum over the original unfiltered events,
which only the upstream production knows. 
``compute_mc_lumifactors.py`` reads this from the ``Runs`` tree of NanoAOD (which carries the unfiltered sum), giving the right answer.

**Partial-dataset runs with absolute-yield comparisons.** If you deliberately process only part of an MC dataset 
-- one era, one file subset, a quick test -- 
runtime computation will normalize to the events you actually ran on. 
Pre-computed lumifactors normalize to the full nominal sample. 
The latter is what you want when comparing absolute yields across runs that processed different subsets, 
or when quoting yields in a note. This is also useful if you have a single MC sample on disk and
you want to check if and how it enters the event-selection.

**PDF and scale-weight normalization without splitting.** 
It is common in many analyses to normalize the PDF weights before any cuts, 
so that these only account for the shape and selection acceptance effects of the uncertainties 
for missing higher orders in the PDF calculation, respectively, 
with the normalisation taken as part of the cross section uncertainty for this sample 
(using the ``"normalize_pdf_uncs": true`` option). 
Pepper can also compute the total pdf variation, rather than returning the individual pdf sources, 
with the ``"split_pdf_uncs": false`` option (note this is not recommended for many analyses, as for any samples with negative weights, 
e.g. from ``aMC@NLO``, this will over-estimate the total variation). 
If using these two options in conjuction, the total ``sumweights`` for these samples must be precomputed using ``compute_mc_lumifactors.py`` with the ``-p`` option.

**Per-event output.** When the processor writes per-event columns to disk (``columns_to_save``), 
the runtime lumifactor is not applied to those files 
-- Pepper has no way to retroactively scale already-written per-event records. 
The per-dataset :math:`\sum w_{\mathrm{gen}}` is saved alongside the output as ``gen_sumws.json`` 
so downstream code can normalize the columns itself, 
but if you want columns to be self-normalizing, use pre-computed lumifactors. 
Pepper logs a warning in this case. In short, for per-event output you need to ensure normalization either
by using ``compute_mc_lumifactors.py`` or manually setting up a script using ``gen_sumws.json``.

In all cases, run ``compute_mc_lumifactors.py`` once per dataset list, commit the resulting JSON alongside the config, 
and use ``exit_on_failed_jobs: "all"`` so that any MC failure stops the run and forces investigation rather than producing a biased result.


See also
--------

- :ref:`configuration-reference` -- detailed documentation of
  ``luminosity``, ``crosssections``, ``mc_lumifactors``,
  ``norm_genweights``, ``split_pdf_uncs``, and ``normalize_pdf_uncs``.
- :ref:`scale-factors-example` -- where MC normalization fits into
  the broader correction workflow.
- :ref:`htcondor` -- the ``exit_on_failed_jobs`` policy and the
  failure modes it interacts with.
- :ref:`bundled-scripts-example` -- when ``compute_mc_lumifactors.py``
  fits in a complete analysis run.