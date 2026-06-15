.. _migrating:

Migrating from Awkward 1 / Coffea 0.7
======================================

Pepper has been updated to use **Awkward Array 2** and **Coffea 2026**.
Most user-facing processor code continues to work without changes, but a handful of
patterns that were common in the Awkward 1 era have changed or been removed.
This page describes the most important differences.

Environment and packages
------------------------

Pepper now requires up-to-date versions of ``awkward``, ``coffea`` and their dependencies.
It is likely easiest to reinstall the environment from scratch, following the instruction in :ref:`installation`.
In case you build the environment on top of LCG, use a newer LCG version that already contains
Awkward Array 2, e.g. ``LCG_109a``.


Lorentz Vector API
------------------

The vector library used for Lorentz vectors has changed.
The attribute and method names that are most likely to appear in user code are listed below.

**Three-momentum magnitude**

+------------------------------------+------------------------------------+
| Awkward 1                          | Awkward 2                          |
+====================================+====================================+
| ``vec.rho``  (3-vector magnitude)  | ``vec.p``                          |
+------------------------------------+------------------------------------+
| ``vec.rho2`` (squared magnitude)   | ``vec.p2``                         |
+------------------------------------+------------------------------------+

Example::

    # Old
    chel = lep_a.dot(lep_b) / lep_a.rho / lep_b.rho
    energy = np.sqrt(mt**2 + top.rho2)

    # New
    chel = lep_a.pvec.dot(lep_b.pvec) / lep_a.p / lep_b.p
    energy = np.sqrt(mt**2 + top.p2)

**Dot product of 3-momenta**

In Awkward 1, calling ``.dot()`` directly on a 4-vector computed the 3-vector (spatial)
dot product. In Awkward 2 the behaviour is different.
Extract the 3-vector first with ``.pvec`` before calling ``.dot()``:

.. code-block:: python

    # Old
    cos_angle = a.dot(b) / a.rho / b.rho

    # New
    cos_angle = a.pvec.dot(b.pvec) / a.p / b.p

Boolean Indexing of Ragged Arrays
----------------------------------

In Awkward 1, boolean indexing ``array[mask]`` would drop ``False`` entries even when the
mask had a ``keepdims``-style shape (e.g. shape ``(n, 1)``). In Awkward 2 the semantics
are stricter: to apply an optional mask and produce ``None`` for failing entries, use
``ak.mask``; then call ``ak.drop_none`` if you want to remove the ``None`` entries
rather than keep them as optional values.

.. code-block:: python

    # Old
    result = ak.zip({...})[has_solution] / normalisation

    # New
    result = ak.drop_none(ak.mask(ak.zip({...}), has_solution)) / normalisation

Removal of ``lazy=True`` in ``set_column``
------------------------------------------

The ``lazy`` keyword argument of ``selector.set_column`` (and the underlying
``pepper.misc.VirtualArrayCopier`` helper) has been removed.

All columns are now evaluated eagerly when ``set_column`` is called.  Simply remove the
``lazy=True`` keyword from any existing calls:

.. code-block:: python

    # Old
    selector.set_column("dilep_pt", self.dilep_pt, lazy=True)

    # New
    selector.set_column("dilep_pt", self.dilep_pt)

``ak.packed`` renamed to ``ak.to_packed``
------------------------------------------

The function ``ak.packed`` was renamed to ``ak.to_packed`` in Awkward 2.
Update any calls in custom code:

.. code-block:: python

    # Old
    array = ak.packed(array)

    # New
    array = ak.to_packed(array)

State-file compatibility
------------------------

The internal executor state file format was bumped from version 3 to version 4.
Old ``pepper_metadata.coffea`` and ``pepper_state.coffea`` files created with Coffea 0.7 / Pepper on Awkward 1
are **not** compatible with the new version.  Delete any existing state files before
running with the updated code, or start a fresh output directory.
