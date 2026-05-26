# Pepper - ParticlE Physics ProcEssoR

[![gitlab repository](https://img.shields.io/badge/gitlab-repo-orange?logo=gitlab)](https://gitlab.cern.ch/cms-analysis/general/pepper)
[![pipeline status](https://gitlab.cern.ch/cms-analysis/general/pepper/pepper/badges/master/pipeline.svg)](https://gitlab.cern.ch/cms-analysis/general/pepper/-/pipelines)
[![documentation](https://img.shields.io/badge/docs-online-blue)](https://cms-pepper.docs.cern.ch/)
<!-- [![license](https://img.shields.io/badge/license-MIT-green)](LICENSE) -->

Pepper is a easy-to-use multi-purpose framework for analysing CMS NanoAOD datasets, built
on [coffea](https://coffea-hep.readthedocs.io/) and
[Awkward Array](https://awkward-array.org/). Originally developed for BSM searches with top quarks in dilepton final states, it is now used for both
searches and measurements across many final states and physics cases. For reference this repository comes with 
the necessary tools for a <img src="https://latex.codecogs.com/gif.latex?\mathrm{t\bar{t}}\rightarrow\mathrm{b\bar{b}}\mathrm{ll\nu\nu}" /> analysis.

📖 **Full documentation, tutorials, and the configuration reference are
available at [cms-pepper.docs.cern.ch](https://cms-pepper.docs.cern.ch/).**

## Installation

It is recommended to use an isolated virtual environment with Pepper. An example environment setup for NAF or LXPLUS 
is included in the repository under `example/environment.sh`. To install Pepper as an editable package and in 
a virtual environment, use the `install.sh` script:

```bash
git clone <repository url> pepper
cd pepper
chmod +x install.sh
./install.sh
# Activate your newly created environment
source environment.sh
```
Everytime you restart the shell, you need to activate the environment again. This just means running the `source environment.sh` command again. You can test your installation by running the included example script:

```bash
cd example
python3 -m pepper.runproc example_processor.py test_config.json
```

For full installation instructions, see the
[installation guide](https://cms-pepper.docs.cern.ch/install.html).

## Quick Start

```bash
cd example
python -m pepper.runproc example_processor.py example_config.json --debug
```

See the [Getting Started](https://cms-pepper.docs.cern.ch/user_guide/getting_started.html)
page for a walkthrough of the example processor and its output.

## Getting Help

Questions and discussion in the
[Pepper Mattermost channel](https://mattermost.web.cern.ch/cms-exp/channels/pepper-users).