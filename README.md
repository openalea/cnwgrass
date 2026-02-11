[![License](https://img.shields.io/badge/license-CeCILL--C-blue )](https://img.shields.io/badge/license-CeCILL--C-blue )

[![Documentation Status](https://readthedocs.org/projects/cn-wheat/badge/?version=latest)](https://cn-wheat.readthedocs.io/en/latest/?badge=latest)

# WheatFspm
WheatFspm is a Functional Structural Plant Model (FSPM) of wheat which fully integrates shoot morphogenesis and the metabolism of carbon (C) and nitrogen (N) at organ scale within a 3D representation of plant architecture. Plants are described as a collection of tillers, each consisting in individual shoot organs (lamina, sheath, internode, peduncle, chaff), a single root compartment, the grains, and a phloem.

WheatFspm simulates:
* Organ photosynthesis, temperature and transpiration from light distribution within the 3D canopy which is provided by Caribu and Adel-Wheat models ((Chelle and Andrieu, 1998) ; (Fournier *et al.*, 2003), respectively)
* Leaf and internode elongation
* Leaf, internode and root growth in mass
* N acquisition, synthesis and allocation of C and N metabolites at organ level and among tiller organs
* Senescence of shoot organs and roots

Model inputs are the climatic conditions (temperature, light, humidity, CO2, wind, soil NO<sub>3</sub><sup>-</sup>) and initial dimensions, mass and metabolic composition of individual organs.

![alt text](_media/Vegetative_stages_topview.gif "Growing canopy")

# Description
WheatFspm consists in a set of sub-models (named submodules in git) which share inputs/outputs through an MTG object:

![alt text](_media/Modular_structure.png "WheatFSPM workflow") 
*Adapted from Gauthier et al. (2020)*

* *Farquhar-Wheat*: Farquhar-based model of photosynthesis, stomatal conductance, organ temperature and transpiration.
* *Elong-Wheat*: regulation of leaf and internode elongation by C and N metabolites, temperature and coordination rules.
* *Growth-Wheat*: growth in biomass of leaves, internodes and roots ; related consumption in C and N metabolites.
* *CN-Wheat*: synthesis and degradation of C and N metabolites at organ level and allocation between tillers' organs. See doc at https://cn-wheat.readthedocs.io/ 
* *Respi-Wheat*: respiratory-costs related to the main biological processes.
* *Senesc-Wheat*: organ senescence and consequences in organ biomass, green area and remoblisation of C and N metabolites.
* *Fspm-Wheat*: is the submodule containing the interfaces (facades) for reading/updating information between each sub-model and the MTG. Also includes the scripts to be run for using all sub-models.


# Table of Contents
- [Table of Contents](#table-of-contents)
- [Installation](#installation)
  * [Prerequisites](#prerequisites)
  * [Installing](#installing)
    + [Users](#users)
    + [Developers](#developers)
- [Usage](#usage)
  * [NEMA](#nema)
  * [Papier_FSPMA2016](#papier-fspma2016)
  * [Vegetative stages](#vegetative-stages)
  * [Scenarios_monoculms](#scenarios-monoculms)
- [Credits](#credits)
  * [Authors](#authors)
  * [Contributors](#contributors)
  * [Funding](#funding)
- [License](#license)


# Installation

*WheatFspm* has been tested on Windows 10 64 bit and Linux Fedora 24 64 bit. *WheatFspm* is now developed in Python 3.

## Prerequisites
*WheatFspm* has the following dependencies (see documentation in the links provided, instructions for their installation are given in [Installing](#installing)):
* To run the model: 
    * [Python](http://www.python.org) >= 3.7
    * [NumPy](http://www.numpy.org/) >= 1.7.2
    * [SciPy](http://www.scipy.org/) >= 0.12.1
    * [Pandas](http://pandas.pydata.org/) >= 0.14.0 
    * [Openalea.MTG](https://github.com/openalea/mtg)
    * [Openalea.Plantgl](https://github.com/openalea/plantgl)
    * [Openalea.Lpy](https://github.com/openalea/lpy)
    * [Alinea.Caribu](https://github.com/openalea-incubator/caribu) 
    * [Alinea.Adel](https://github.com/openalea-incubator/adel)
    * [Alinea.Astk](https://github.com/openalea-incubator/astk)
    
* To run the tools: [Matplotlib](http://matplotlib.org/) >= 1.3.1 
* To build the documentation: [Sphinx](http://sphinx-doc.org/) >= 1.1.3
* To run the tests: [Nose](http://nose.readthedocs.org/) >= 1.3.0 

## Installing
*WheatFspm* has to be installed in a conda environment containing all dependencies.

* Install Miniconda 2 or 3 for Python 3.7: https://docs.conda.io/en/latest/miniconda.html
* Create a new environment in an Anaconda prompt:
   `conda create -n WheatFspm python=3.7 openalea.mtg openalea.plantgl openalea.lpy alinea.caribu alinea.astk coverage nose sphinx statsmodels -c conda-forge -c fredboudon`
* Activate the conda environment:
    `conda activate WheatFspm`
* Install Adel-Wheat:
    Download and extract the following archive https://github.com/rbarillot/adel/archive/python3.zip (temporary branch for python 3 compatibility)
    cd to your local reposoitory of adel and install it in your conda environment `python setup.py develop`
* Install FspmWheat:
    'todo'
   
### Users
* Download WheatFspm package from https://github.com/openalea-incubator/WheatFspm/archive/master.zip
* Extract archive
* Open a command line interpreter in your local copy of *WheatFspm*,
* Run command: `python setup.py install --user` 

### Developers

This  package contains Git [submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) ([in French](https://git-scm.com/book/fr/v2/Utilitaires-Git-Sous-modules)).

Development is done in the different submodules.

#### Cloning

To clone the project, please use:

    git clone --recurse-submodules https://github.com/openalea-incubator/WheatFspm

#### Updating submodules

If you want to update all submodules:

    git submodule update --remote

Otherelse, update each one to a specific version, branch or tag, do:

    cd mypackage
    git fetch
    git merge origin/master

#### Install packages in develop mode

* Open a command line interpreter into each copy of submodels
* Run command: `python setup.py develop --user`

# Usage

To date, *WheatFspm* has been used in four main contexts described below. 

The scripts to run *WheatFSPM* are located in:
* `WheatFspm\fspm-wheat\example\NEMA`
* `WheatFspm\fspm-wheat\example\Papier_FSPMA2016`
* `WheatFspm\fspm-wheat\example\Vegetative_stages`
* `WheatFspm\fspm-wheat\example\Scenarii_monoculms`

## NEMA
This example deals with the post-flowering stages of wheat developement under 3 nitrogen fertilisation regies (H0, H3 and H15). The main processes described are leaf senescence, C and N remobilisation, grain filling). During that stages, all vegetative organs have completed their growth. 
This work led to the research articles [Barillot *et al.* (2016a)](https://doi.org/10.1093/aob/mcw143) and [Barillot *et al.* (2016b)](https://doi.org/10.1093/aob/mcw144).
This example has been maintained in the current version ; results of above papers were generated using the tag *release-1.0* of *CN-Wheat*. 

To run the example:
* Open a command line interpreter in `WheatFspm\fspm-wheat\example\NEMA`
* Run each script called *main.py* located in *NEMA_H0*, *NEMA_H3*, *NEMA_H15*: `python main.py`

## Papier_FSPMA2016
This example deals with the effects of leaf inclination, radiations regimes, plant density and sowing patterns on plant metabolism and grain filling during the post-flowering stages.   
This work led to the research article [Barillot *et al.* (2019)](https://doi.org/10.1093/aob/mcy208).
The scripts have not been maintained in the current version but are available using tags *paper_FSPMA16* or each submodule listed in [Description](#description).

To run the example:
* Open a command line interpreter in `WheatFspm\fspm-wheat\example\Papier_FSPMA2016`
* Run each script called *main.py* located in the different sub-directories: `python main.py`

## Vegetative stages
This example deals with the early vegetative stages of wheat development. It mainly covers the processes of leaf, internode and roots growth.
Tillering is simplified: tiller emergence is a model input while tiller metabolism and growth is approximated from that  of the main stem.
This work led to the research article [Gauthier *et al.* (2020)](https://doi.org/10.1093/jxb/eraa276). Results were obtained from the tag [paper_JXBot_2020](https://github.com/openalea-incubator/WheatFspm/releases/tag/paper_JXBot_2020).
To run the model used for the paper, please download the code archives at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3759714.svg)](https://doi.org/10.5281/zenodo.3759714)
 
To run the example:
* Open a command line interpreter in `WheatFspm\fspm-wheat\example\Vegetative_stages`
* Run script *main.py*: `python main.py`

## Scenarios_monoculms
This example deals with the plasticity of leaf growth during the vegetative stages of wheat development. The growth of wheat monoculms was simulated for highly contrasting conditions of soil nitrogen concentration, incident light and planting density.
This work led to the research article [Gauthier *et al.* (2021)](https://doi.org/10.1093/insilicoplants/diab034). Results were obtained from the following sources: [cn-wheat](https://github.com/openalea-incubator/cn-wheat/releases/tag/v2.0), [elong-wheat](https://github.com/openalea-incubator/elong-wheat/releases/tag/v3.0), [growth-wheat](https://github.com/openalea-incubator/growth-wheat/releases/tag/v2.0), [farquhar-wheat](https://github.com/openalea-incubator/farquhar-wheat/releases/tag/v2.0), [senesc-wheat](https://github.com/openalea-incubator/senesc-wheat/releases/tag/2.0), [respi-wheat](https://github.com/openalea-incubator/respi-wheat/releases/tag/v2.0) and [fspmwheat-wheat](https://github.com/openalea-incubator/fspm-wheat/releases/tag/release-2.0).

To run the model used for the paper, please download the code archives at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5503312.svg)](https://doi.org/10.5281/zenodo.5503312)
 
To run the example:
* Open a command line interpreter in `WheatFspm\fspm-wheat\example\Scenarios_monoculms`
* For a single scenraio, run the script *main.py*: `python main.py`
* The whole set of scenarios was run in the high-performance computing center [MESO@LR](https://meso-lr.umontpellier.fr/) (Université de Montpellier, France) 

# Credits
## Authors
* Romain BARILLOT - model designing, development and validation - [rbarillot](https://github.com/rbarillot)
* Marion GAUTHIER - model designing, development and validation - [mngauthier](https://github.com/mngauthier)
* Camille CHAMBON - software designing, development, deployment and optimization - [cachambon](https://github.com/cachambon)
* Bruno ANDRIEU - model designing and validation, scientific project management - [bandrieu](https://orcid.org/0000-0002-7933-9490)

## Contributors
* Christian FOURNIER - [christian34](https://github.com/christian34)
* Christophe PRADAL - [pradal](https://github.com/pradal)

## Funding
* [INRAE](https://www.inrae.fr/): salaries of permanent staff 
* French Research National Agency: projects [Breedwheat](https://breedwheat.fr/) (ANR-10-BTBR-03) and [Wheatamix](https://www6.inrae.fr/wheatamix/) (ANR-13-AGRO0008): postdoctoral research of R.Barillot
* [itk](https://www.itk.fr/en/) company and [ANRT](http://www.anrt.asso.fr/fr): funded the [Cifre](http://www.anrt.asso.fr/fr/cifre-7843) PhD thesis of M.Gauthier

# License
This project is licensed under the CeCILL-C License - see file [LICENSE](LICENSE) for details
 
