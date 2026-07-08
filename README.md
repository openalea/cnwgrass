# CNW-Grass
[![License](https://img.shields.io/badge/license-CeCILL--C-blue )](https://img.shields.io/badge/license-CeCILL--C-blue )

[![OpenAlea CI](https://github.com/openalea/cnwgrass/actions/workflows/openalea_ci.yml/badge.svg?branch=master)](https://github.com/openalea/cnwgrass/actions/workflows/openalea_ci.yml)
[![Documentation Status](https://readthedocs.org/projects/cnwgrass/badge/?version=latest)](https://cnwgrass.readthedocs.io/en/latest/?badge=latest)

[![Platform](https://anaconda.org/openalea3/openalea.cnw-grass/badges/version.svg)](https://anaconda.org/openalea3/openalea.cnw-grass)
[![Platform](https://anaconda.org/openalea3/openalea.cnw-grass/badges/platforms.svg)](https://anaconda.org/openalea3/openalea.cnw-grass)
[![Python Version](https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org/downloads/)

[![Anaconda-Server Badge](https://anaconda.org/openalea3/openalea.cnw-grass/badges/latest_release_relative_date.svg)](https://anaconda.org/openalea3/openalea.cnw-grass)

[![Downloads](https://anaconda.org/openalea3/openalea.cnw-grass/badges/downloads.svg)](https://anaconda.org/openalea3/openalea.cnw-grass)

## About

CNW-Grass is a Functional Structural Plant Model (FSPM) of grasses which fully integrates shoot morphogenesis and the 
metabolism of carbon (C) and nitrogen (N) at organ scale within a 3D representation of plant architecture. 
Plants are described as a collection of tillers, each consisting in individual shoot organs (lamina, sheath, internode, peduncle, chaff), 
a single root compartment, the grains, and a phloem.
CNW-Grass also includes a hydraulic model allowing to compute water flow in the plant and the co-regulation of 
leaf growth by metabolic and hydraulic processes. In this case, the plants also include a xylem compartment.

CNW-Grass simulates:
* Organ photosynthesis, temperature and transpiration from light distribution within the 3D canopy.
* Leaf and internode elongation.
* Leaf, internode and root growth in mass.
* N acquisition, synthesis and allocation of C and N metabolites at organ level and among tiller organs.
* Senescence of shoot organs and roots.
* Water fluxes and water potentials.

Model inputs are the pedoclimatic conditions (temperature, light, humidity, CO<sub>2</sub>, wind, 
soil NO<sub>3</sub><sup>-</sup>, Soil Relative Water Content) and initial dimensions, mass and metabolic composition of individual organs.

![Growing canopy](https://github.com/openalea/cnwgrass/blob/master/doc/_static/Vegetative_stages_topview.gif?raw=true "Growing canopy")

# Description
CNW-Grass consists in a set of sub-models (named submodules in git) which share inputs/outputs through an MTG object:

![CNW-Grass workflow](https://github.com/openalea/cnwgrass/blob/master/doc/_static/Modular_structure.png?raw=true "CNW-Grass workflow") 
*Adapted from Gauthier et al. (2020)*

* *Gas-Exchange*: Farquhar-based model of photosynthesis, stomatal conductance, organ temperature and transpiration.
* *Morphogenesis*: regulation of leaf and internode elongation by C and N metabolites, temperature and coordination rules.
* *Growth*: growth in biomass of leaves, internodes and roots ; related consumption in C and N metabolites.
* *CN-Metabolism*: synthesis and degradation of C and N metabolites at organ level and allocation between tillers' organs. 
* *Hydraulics*: water fluxes, organ water potential and co-regulation of leaf growth with Morphogenesis.
* *Respiration*: respiratory-costs related to the main biological processes.
* *Senescence*: organ senescence and consequences in organ biomass, green area and remobilisation of C and N metabolites.
* *Integration*: is the submodule containing the interfaces (facades) for reading/updating information between each sub-model and the MTG. Also includes the scripts to be run for using all sub-models.

Full documentation of each submodule is available at https://cnwgrass.readthedocs.io/

# Table of Contents
- [Table of Contents](#table-of-contents)
- [Installation](#installation)
  * [Prerequisites](#prerequisites)
  * [Installing](#installing)
    + [Users](#users)
    + [Developers](#developers)
- [Documentation](#installation)
- [Usage](#usage)
  * [NEMA](#nema)
  * [Vegetative stages](#vegetative-stages)
  * [Scenarios monoculms](#scenarios-monoculms)
- [Credits](#credits)
  * [Authors](#authors)
  * [Contributors](#contributors)
  * [Funding](#funding)
- [License](#license)


# Installation

## Prerequisites
*CNW-Grass* has the following dependencies, which are automatically installed (see [Installing](#installing)):

* [openalea.MTG](https://github.com/openalea/mtg)
* [openalea.Plantgl](https://github.com/openalea/plantgl)
* [openalea.Lpy](https://github.com/openalea/lpy)
* [openalea.Caribu](https://github.com/openalea/caribu) 
* [openalea.Adel](https://github.com/openalea/adel)
* [openalea.Astk](https://github.com/openalea/astk)
    
## Installing
For general information about OpenAlea installation, see https://openalea.readthedocs.io/en/latest/install.html 
### Users

```shell
conda create -n cnwgrass openalea.cnw-grass -c conda-forge -c openalea3
```
To activate the environment: `conda activate cnwgrass`

### Developers

1) To clone the project, please use:
```commandline
git clone https://github.com/openalea/cnwgrass
```

2) Move to the cloned directory, then create and activate a conda environment with dependencies:
```commandline
conda env create -n cnwgrass -f conda/environment.yml 
activate cnwgrass
```
# Documentation
https://cnwgrass.rtfd.io

# Usage

To date, *CNW-Grass* has been used in four main contexts described below. 

The scripts to run *CNW-Grass* are located in:
* `cnwgrass\example\NEMA`
* `cnwgrass\example\Papier_FSPMA2016`
* `cnwgrass\example\Vegetative_stages`
* `cnwgrass\example\Scenarii_monoculms`

## external soil model
This example illustrates the coupling of CNW-Grass with an external soil model, which provides root N uptake as an input.
The initial conditions and weather data are the same as "Vegetative_stages".
No coupling with hydraulics implemented.

In this example, to mimic the coupling with an external soil model, the [main.py](example/Farquhar_standalone/main.py)root N uptake at each time step 
is provided as an input in `cnwgrass\example\external soil model\inputs\nitrates_uptake_forcings.csv` 

To run the example:
* Open a command line interpreter in `cnwgrass\example\external soil model`
* To run the simulation, use : `python main.py`

## Farquhar_standalone
An example to show how to initialize and run the model Gas-Exchange in a standalone version.
The example runs Gas-Exchange with two different options : either with or without a coupling with a hydraulic model.

To run the example:
* Open a command line interpreter in `cnwgrass\example\Farquhar_standalone`
* To run the simulation, use : `python main.py`

## NEMA
This example deals with the post-flowering stages of wheat development under 3 nitrogen fertilisation regimes (H0, H3 and H15). The main processes described are leaf senescence, C and N remobilisation, grain filling). During that stages, all vegetative organs have completed their growth. 
This work led to the research articles [Barillot *et al.* (2016a)](https://doi.org/10.1093/aob/mcw143) and [Barillot *et al.* (2016b)](https://doi.org/10.1093/aob/mcw144).

To run the example:
* Open a command line interpreter in `cnwgrass\example\NEMA`
* To run the 3 scenarios, use : `python main.py`

## Scenarios monoculms
This example explores the plasticity of leaf growth during the vegetative stages of wheat development. 
The growth of wheat monoculms was simulated for highly contrasting conditions of soil nitrogen concentration, incident light and planting density.
The list of scenarios and their characteristics are specified in  `cnwgrass\example\Scenarios_monoculms\inputs\scenarios_list.csv`
This work led to the research article [Gauthier *et al.* (2021)](https://doi.org/10.1093/insilicoplants/diab034).
The original outputs as well as a singularity container with the code version used for the paper can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5503312.svg)](https://doi.org/10.5281/zenodo.5503312)
 
To run the example:
* Open a command line interpreter in `cnwgrass\example\Scenarios_monoculms`
* List the scenarios you want to run in the script *main.py* (scenario id should match with those listed in the *scenarios_list.csv* file) 
* Run the script *main.py*: `python main.py`
* The whole set of scenarios was run in the high-performance computing center [MESO@LR](https://meso-lr.umontpellier.fr/) (Université de Montpellier, France) 

## Vegetative_stages
This example simulates the early vegetative stages of wheat growth as measured from a field experiment conducted in 1998-99 in Grignon (France). It mainly covers the processes of leaf, internode and roots growth.
Tillering is simplified: tiller emergence is a model input while tiller metabolism and growth is approximated from that  of the main stem.
This work led to the research article [Gauthier *et al.* (2020)](https://doi.org/10.1093/jxb/eraa276). Results were obtained from the tag [paper_JXBot_2020](https://github.com/openalea/cnwgrass/releases/tag/paper_JXBot_2020).
The simulation starts at leaf 4 emergence on December 1998 and finishes on April 1999 at the beginning of internode elongation.
The model only accounts for C-N relations, there is no effect of water (hydraulics=False).
 
To run the example:
* Open a command line interpreter in `cnwgrass\example\Vegetative_stages`
* Run script *main.py*: `python main.py`

## Vegetative_stages_hydraulics
This example simulates a field experiment on winter wheat (cv Soissons) described in Ljutovac (2002).
The simulation starts at leaf 4 emergence on December 1998 and finishes on April 1999 at the beginning of internode elongation.
This example integrates the co-regulation of leaf growth by the trophic and hydraulic status.
This work led to the research article [Acker *et al.* (2026)](https://doi.org/10.1093/jxb/erag248).
The original outputs as well as a singularity container with the code version used for the paper can be found at https://doi.org/10.57745/W850YH
To further illustrate the model behaviour, the example also described how the user can trigger a drought event.

To run the example:
* Open a command line interpreter in `WheatFspm\example\Vegetative_stages_hydraulics`
* Run script *main.py*: `python main.py`


# Credits
## Authors
* Romain BARILLOT - model designing, development and validation - [rbarillot](https://github.com/rbarillot)
* Marion GAUTHIER - model designing, development and validation - [mngauthier](https://github.com/mngauthier)
* Victoria ACKER - model designing, development and validation - [victoriaacker](https://github.com/victoriaacker)
* Camille CHAMBON - software designing, development, deployment and optimization - [cachambon](https://github.com/cachambon)
* Bruno ANDRIEU - model designing and validation, scientific project management - [bandrieu](https://orcid.org/0000-0002-7933-9490)

## Contributors
* Christian FOURNIER - [christian34](https://github.com/christian34)
* Christophe PRADAL - [pradal](https://github.com/pradal)

## Funding
* [INRAE](https://www.inrae.fr/): salaries of permanent staff 
* French Research National Agency: projects [Breedwheat](https://breedwheat.fr/) (ANR-10-BTBR-03) and [Wheatamix](https://www6.inrae.fr/wheatamix/) (ANR-13-AGRO0008): postdoctoral research of R.Barillot
* [itk](https://www.itk.fr/en/) company and [ANRT](http://www.anrt.asso.fr/fr): funded the [Cifre](http://www.anrt.asso.fr/fr/cifre-7843) PhD thesis of M.Gauthier
* The Région Nouvelle-Aquitaine (France) and the [AgroEcoSystèmes](https://www.inrae.fr/departements/agroecosystem) Department of INRAE funded the PhD of V.Acker (Convention no. AAPR2021A-2020-11767810). 
  

# License
This project is licensed under the CeCILL-C License - see file [LICENSE](LICENSE) for details
 
