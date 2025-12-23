# Data and R code for *Sixty years of plant community change in Europe indicate a shift towards nutrient-richer andshadier vegetation*

---

This repository includes the **data** and **R code** used in our study. It allows for the reproduction of the main analyses and supplementary materials. Due to their substantial file size, model outputs and predictions are not included in this repository but are currently available on Zenodo (<LINK>).

We studied temporal trends (1960-2020) for five EIVs (light, temperature, soil moisture, soil nitrogen, soil reaction) using data from 13,877 vascular plant taxa across 640,751 European vegetation plots. We calculated plot-level mean indicator values (CM<sub>EIV</sub>) and then applied Random Forests to interpolate their spatiotemporal dynamics.


## Table of Contents

* [Authors and Contact](#authors-and-contact)
* [Source Data (raw)](#source-data-raw-data-sources-not-deposited-here)
* [Data (data folder)](#data-cleaned-data-sources-data-folder)
* [R Code (src folder)](#r-code-r-scripts-src-folder)
    * [1. Random Forests tuning and final model fit](#1-random-forests-tuning-and-final-model-fit)
    * [2. Random Forests Cross-Validation](#2-random-forests-cross-validation)
    * [3. Random Forests diagnostics and evaluation](#3-random-forests-diagnostics-and-evaluation)
    * [4. Analyses on actual time series data (ReSurveyEurope)](#4-analyses-on-actual-time-series-data-resurveyeurope)
    * [5. Predict Random Forests (Interpolate ΔCMEIV)](#5-predict-random-forests-interpolate-%ce%b4cmeiv)
    * [6. Visualize interpolation results](#6-visualize-interpolation-results)
    * [7. Sensitivity analyses (CWM and NOTREES)](#7-sensitivity-analyses-using-community-weighted-means-cwm-and-excluding-treeshurb-species-notrees)
    * [Additional scripts](#additional-scripts)
* [License](#license)
* [Citation](#citation)

## Authors and Contact

### Contact
**Gabriele Midolo**  
Department of Spatial Sciences  
Faculty of Environmental Sciences  
Czech University of Life Sciences Prague, Praha-Suchdol, Czech Republic  
ORCID: [0000-0003-1316-2546](https://orcid.org/0000-0003-1316-2546)  
Email 1: `gabriele.midolo [at] gmail [dot] com`
Email 2: `midolo [at] fzp.czu [dot] cz`

### Authors and data contributors
**Gabriele Midolo** <a href="https://orcid.org/0000-0003-1316-2546" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Adam Thomas Clark <a href="https://orcid.org/0000-0002-8843-3278" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Milan Chytrý <a href="https://orcid.org/0000-0002-8122-3075" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Franz Essl <a href="https://orcid.org/0000-0001-8253-2112" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Stefan Dullinger <a href="https://orcid.org/0000-0003-3919-0887" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Ute Jandt <a href="https://orcid.org/0000-0002-3177-3669" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Helge Bruelheide <a href="https://orcid.org/0000-0003-3135-0356" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Jürgen Dengler <a href="https://orcid.org/0000-0003-3221-660X" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Irena Axmanová <a href="https://orcid.org/0000-0001-9440-7976" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Svetlana Aćić <a href="https://orcid.org/0000-0001-6553-3797" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Olivier Argagnon <a href="https://orcid.org/0000-0003-2069-7231" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Idoia Biurrun <a href="https://orcid.org/0000-0002-1454-0433" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Gianmaria Bonari <a href="https://orcid.org/0000-0002-5574-6067" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Alessandro Chiarucci <a href="https://orcid.org/0000-0003-1160-235X" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Renata Ćušterevska <a href="https://orcid.org/0000-0002-3849-6983" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Pieter De Frenne <a href="https://orcid.org/0000-0002-8613-0943" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Michele De Sanctis <a href="https://orcid.org/0000-0002-7280-6199" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Jan Divíšek <a href="https://orcid.org/0000-0002-5127-5130" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Jiří Doležal <a https://orcid.org/0000-0002-5829-4051" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Tetiana Dziuba <a href="https://orcid.org/0000-0001-8621-0890" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Rasmus Ejrnæs <a href="https://orcid.org/0000-0003-2538-8606" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Emmanuel Garbolino <a href="https://orcid.org/0000-0002-4954-6069" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Anke Jentsch <a href="https://orcid.org/0000-0002-2345-8300" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Borja Jiménez-Alfaro <a href="https://orcid.org/0000-0001-6601-9597" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Jonathan Lenoir <a href="https://orcid.org/0000-0003-0638-9582" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Jesper Erenskjold Moeslund <a href="https://orcid.org/0000-0001-8591-7149" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Francesca Napoleone <a href="https://orcid.org/0000-0002-3807-7180" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Sabine Rumpf <a href="https://orcid.org/0000-0001-5909-9568" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Jens-Christian Svenning <a href="https://orcid.org/0000-0002-3415-0862" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Grzegorz Swacha <a href="https://orcid.org/0000-0002-6380-2954" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Irina Tatarenko <a href="https://orcid.org/0000-0001-6835-2465" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Martin Večeřa <a href="https://orcid.org/0000-0001-8507-791X" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Denys Vynokurov <a href="https://orcid.org/0000-0001-7003-6680" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
Petr Keil <a href="https://orcid.org/0000-0003-3017-1858" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>


## Source Data: Raw Data Sources (*Not Deposited Here*)

The primary data sources used in this project are not directly stored in this repository but can be accessed as follows:

-   **Complete (species) vegetation-plot data:** These data are stored within the EVA database repository (<https://doi.org/10.58060/250x-we61>). Access can be requested through the EVA Coordinating Board (see <https://euroveg.org/eva-database/>).
-   **Ecological Indicator Values for Europe (EIVE 1.0):** These data are publicly available in Dengler et al. 2023: <https://doi.org/10.3897/VCS.98324.suppl8>.

## Data: Cleaned Data Sources (`data` folder)

The `data` folder contains cleaned and processed data files used for the analyses:

-   [`data/EUNIS_ESy2_habitat.names.txt`](data/EUNIS_ESy2_habitat.names.txt): A text file providing full names for the EUNIS-ESy level 2 habitat classification (for more details on EUNIS classification, refer to <https://floraveg.eu/habitat/>).
-   [`data/EVA.csv.xz`](data/EVA.csv.xz): A xz-compressed CSV file containing selected vegetation plots from the EVA database. 
-   [`data/ReSurveyEU_clean.csv.xz`](data/ReSurveyEU_clean.csv.xz): A xz-compressed CSV file with selected plots from the ReSurveyEurope dataset.
-   [`data/EVA_ReSu_CWM.csv.xz`](data/EVA_ReSu_CWM.csv.xz): A xz-compressed CSV file containing a subset of plots for which the community-weighted mean (CWM) of EIVs were calculated (used for sensitivity analyses).
-   [`data/EVA_ReSu_NOTREES.csv.xz`](data/EVA_ReSu_NOTREES.csv.xz): An xz-compressed CSV file containing community means (CM) of EIVs, calculated with trees and shrubs excluded (used for sensitivity analyses).

Overview of the EVA and ReSurveyEurope data columns stored in this repository:

| Column Name | Description |
| :---------- | :----------------------------------------------------------------- |
| `database` | Database (`EVA` or `ReSurveyEU`) |
| `resurv_id` | Unique identifier for resurvey plot |
| `plot_id` | Unique identifier for plot observation (anonymized) |
| `ReSur_type` | Resurvey plot type (either `Permanent` or `Resampling`) |
| `dataset` | Dataset name (anonymized) |
| `ESy` | EUNIS habitat code (level 3, where possible) |
| `habitat` | Level 1 EUNIS habitat category (`Forest`, `Grassland`, `Scrub`, or `Wetland`) |
| `lon` | Longitude in WGS84 |
| `lat` | Latitude in WGS84 |
| `x` | Longitude in EPSG:25832 |
| `y` | Latitude in EPSG:25832 |
| `elev` | Elevation above sea level (m) |
| `year` | Year of sampling |
| `plot_size` | Plot size in squared meters |
| `n` | Number of species (species richness) |
| `x_mean` | Centroid longitude (in EPSG:25832) of plot observations assigned to the same `resurv_id` |
| `y_mean` | Centroid latitude (in EPSG:25832) of plot observations assigned to the same `resurv_id` |
| `max_dist_m` | Maximum distance (in meters) between plot observations assigned to the same `resurv_id` |
| `n.EIV_M` | Proportion of species with available EIV for moisture |
| `n.EIV_N` | Proportion of species with available EIV for nitrogen |
| `n.EIV_R` | Proportion of species with available EIV for reaction |
| `n.EIV_L` | Proportion of species with available EIV for light |
| `n.EIV_T` | Proportion of species with available EIV for temperature |
| `cm.EIV_M` | Community mean EIV for moisture |
| `cm.EIV_N` | Community mean EIV for nitrogen |
| `cm.EIV_R` | Community mean EIV for reaction |
| `cm.EIV_L` | Community mean EIV for light |
| `cm.EIV_T` | Community mean EIV for temperature |

\* = N.B. plot_id in the shared data is not the same as the one stored in models outputs in the `models` and `preds` folder! The former was generated to anonymize original data in EVA and ReSurvey Europe, but results from models and predictions stores the original values.

In addition, [`data/EVA_ReSu_CWM.csv.xz`](data/EVA_ReSu_CWM.csv.xz) contains the following columns (community weighted means calculated using species cover in each plot):

| Column Name | Description |
| :---------- | :----------------------------------------------------------------- |
| `cwm.EIV_M` | community weighted mean EIV for moisture |
| `cwm.EIV_N` | community weighted mean EIV for nitrogen |
| `cwm.EIV_R` | community weighted mean EIV for reaction |
| `cwm.EIV_L` | community weighted mean EIV for light |
| `cwm.EIV_T` | community weighted mean EIV for temperature |

## R Code: R Scripts (`src` folder)

The `src` folder contains the R scripts organized by their analytical purpose:

### 1. Random Forests tuning and final model fit

- [`1_tuning.R`](src/1_tuning.R): Script for tuning and fitting the Random Forest model using all available observations.
- [`1_tuning_EUNIS.lev.2.R`](src/1_tuning_EUNIS.lev.2.R): Script for tuning and fitting the Random Forest model exclusively on plot data with available EUNIS-ESy level 2 habitat classifications.

### 2. Random Forests Cross-Validation

- [`2_kfoldcv.R`](src/2_kfoldcv.R): Performs 10-fold cross-validation on the training data.

### 3. Random Forests diagnostics and evaluation

- [`3_modeldiagnostic.R`](src/3_modeldiagnostic.R): Generates figures and tables for model evaluation, including tuning results and other statistical summaries.

### 4. Analyses on actual time series data (ReSurveyEurope)

- [`4_lme_resurvey.R`](src/4_lme_resurvey.R): Fit and visualize linear mixed-effects models on the time series data.
- [`4_validation.R`](src/4_validation.R): Validates the interpolation approach using the time series data.

### 5. Predict Random Forests (Interpolate ΔCM<sub>EIV</sub>)

- [`5_interpolation.R`](src/5_interpolation.R): Interpolates changes in community-mean ecological indicator values (ΔCM<sub>EIV</sub>). Uses prediction ensemble (mean) across all trees in the random forests.
- [`5_raw_interpolation.R`](src/5_raw_interpolation.R): Interpolates average ΔCM<sub>EIV</sub>s across all trees in the random forests and store summary stats for predictions across all plots for each year (1960-2020) and habitat type.

### 6. Visualize interpolation results

- [`6_plot_trends.R`](src/6_plot_trends.R): Visualize histograms and partial plots illustrating the interpolated dynamics of CM<sub>EIV</sub>.
- [`6_plot_trends_EUNIS.lev.2.R`](src/6_plot_trends_EUNIS.lev.2.R): Visualzie partial plots illustrating the interpolated dynamics of CM<sub>EIV</sub> across different EUNIS-ESy level 2 habitats.
- [`6_map.geo.R`](src/6_map.geo.R): Visualize geographical maps displaying the average ΔCM<sub>EIV</sub> across Europe.

### 7. Sensitivity analyses using Community Weighted Means ('CWM') and excluding tree/shurb species ('NOTREES')

- [`6_plot_trends.R`](src/6_plot_trends.R): Visualize histograms and partial plots illustrating the interpolated dynamics of CM<sub>EIV</sub>.
- [`6_plot_trends_EUNIS.lev.2.R`](src/6_plot_trends_EUNIS.lev.2.R): Visualzie partial plots illustrating the interpolated dynamics of CM<sub>EIV</sub> across different EUNIS-ESy level 2 habitats.

### Additional scripts

- [`0_anonymize_and_clean.R`](src/0_anonymize_and_clean.R): An internal script used to anonymize the `plot_id` column from EVA/ReSurveyEurope data shared in this repository (no need to run this script).
- [`0_helpfunctions.R`](src/0_helpfunctions.R): Contains R functions sourced by other analysis scripts.

## License

**Data** are available under the terms of the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International license (CC BY-NC-ND 4.0) (<https://creativecommons.org/licenses/by-nc-nd/4.0/>).

**Code** is available under the terms of the GNU General Public License v3.0 (GPL-3.0) (<https://www.gnu.org/licenses/gpl-3.0.html>).

## Citation

*This repository is not linked to any publication yet.*

---
