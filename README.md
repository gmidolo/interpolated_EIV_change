# Data and R code for *Sixty years of plant community change in Europe indicate a shift towards nutrient‑richer and denser vegetation*

---

This repository includes the **data** and **R code** used in our study. It allows for the reproduction of the main analyses and supplementary materials. Due to their substantial file size, model outputs and predictions are not included in this repository but can be downloaded from Zenodo (DOI: [10.5281/zenodo.18266902](https://doi.org/10.5281/zenodo.18266902)).

We studied temporal trends (1960-2020) for five EIVs (light, temperature, soil moisture, soil nitrogen, soil reaction) using data from 13,874 vascular plant taxa across 692,393 vegetation plots in Europe. We calculated plot-level mean indicator values (CM<sub>EIV</sub>) and then applied Random Forests to interpolate their plot-level spatiotemporal dynamics.

The [interactive map to explore interpolated spatiotemporal changes of CM<sub>EIV</sub>s](https://gmidolo.shinyapps.io/interpolated_EIV_change_app) is deposited here: [https://github.com/gmidolo/interpolated_EIV_change_app](https://github.com/gmidolo/interpolated_EIV_change_app). 

## Table of Contents

* [Authors and Contact](#authors-and-contact)
* [Source Data (raw)](#source-data-raw-data-sources-not-deposited-here)
* [Data (data folder)](#data-cleaned-data-sources-data-folder)
* [R Code (src folder)](#r-code-r-scripts-src-folder)
    * [1. Random Forests tuning and final model fit](#1-random-forests-tuning-and-final-model-fit)
    * [2. Random Forests Cross-Validation](#2-random-forests-cross-validation)
    * [3. Random Forests diagnostics and evaluation](#3-random-forests-diagnostics-and-evaluation)
    * [4. Analyses on actual time series data (ReSurveyEurope)](#4-analyses-on-actual-time-series-data-resurveyeurope)
    * [5. Predict Random Forests (interpolate CM<sub>EIV</sub>)](#5-predict-random-forests-interpolate-%ce%b4cmeiv)
    * [6. Visualize interpolation results](#6-visualize-interpolation-results)
    * [7. Sensitivity analyses](#7-sensitivity-analyses-using-community-weighted-means-cwm-and-excluding-treeshurb-species-notrees)
    * [Additional scripts](#additional-scripts)
    * [Scripts for raw data processing](#raw-data-processing-(not-reproducible))
* [Metadata](#metadata)
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
Gabriele Midolo <a href="https://orcid.org/0000-0003-1316-2546" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
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
Jiří Doležal <a href="https://orcid.org/0000-0002-5829-4051" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" class="is-rounded" width="15"/></a>,
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

-   [`data/EUNIS_ESy2_habitat.names.txt`](data/EUNIS_ESy2_habitat.names.txt): A text file providing full names for the level-2 EUNIS-ESy habitat classification (for more details on EUNIS classification, refer to <https://floraveg.eu/habitat/>).
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

- [`1_tuning.R`](src/1_tuning.R): Tune and fit Random Forest models using all available observations.
- [`1_tuning_EUNIS.lev.2.R`](src/1_tuning_EUNIS.lev.2.R): Tune and fit Random Forest models using level-2 EUNIS-ESy habitats as predictors.

### 2. Random Forests Cross-Validation

- [`2_kfoldcv.R`](src/2_kfoldcv.R): Performs 10-fold cross-validation on the training data.

### 3. Random Forests diagnostics and evaluation

- [`3_modeldiagnostic.R`](src/3_modeldiagnostic.R): Generates figures and tables for model evaluation, including tuning results and other statistical summaries.

### 4. Analyses on actual time series data (ReSurveyEurope)

- [`4_lme_resurvey.R`](src/4_lme_resurvey.R): Fit and visualize linear mixed-effects models on the time series data.
- [`4_validation.R`](src/4_validation.R): Validates the interpolation approach using the time series data.

### 5. Predict Random Forests (Interpolate ΔCM<sub>EIV</sub>)

- [`5_interpolation.R`](src/5_interpolation.R): Interpolate changes in community-mean ecological indicator values (ΔCM<sub>EIV</sub>). Uses prediction ensemble (mean) across all trees in the random forests.
- [`5_raw_interpolation.R`](src/5_raw_interpolation.R): Interpolate average CM<sub>EIV</sub>s across all trees in the random forests and store summary stats for predictions across all plots for each year (1960-2020) and main habitat type.
- [`5_raw_interpolation_EUNIS.lev.2.R`](src/5_raw_interpolation_EUNIS.lev.2.R): Interpolate average CM<sub>EIV</sub>s across all trees in the random forests and store summary stats for predictions for level-2 EUNIS-ESy habitats.

### 6. Visualize interpolation results

- [`6_plot_trends.R`](src/6_plot_trends.R): Visualize histograms and partial plots illustrating the interpolated dynamics of CM<sub>EIV</sub>.
- [`6_plot_trends_EUNIS.lev.2.R`](src/6_plot_trends_EUNIS.lev.2.R): Visualize partial plots illustrating the interpolated dynamics of CM<sub>EIV</sub> across different level-2 EUNIS-ESy habitats.
- [`6_map.geo.R`](src/6_map.geo.R): Visualize geographical maps displaying the average ΔCM<sub>EIV</sub> across Europe.

### 7. Sensitivity analyses using Community Weighted Means ('CWM') and excluding tree/shurb species ('NOTREES')

- [`7_tuning_CWM.R`](src/7_tuning_CWM.R): Tune and fit separate Random Forest models for both community means (CM<sub>EIV</sub>) and community weighted means(CWM<sub>EIV</sub>).
- [`7_tuning_NOTREES.R`](src/7_tuning_NOTREES.R): Tune and fit Random Forest models for CM<sub>EIV</sub>s calculated by excluding tree and shrub species.
- [`7_interpolate_and_plot_CWM.R`](src/7_interpolate_and_plot_CWM.R): Interpolate and plot results from [`7_tuning_CWM.R`](src/7_tuning_CWM.R).
- [`7_interpolate_and_plot_TREES.R`](src/7_interpolate_and_plot_TREES.R): Interpolate and plot results from [`7_tuning_NOTREES.R`](src/7_tuning_NOTREES.R).

### Additional scripts

- [`0_anonymize_and_clean.R`](src/0_anonymize_and_clean.R): An internal script used to anonymize the `plot_id` column from EVA/ReSurveyEurope data shared in this repository (no need to run this script).
- [`0_helpfunctions.R`](src/0_helpfunctions.R): Contains R functions sourced by other analysis scripts.

### Raw data processing (not reproducible)

**N.B. This part of the code is used to prepare the raw data from EVA and ReSurveyEurope. It is not directly reproducible unless a data request is made to the EVA and ReSurveyEurope governing board.**

- [`raw_data_processing/1_prepare.data.R`](src/raw_data_processing/1_prepare.data.R): Main script to preprocess raw data retrieved in EVA/ReSurveyEurope proj. no. 222 (DOI: [ 10.58060/250x-we61](https://doi.org/10.58060/250x-we61))
- [`raw_data_processing/2_search.and.remove.duplicates.R`](src/raw_data_processing/2_search.and.remove.duplicates.R): Remove presumed or actual duplicate plots within and across EVA and ReSurveyEurope data (plots with the same year of sampling, geographic coordinates, and species composition)
- [`raw_data_processing/3_clean_ReSurveyEU.R`](src/raw_data_processing/3_clean_ReSurveyEU.R): Remove plots in ReSurveyEurope with 'uncertain' plots location (i.e., $\ge$ 100 m distance between observations of the same plots)
- [`raw_data_processing/4_finalcleanforCWM.R`](src/raw_data_processing/4_finalcleanforCWM.R): Correct plot selection made in previous steps for the subset of data containing community weighted means (CWM)
- [`raw_data_processing/5_noTREES.and.rare.species`](src/raw_data_processing/5_noTREES.and.rare.species): Recalculate CM<sub>EIV</sub>s by excluding trees and shurb species; compare sensitivity of CM<sub>EIV</sub>s when "rare" species are excluded

## Metadata

### 1. General Dataset Information

|  | **Description** |
|-------------------------------|-----------------------------------------|
| **Dataset names** | European Vegetation Archive (EVA); ReSurveyEurope |
| **Version** | Version 2024-09-19 (DOI: <https://doi.org/10.58060/hgrb-sw46>) |
| **Project name** | "EVA project \# 222 – 2024-09-12 Interpolated dynamics of local plant diversity in European vegetation - G. Midolo \| SELECTION 2024-10-31" (DOI: <https://doi.org/10.58060/250x-we61>) |
| **Date of creation** | Data selection date for project #222: 2024-10-31. The EVA database is in development since 2012 and first made available for use in research projects in 2014. The first data call for ReSurveyEurope was announced in 2020. |
| **Citation** | Chytrý, M., Hennekens, S. M., Jiménez‐Alfaro, B., Knollová, I., Dengler, J., Jansen, F., ... & Yamalov, S. (2016). European Vegetation Archive (EVA): an integrated database of European vegetation plots. *Applied Vegetation Science*, *19*(1), 173-180. <https://doi.org/10.1111/avsc.12191> <br> Knollová, I., Chytrý, M., Bruelheide, H., Dullinger, S., Jandt, U., Bernhardt‐Römermann, M., ... & Essl, F. (2024). ReSurveyEurope: A database of resurveyed vegetation plots in Europe. *Journal of Vegetation Science*, *35*(2), e13235. <https://doi.org/10.1111/jvs.13235> |
| **Data curators** | [The EVA Coordinating Board and EVA Council](https://euroveg.org/eva-database/who-we-are) <br> [The ReSurveyEurope Board](https://euroveg.org/resurvey/) |

------------------------------------------------------------------------

### 2. Dataset Description

|   | **Raw data (\*)** | **Data in the repository** |
|-----------------|----------------------|----------------------------------|
| **Summary** | EVA contains vegetation plot data across Europe, including species composition and plot metadata. ReSurveyEurope contains data in a similar format, but for plots with repeated measurements over time. | Subset restricted to plot observations relevant for analyzing temporal changes in community means of ecological indicator values of vascular plants across Europe. Includes processed and harmonized species and site-level data used for modelling community-mean ecological indicator value change. |
| **Provenance** | Compiled from 308 datasets contributed by data owners under EVA and ReSurveyEurope governance. | Compiled from 270 datasets contributed by data owners under EVA and ReSurveyEurope governance. |
| **Temporal coverage (range)** | 1873–2023 | 1945–2023 (model training); 1960–2020 (predictions/interpolation) |
| **Geographical coverage (range)** | Longitude (WGS84): -180.00 – 64.84; Latitude (WGS84): -90.00 – 80.15 (includes geographic outliers) | Longitude (WGS84): -10.52 – 38.78; Latitude (WGS84): 34.80 – 71.12 |
| **Sampling frame** | Vegetation plots categorized as various European habitats: marine, coastal, inland water, wetland, grassland, scrub, forest, inland habitats with little soil, and vegetated man-made habitats, as defined in the EUNIS Habitat Classification System. | Vegetation plots categorized exclusively as forest, grassland, scrub, and wetland vegetation, as defined in the EUNIS Habitat Classification System. |
| **Number of records and plots** | No. of plots: 1,745,721; core EVA: 1,676,182; ReSurveyEurope plot observations: 103,397 | No. of core EVA: 622,906 vegetation plots; core ReSurveyEurope plot observations: 69,487 (from 21,618 resurvey plots) |
| **Number of variables** | 46 fields in the header metadata (raw) | 28 fields in the shared data (after merging EVA and ReSurveyEurope) |

(\*) These values refer to the raw data released for project \# 222 (<https://doi.org/10.58060/250x-we61>)

------------------------------------------------------------------------

### 3. Access and Licensing

|   | **Raw data (\*)** | **Data in the repository** |
|-----------------|----------------------|----------------------------------|
| **Access conditions** | Data access follows either free-access, semi-restricted, or restricted models, depending on the availability regime assigned by the database custodians. In any case, **the data are only accessible through a data request to the EVA and ReSurveyEurope Governing Board.** | Data access follows either free-access, semi-restricted, or restricted models, depending on the availability regime assigned by the database custodians. The data stored in this repository can be used to reproduce the main analyses (models, predictions, figures). **Use in other works or publications requires Governing Board approval.** |
| **How to obtain access** | Data can be accessed through a request to the Governing Board, following [Article 5 of the EVA rules](https://euroveg.org/download/eva-rules.pdf). More info: <https://euroveg.org/eva-database/obtaining-data> | Through this repository ([`input`](data/input)). **To use these data elsewhere, a data request must be submitted to the Governing Board**. |
| **License / Terms of Use** | Depending on the availability regime assigned by the database custodians and EVA/ReSurveyEurope Governing Board. | Derived data follow the same access constraints. Metadata and code in this repository are licensed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). |
| **Data availability statement** | Not publicly accessible due to third-party ownership. | All data-processing scripts and data that are essential to reproduce the analyses are publicly available in this repository. |

(\*) These values refer to the raw data released for project \# 222 (<https://doi.org/10.58060/250x-we61>)

------------------------------------------------------------------------

### 4. Methods and Processing

|   | **Raw data (\*)** | **Data in the repository** |
|-----------------|----------------------|----------------------------------|
| **Data collection methods** | Records of the abundance and/or occurrence of plant species found in vegetation plots collected in the field. | Not available (data are elaborated based on the raw data). |
| **Inclusion / exclusion criteria** | Includes different vegetation types, plot sizes, quality levels, and time periods available in the raw data. Eligibility of datasets to be included in EVA and ReSurveyEurope is regulated by the Governing Board. | Vegetation plots with full lists of vascular plant species categorized as forest, grassland, scrub, and wetland, with plot size 1–1000 m², valid geographic coordinates, and sampled between 1945–2023. Detailed selection criteria are in the main manuscript. |
| **Data harmonization** | Managed by the EVA and ReSurveyEurope curators (species abundances and nomenclature). | Same as 'raw data'. |
| **Data processing (main) steps** | Managed by EVA and ReSurveyEurope curators. | 1\. Application of inclusion/exclusion criteria to the raw data. <br> 2. Calculation of community-mean ecological indicator value (CM<sub>EIV</sub>) for each plot using the species list data from the raw data. |
| **Quality control** | Conducted by data providers and EVA/ReSurveyEurope curators. | Conducted by Gabriele Midolo with help from database custodians. Included coordinate plausibility checks (manual and automated), species richness plausibility checks (manual), and duplicate removal within and between EVA and ReSurveyEurope (automated). |
| **Processing code** | Not available. | All R scripts documenting processing steps and quality control are stored in this repository ([src/raw_data_processing](src/raw_data_processing)). |
| **Software environment** | Varies by data provider. | R version 4.4.2. R packages for data processing: `tidyverse`, `sf`, `terra`. |
| **Reproducibility** | Not reproducible without subset definition. | Fully reproducible for EVA-authorized users using the R code stored in this repository. |

(\*) These values refer to the raw data released for project \# 222 (<https://doi.org/10.58060/250x-we61>)

## License

**Data** are available under the terms of the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International license (CC BY-NC-ND 4.0) (<https://creativecommons.org/licenses/by-nc-nd/4.0/>).

**Code** are available under the terms of the GNU General Public License v3.0 (GPL-3.0) (<https://www.gnu.org/licenses/gpl-3.0.html>).

## Citation

*This repository is not linked to any publication yet.*