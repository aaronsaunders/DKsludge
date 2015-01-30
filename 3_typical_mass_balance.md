# Typical Mass balance calculation


```r
source("R/functions.R")
```

```
## Loading required package: phyloseq
## Loading required package: ade4
## Loading required package: picante
## Loading required package: ape
## Loading required package: vegan
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.0-10
## 
## Attaching package: 'vegan'
## 
## The following object is masked from 'package:ade4':
## 
##     cca
## 
## Loading required package: nlme
## Note: the specification for S3 class "AsIs" in package 'BiocGenerics' seems equivalent to one from package 'RJSONIO': not turning on duplicate class definitions for this class.
## Note: the specification for S3 class "connection" in package 'BiocGenerics' seems equivalent to one from package 'RJSONIO': not turning on duplicate class definitions for this class.
## Note: the specification for S3 class "file" in package 'BiocGenerics' seems equivalent to one from package 'RJSONIO': not turning on duplicate class definitions for this class.
## Note: the specification for S3 class "pipe" in package 'BiocGenerics' seems equivalent to one from package 'RJSONIO': not turning on duplicate class definitions for this class.
## Note: the specification for S3 class "textConnection" in package 'BiocGenerics' seems equivalent to one from package 'RJSONIO': not turning on duplicate class definitions for this class.
## Loading required package: gridExtra
## Loading required package: grid
## 
## Attaching package: 'dplyr'
## 
## The following object is masked _by_ '.GlobalEnv':
## 
##     n
## 
## The following object is masked from 'package:nlme':
## 
##     collapse
## 
## The following objects are masked from 'package:plyr':
## 
##     arrange, desc, failwith, id, mutate, summarise
## 
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```


  - inf = influent
  - eff = effluent
  - WAS = waste activated sludge
  - AS  = activated sludge
  - WW  = wet weight
  - DW  = dry weight
  - OM  = organic matter
  - OM  = DW  -> inorganic matter assumed negligible
  - COD = Chemical oxygen demand
  - Vp   = total plant volume
  - Xp   = total plant biomass

## Define Constants


```r
# Nominal plant volume
Vp <- 1
# Hydrolyic retension time
HRT <- 1  # days
# Solids retension time
SRT <- 30  # days
# Assumed yield
Y_assumed <- 0.28

# unit conversion factors
um3_per_cm3 <- 1e+12
g_per_kg <- 1000
cm3_per_m3 <- 1e+06
```


## Dry weather influent


```r
# COD to organic matter ratio
COD_per_OM_inf <- 1.42  # kg-COD kg-OM^-1 [2]
# mean influent COD
COD_inf <- 0.582  # kg-COD m⁻³        [1] 
# mean influent cell concentration
cells_inf <- 1e+08  # cells cm⁻³        [1, 4] 
# influent flow
flow_inf <- Vp * HRT  # V day⁻¹

cells_inf <- cells_inf * cm3_per_m3  # cells m^-3 
OM_inf <- COD_inf/COD_per_OM_inf  # kg-OM m⁻³

# calculation of cell fraction of organic matter cell mass = 95 fg cell⁻¹
# [3] cell mass = 9.5e-13 g cell⁻¹

cell_volume <- 0.285  # µm³-cell cell⁻¹ [1] 
cell_density <- 1.1  # g-cell cm⁻³     [1, 3] 
OM_per_cell_WW <- 0.22  # g-OM g-cell⁻¹   [1, 3] 

cell_OM <- calculate_cell_OM(cell_volume, cell_density, OM_per_cell_WW)
cell_OM_inf <- cells_inf * cell_OM
# g-OM m⁻³ cell m⁻³ * g-OM cell⁻¹
cell_OM_inf <- cell_OM_inf/g_per_kg
# kg-OM m⁻³
cell_OM_percent <- round(cell_OM_inf/OM_inf * 100, 1)
```


## Waste activated sludge


```r
# mean AS cell concentration
cells_AS <- 1e+10  # cells cm⁻³  
cells_AS <- cells_AS * cm3_per_m3  # cells m⁻³ 
# WAS OM & COD
OM_AS <- 5  # kg-OM  m⁻³

COD_per_VSS_AS <- 1.42  # kg-COD kg-VSS⁻¹
VSS_per_OM_AS <- 0.75  # kg-VSS kg-OM⁻¹
COD_AS <- OM_AS * VSS_per_OM_AS * COD_per_VSS_AS  # kg-COD m⁻³

# WAS thickening
WAS_thick <- 3
OM_WAS <- OM_AS * WAS_thick  # kg-OM  m⁻³
COD_WAS <- COD_AS * WAS_thick  # kg-COD m⁻³

cells_WAS <- cells_AS * WAS_thick

# WAS flow
flow_WAS <- Vp * (1/SRT)/WAS_thick  # V day⁻¹
```


## Dry weather effluent


```r
# effluent flow
flow_eff <- Vp - flow_WAS  # V day⁻¹
# effluent OM & COD
OM_eff <- 20  #  g-OM m⁻³
OM_eff <- OM_eff/g_per_kg  # kg-OM m⁻³
COD_eff <- OM_eff * COD_per_OM_inf  # kg-COD m⁻³
# mean effleunt cell concentration
cells_eff <- 1e+07  # cells cm⁻³  [5] measured
# cells_eff <- 1.0e7 # cells cm⁻³ [4] measured
cells_eff <- cells_eff * cm3_per_m3  # cells m⁻³ 
```


## Mass balance


```r
# Cell budget
ncells_in <- cells_inf * flow_inf
ncells_eff <- cells_eff * flow_eff
ncells_AS <- cells_AS * Vp
ncells_WAS <- cells_WAS * flow_WAS
ncells_out <- ncells_eff + ncells_WAS

pprint <- function(x) {
    paste0(round(x * 100, 1), "%")
}
# % of out cells in effluent
pprint(ncells_eff/ncells_out)
```

```
## [1] "2.9%"
```

```r
# % of out cells in waste
pprint(ncells_WAS/ncells_out)
```

```
## [1] "97.1%"
```

```r
# fraction of out cells coming with WW
pprint(ncells_in/ncells_out)
```

```
## [1] "29.1%"
```

```r
# Daily effluent fraction of AS
pprint(ncells_eff/ncells_AS)
```

```
## [1] "0.1%"
```

```r
# Daily waste fraction of AS
pprint(ncells_WAS/ncells_AS)
```

```
## [1] "3.3%"
```



## References
  [1]: Vollertsen et al. 2001. Comparison of methods for determination of microbial biomass in wastewater
         
  [2]: Henze et al. 2002. 
  
  [3]: Loferer-Krössbacher et al. 1998. Determination of bacterial cell dry mass by transmission electron microscopy and densitometric image analysis.
         
  [4]: Snaidy. 2009. PhD Thesis. Detection and Enumeration of E. coli and Campylobacter in wastewater and treated water by CARD-FISH.
         
  [5]: Morgan-Sagastume 2008
  
  [6]: Frølund 1996
