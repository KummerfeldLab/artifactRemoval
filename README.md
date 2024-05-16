


## This is a python package for the purpose of tissue artifacts detection (Edge, Border, abnormal spots).




## Current version 

## Installation 

### Python
```
pip install git+https://github.com/KummerfeldLab/artifactRemoval
```


### R
Due to many existing software are already developed in R. We provided a tutorial of applying package **'reticulate'** in R for run the same procedure as in python:  
```
install.packages("reticulate")
library(reticulate)
py_install("git+https://github.com/KummerfeldLab/artifactRemoval")
```



## Usage Tutorial


Here is an example of we read in a tissue sample and produce a cleaned up "tissue_position.csv" file. 

#### Python Version:
```
dir = "your tissue directory"
test = Artifact_remove(dir = dir)

test.remove_border() # remove border spots
test.remove_edge(distance=3) # remove edge spots, distance means all spots within 3 spots from edge
test.remove_malfunction() # remove malfunction spots
test.remove_edge(distance=6) # for illustration for repeatedly remove more edge spots.


test.review_removing() # review what we have done
test.simple_cleanser() # back to original tissue sample
test.save(dir) # save the cleaned tissue "tissue_position.csv" file to the director you set

```

#### R Version:
```
use_python("~/.virtualenvs/r-reticulate/bin/python")
use_virtualenv("r-reticulate")
py_available()
py_config()
py_install("git+https://github.com/KummerfeldLab/artifactRemoval")
py_install("scipy")
py_install("numpy")
py_install("pickle")
#version <- "3.12.3"
#install_python(version)

artifactRemoval <- reticulate::import("artifactsRemoval")
dir = "your tissue directory"


test <- artifactRemoval$Artifact_remove(dir)
test$remove_border() # remove border spots
test$remove_edge(distance=3) # remove edge spots, distance means all spots within 3 spots from edge
test$remove_malfunction() # remove malfunction spots
test$remove_edge(distance=6) # for illustration for repeatedly remove more edge spots.


test$review_removing() # review what we have done
test$simple_cleanser() # back to original tissue sample
test$save(dir) # save the cleaned tissue "tissue_position.csv" file to the director you set
```