# Dominant Variant Finder

## Version History
### v0.1 - First release
* Proof of concept using 'imaginary data':
![Graph of imaginary data](dominant_variant_finder/Imaginary_Data_Graph.png)
* Imports dataset csv file
* Determines percentage of daily cases related to each variant (column)
* Calculate 7-day change in percentage of cases related to variant
* Identify variants of interest by two possible definitions/triggers
* Reports variants of interest with date that trigger was met
* Currently outputs graphs of variant cases Vs time 

![image](https://user-images.githubusercontent.com/81530279/113513443-afe02200-9561-11eb-8eb7-c383b8e346d7.png)
* Percentage contribution of variant Vs time

![image](https://user-images.githubusercontent.com/81530279/113513487-ee75dc80-9561-11eb-91ea-d0b39047eae3.png)
* 7day percentage change vs time

![image](https://user-images.githubusercontent.com/81530279/113513504-fa619e80-9561-11eb-93c6-d5d1fe26b49c.png)
* Prints results:
```
             Rapid Trigger At Threshold Met At  Alert Count
var_dominant       2020-05-14       2020-05-27            2
var_original              NaT       2020-03-18            1
```

## Known Issues/Further work
### Known Issues
* Definitions/triggers have been arbitrarily assigned after visual analysis of 'imaginary' data - would benefit from expert input!
### Further Work
* Import real dataset split by date and variant (possibly location too?)
* Adjust definitions/triggers to values more suited to real data (v0.1 output graphs should enable this analysis/inspection)
