# WWF Waterhole mapping project in KAZA

Margaret Swift, Robin Naidoo, Steve Osofsky, Shirley Atkinson

# Project steps
The general steps for the project. For 1,2, and 4 I anticipate finding a "good enough" method that allows me to proceed with the next step; I'll then return later to fine-tune the (1) segmentation, (2) cloud masking, and (4) maximum-wetness water mask.

## 1.	Segment KAZA into subregions
The KAZA region is home to many different habitats and landscape types that will determine how waterholes are best mapped in that subregion. Using an algorithm TBD, I will break KAZA up into subregions, ideally contiguous, where in-group autocorrelation is maximized and between-group autocorrelation is minimized.

## 2.	Cloud masking
Clouds are a significant issue in most remote sensing applications, but especially when trying to map small water bodies. I need to find the optimal cloud-removal method for these data (Sentinel-2)

## 3. Calculating Otsu thresholds for each time period
As in 2022, I'll break up the year into six time periods: Each month of June, July, and August; tri-monthly from Sept-May (Sept-Nov; Dec-Feb; March-May). For each period _P_ and within each subregion _R_, I will take the best cloud-free imagery, add AWEI, and create a single median image across all years in the collection. I will then use the Otsu method on this across-year median image to find an optimal AWEI threshold, and save these thresholds for each {_P, R_}.

## 4.	Creating maximum wetness water mask
We want to track waterhole size over time, so creating a binary mask of waterhole-notwaterhole will be important. First, I will pull Sentinel-2 images from February - May of each year, add AWEI, and create a median image for each year. Using the threshold generated in step 3, I will then threshold each within-year median image, then sum those images across years. This summed image will have values from 0 to N where N is the number of years. After all of this processing, I will save (1) the summed image for future reference as a measure of water frequency, and (2) a binary layer where every value > 0 is labeled as potentially filled with water.

## 5. Mapping water across years and regions
Finally, for each region _R_, time period _P_, and year _Y_ I will create median AWEI images, mask these images using the binary mask in step 4, and threshold using the period-specific Otsu threshold developed in step 3. 

## 6. Testing and validation
Once I have all the maps from part 5, I can start validating how well this process works against real data (validation dataset created from Planet or some other high-resolution source). Statistics will include:
1. Commission error rate.
2. Omission error rate.


        
