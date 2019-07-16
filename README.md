# WaterRegionEstimationLandsat
Land cover classification and investigation of temporal changes are considered to be common applications of remote sensing. Water/non-water region estimation is one of the most fundamental classification tasks, analyzing the occurrence of water on the Earth’s surface. However, common remote sensing practices such as thresholding, spectral analysis, and statistical approaches are not sufficient to produce a globally adaptable water classification. The aim of this study is to develop a formula with automatically derived tuning parameters using perceptron neural networks for water/non-water region estimation, which we call the Perceptron-Derived Water Formula (PDWF), using Landsat-8 images.

Details of the work can be found at https://doi.org/10.3390/s18124333

**Steps:-**
1. Top of atmospheric reflectance 
2. Generate indices and water region estimation using  Perceptron derived parameters; 
 **water** = [0.989465, 1.14267147, 0.78721398, −0.93026412, −0.57805818], bias = [0.8181203]
 **non-water** = [−1.04869103, −1.17793739, −0.73774189, 1.03303862, 0.65516961], bias = [0.88329011]

3. Sunglint correction using sun zenith, sensor zenith, sun azimuth, and sensor azimuth 
4. Snow correction usins snow index

**Notes:**
Landsat-8 data can be downloaded, unzip and kept in the defined loc folder (here it is 'Landsat_data') also a text file represent the the list of the datas kep at loc folder (here it as 'list_.txt')

**Citation:**
Vinayaraj, P.; Imamoglu, N.; Nakamura, R.; Oda, A. Investigation on Perceptron Learning for Water Region Estimation Using Large-Scale Multispectral Images. Sensors 2018, 18, 4333.

