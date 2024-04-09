# Trend-screening

This repository contains R-scripts and functions to visualise trend for multiple series simultaneously. Main code is found in 'Screening_plots_for_article.Rmd' where the dataset 'vdr 88-17 ANC300.txt' is analysed. The available functions are found in 'sourcecode_screening2020.R'. An additional analysis on a smaller regional scale is described in 'Hierachical GAM for the region of Värmland.Rmd'. 

Update 9 April 2024: Codes are updated to reflect changes in the package gratia (version 0.9.0).

Update 15 September 2023: A new sourcecode file is available to match changes in the gratia package: sourcecode_screening2023.R. Also a sourcecode file that works with annual rather than seasonal data is available (sourcecode_screening2023_annual_data.R. 

The following functions are available: 

screeningmodeling()

  The function screeningmodeling() fits generalized additive models to individual series in the dataset. Data must be provided in long format with different stations (here     referred to as id) and variables (referred to as name) presented consecutively in a single column.           

  screeningmodeling(.data, datevar, values, link = "identity", autocor = TRUE, conf.type = "conf", tdist = FALSE, beep = FALSE,…)

  Arguments: 
	  values  	the variable that contains the measurements
    datevar	  a variable containing dates in an approved date format, e.g. as decimal numbers
	  link		  the chosen link function (only Gaussian family is supported)
    autocor	  indicates whether autocorrelation should be modelled (TRUE) or not (FALSE)
    conf.type	confidence interval computed as pointwise ("confidence") or simultaneous ("simultaneous")  following the package gratia
    tdist	    indicates whether the t-distribution should be used to model the error distribution (TRUE), employing a Gaussian distribution as default; this can be done only if the autocor=FALSE option is chosen
    beep	    whether  an audible signal should be emitted  when computations are     complete  (TRUE) or not (FALSE)
	  ….		    variables to be passed to the program, e.g. site id and variable names


  The output data set provided contains model fits with all necessary information nested within each of the passed variables (e.g., id and name). This output is needed to run the following plot functions.  


plot_individual_trend()

  The function plot_individual_trend(.output, y=NULL, title=NULL) uses the output of the function screeningmodeling() to plot individual series including color coding for  periods with significant trends. To apply this function, a single series (site, variable) must be selected before the function is run on that subset. The y-axis title and the plot title can be specified.


plot_screeningtrends() and plot_screeningtrends_pvalues

  The function plot_screeningtrends() plots trends using color coding (i.e., yellow [non-significant], red [increasing], or blue [decreasing]) for multiple series ordered by the  variable given in sorting=. For example, id, name or intercept can be used as sorting variable. The variable intercept is produced by screeningmodeling() and contains the      intercept of the individual models. If the model used contains only spline terms and no parametric functions for explanatory variables, the intercept corresponds to the overall mean of the individual series. It can then be used to sort the plot by mean levels of the response. New variables can also be created outside this function, added to the output, and used as sorting variables. One value per station has to be given. The same rules apply to the wrapper variable (wrappingvar=) that is used to define blocks of trends (e.g., variable names as used in Figs. 3–5).

  plot_screeningtrends(.output, sorting = id, y_id = id, wrappingvar = name) 

  Arguments: 
	  sorting	  	a variable according to which the series should be sorted
	  y_id		    an identification variable printed along with each series
	  wrappingvar a variable that defines groups that belong together


  A similar plot can be produced with color coding according to computed p-values. The input variables are identical. 

  plot_screeningtrends_pvalues(.output, sorting = id, y_id = id, wrappingvar = name)


plot_proportions()

  A proportion plot enables summarization of the results of many series in a single plot that shows at any specific time point how many of the trends are increasing, decreasing, or non-significant. The wrapping variable works in the same way as in the screening plots.

  plot_proportions(.output, wrappingvar = name)


plot_screeningtrends_reference() and plot_screeningtrends_relative()

  Two different functions are available for visualizing the magnitude of the trend. The function compares point-wise levels of the estimated trend function with the mean of a reference period. At present, the reference period is set as the first three years of the observed period for each site respectively. Input variables are the same as described for other screening plots: 
  
  plot_screeningtrends_reference(.output, sorting = id, y_id = id, wrappingvar = name)

  The function plot_screeningtrends_relative() can be applied to compare estimated magnitudes of change at a certain time point with the estimated level of the variable at the same time point. This function can be useful to identify episodes of high values and potentially also for data quality control. The input variables are again the same as above:   
  
plot_screeningtrends_relative(.output, sorting = id, y_id = id, wrappingvar = name)
