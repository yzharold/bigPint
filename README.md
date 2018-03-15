# bigPint

**Motivation:**

It is important to plot subsets of variables in order to examine variable associations in a dataset. Traditional modeling approaches without plotting the data are problematic because models are replete with assumptions that they alone cannot call into question. By plotting the data, analysts can improve modeling; they can iterate between visualizations and modeling to enhance the models based on feedback from the visuals.

Large multivariate datasets are common across numerous disciplinary fields. Unfortunately, the most popular visualization techniques for such data are often inadequate, if not misleading. The best approaches for looking at quantitative multivariate data are scatterplots of all pairs of variables, often laid out in a matrix format; parallel coordinate plots; and replicate line plots. Each of these plots enable assessing the association between multiple variables.

However, these plots are ineffective with large quantities of data: Overplotting can obscure important structure, and the plots can be slow to render if every observation is mapped to a graphical element. In this project, we aim to develop more useful visualization techniques for large multivariate datasets by incorporating appropriate summaries and using interactivity. This project will explore these new visualization techniques on large RNA-sequencing datasets with plans to explore them further using large multivariate datasets from additional fields of study.

**Description:** 

To address the shortcomings of popular plotting techniques, we are focusing on expanding upon three other types of plots - the parallel coordinate plot, pairwise scatterplot matrices, and replicate line plots. There are already popular tools to visualize parallel coordinate plots and scatterplot matrices (such as, respectively, the ggparcoord() and scatmat() functions in the GGally package). However, these functions are only adequate when working with small datasets; they cannot be used with large datasets due to overplotting and time constraints. As a result, our goal is to improve upon and tailor these three types of plots so they may be useful for large multivariate data.

**Installation:**

* The latest released version: `install.packages("bigPint")`
* The latest development version: `install_github("lrutter/bigPint")`

**Resources:**

Installation of the package will automatically download a vignette, which contains a more thorough explanation of the available methods, and example code.

**License:**

GPL