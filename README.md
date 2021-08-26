# UniformProbabilityDensityAnalysis
This is an R implementation of the Uniform Probability Density Analysis approach published by Scott Ortman (2016) and also implemented into the cyberSW web platform.

See the [cyberSW web platrorm](https://www.cybersw.org) for more.

## Description of Method

This procedure using an empirical Bayesian approach to model the population history of settlements based on their dated ceramic assemblages. Details of this approach are provided by Ortman (2016) but in brief this procedure uses the frequency of ceramic types and the associated date ranges of those types to model the probability that a site was occupied in a given period. The steps are as follows:

1) Modeling periods are first created to respresent the minimum period of overlap between type date ranges. For example, if a site has two ceramic types dated to AD 1000-1150 and 1100-1275 this would result in three modeling periods (AD 1000-1099, 1100-1149, and 1150-1275). Ceramic date assignments can also be rounded to the nearest 1...n years using a user controlled parameter to avoid very short intervals. There is also another user controlled parameter which allows the user to merge periods shorter than a selected length.

2) Next we create a prior based on the uniform distribution for each type such that each sherd has the same probability of having been deposited in any year in date range of that type and mulitply this prior distribution by the count for each type.

3) We then create a conditional probability (described in detail by Ortman 2016) as the mean probability profile across all chronologically sensitive types (as defined by the cyberSW team to exclude types like "undifferentiated gray ware" or very long lived non-specific types like "Middle Gilla Buff-ware type"). This procedure essentially has the effect of loading probability toward periods with more overlap across multiple types. For example, let's consider a hypothetical site with 4 ceramic types with date ranges of AD 1000-1100, 1000-1100, 1050-1100, and 1075-1400. The last type has a much longer range than the others that is mostly outside of the range of the other three types. This conditional procedure would load probability such that we would assume a higher probability that sherds of type 4 were deposited during the porrtion of the long interval of production that overlaps with most other date ranges at the site (in other words, the 1075-1100 period is more likely than the 1100-1400 period). 

4) We combine the prior and conditional into a posterior estimate which provides an estimate of the probability that a site was occupied in each ceramic modeling period.

5) Finally, we can trim the probability of occuption using the user controlled "cutoff" parameter which defines the minimum probability a modeling period must reach to be included in the date range for that site (see Mills et al. 2018 for an example). 

## Repository Contents

UPDA_functions.R - This file contains the primary functions for creating and plotting UPDA population curves.
Example_data.csv - This file contains an example data file formatted to work with this script. The file should include the following columns with labels exactly as they are presented here to reproduce this example. These data are in "long" format such that each row represents a unique combination of site, type, count information. 
  * Site - site name or other designation
  * Type - ceramic type name 
  * Count - count of ceramic type for the given site
  * BegDate - the beginning date associated with the ceramic type
  * EndDate - the end date associated with the ceramic type
  * Chron - this is a binary variable defined as either 0 or 1 which defines whether or not a type should be considered chronologically sensitive (1 = chronologically sensitive, 0 = not chronologically sensitive).

## References Cited:

Mills, Barbara J., Matthew A. Peeples, Leslie D. Aragon, Benjamin A Bellorado, Jeffery J. Clark, Evan Giomi, and Thomas C. Windes
2018
Evaluating Chaco Migration Scenarios Using Dynamic Social Network Analysis. Antiquity 92(364):922-939.


Ortman, Scott G.
2016
Uniform Probability Density Analysis and Population History in the Northern Rio Grande. Journal of Archaeological Method and Theory 23(1):95–126.

