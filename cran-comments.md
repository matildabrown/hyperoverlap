## Test environments
* local Ubuntu 18.04 R 4.1.0
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTES:

* checking package dependencies ... NOTE
  Package suggested but not available for checking: ‘webshot2’

  webshot2 is not on CRAN, I have added Additional_repositories: https://dmurdoch.github.io/drat to DESCRIPTION as suggested by Duncan Murdoch via email
  
 * checking installed package size ... NOTE
    installed size is  6.5Mb
    sub-directories of 1Mb or more:
      doc   6.3Mb
      
This is caused by the vignette html file. 

