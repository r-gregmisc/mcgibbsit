## Current Submission (1.2.2):

-   Local R CMD check results

    0 errors \| 0 warnings \| 0 note

-   This is a resubmission of `mcgibbsit` with back quotes(\`\`) changed to single quotes ('') in the DESCRIPTION title field as requested by Uwe Ligges.

## Previous Submission (1.2.1):

-   Local R CMD check results

    0 errors \| 0 warnings \| 0 note

-   R-Hub builder R CMD check results:

    -   Two NOTEs are present under `R devel` on all three of Windows, MacOS, and Ubuntu Linux

        1.  `Package was archived on CRAN`: This submission resolves all of the issues listed on the the CRAN package check page (<https://cran-archive.r-project.org/web/checks/2022/2022-06-23_check_results_mcgibbsit.html>).

        2.  `Possibly misspelled words in DESCRIPTION`: All of these are author or method names.

    -   The R-Hub builder reports some NOTEs that appear to be issues with the particular R installation, rather than with the `mcgibbsit` package itself.

        -   Windows:
        
            ```         
            * checking for non-standard things in the check directory ... NOTE
            Skipping checking math rendering: package 'V8' unavailable
            Found the following files/directories:
              ''NULL''
            * checking for detritus in the temp directory ... NOTE
            Found the following files/directories:
              'lastMiKTeXException'
            ```

        -   MacOs and Ubuntu:

            ```         
            Found the following (possibly) invalid URLs:
              URL: https://stat.uw.edu/sites/default/files/files/reports/2001/tr395.pdf
                From: man/mcgibbsit.Rd
                Status: Error
                Message: libcurl error code 60:
                    SSL certificate problem: unable to get local issuer certificate
                    (Status without verification: OK)
                    
            * checking HTML version of manual ... NOTE
            Skipping checking HTML validation: no command 'tidy' found
            Skipping checking math rendering: package 'V8' unavailable
            ```
