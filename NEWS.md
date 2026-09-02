# funresMech 1.0.1

## Improvements and Fixes

- Resubmission replacing version 1.0.0 currently under review.
- Fixed issues affecting the saving of graphics in JPG and PDF formats within the Shiny application.
- Removed all non-ASCII characters from package code to ensure full portability.
- Added missing ggplot2 imports and declared required global variables to eliminate NOTES during R CMD check.
- Improved handling of dynamically generated variables in plotting functions.
- Updated DESCRIPTION to correctly declare all dependencies.
- Enhanced stability of the Shiny application and improved consistency of plot rendering.
- Ensured that R CMD check passes with 0 errors, 0 warnings, and 0 notes.

## Internal Enhancements

- Refactored server logic for clarity and robustness.
- Improved naming consistency for density-species interaction variables.
- Updated parallel execution handling using `future` to ensure predictable behavior.
- Minor code cleanups and formatting improvements.
