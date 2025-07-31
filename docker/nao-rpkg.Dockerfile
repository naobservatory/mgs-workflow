FROM rocker/tidyverse:latest

RUN Rscript -e 'install.packages("optparse", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("jsonlite", repos="https://cloud.r-project.org")'
