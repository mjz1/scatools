FROM zatzmanm/rstudio:latest

# install scatools
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'devtools::install_github("https://github.com/mjz1/scatools")'
