FROM zatzmanm/rstudio:v0.2

# install scatools
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'devtools::install_github("https://github.com/mjz1/scatools", quiet = TRUE)'