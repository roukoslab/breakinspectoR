# shiny-server / CD 26062024
# Uses bioconductor image as base
FROM docker.io/bioconductor/bioconductor_docker:RELEASE_3_19

# Set noninteractive to avoid prompts during build
ENV DEBIAN_FRONTEND noninteractive

# Update system and install packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils gdebi-core && \
    apt-get full-upgrade -y && \
    curl "https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-$(curl https://download3.rstudio.org/ubuntu-18.04/x86_64/VERSION)-amd64.deb" -o ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    apt-get clean && \
    rm -f ss-latest.deb && \
    rm -rf /var/lib/apt/lists/* /etc/services.d/rstudio /etc/init.d/rstudio-server /etc/cont-init.d

# Install bluntPred package and all its dependencies
RUN R --slave -e 'install.packages(c("RCurl","jsonlite", "devtools", "shiny", "shinydashboard", "shinyjs", "ggplot2", "cowplot"))' && \
    R --slave -e 'install.packages("h2o", type="source", repos="http://h2o-release.s3.amazonaws.com/h2o/rel-zumbo/2/R")' && \
    R --slave -e 'devtools::install_github("omarwagih/ggseqlogo")' && \
    R --slave -e 'devtools::install_github("roukoslab/breakinspectoR")' && \
    R --slave -e 'file.copy(system.file("shiny_BTmotif", package="breakinspectoR"), "/srv/shiny-server/", recursive=TRUE)' && \
    R --slave -e 'file.copy(system.file("shiny_bluntPred", package="breakinspectoR"), "/srv/shiny-server/", recursive=TRUE)' && \
    sed -i 's|#{{Impressum-placeholder}}|, p(style="text-align: right;", a("the Institute of Molecular Biology (IMB) in Mainz", href="https://imb.de/", target="_blank")), p(style="text-align: right;", a("Impressum - Imprint", href="https://imb.de/impressum-imprint", target="_blank"))|' /srv/shiny-server/shiny_BTmotif/app.R && \
    sed -i 's|#{{Impressum-placeholder}}|, p(style="text-align: right;", a("the Institute of Molecular Biology (IMB) in Mainz", href="https://imb.de/", target="_blank")), p(style="text-align: right;", a("Impressum - Imprint", href="https://imb.de/impressum-imprint", target="_blank"))|' /srv/shiny-server/shiny_bluntPred/app.R

# Setup permissions
RUN chown -R shiny:shiny /var/lib/shiny-server && \
    sed -i 's/directory_index on/directory_index off/' /etc/shiny-server/shiny-server.conf && \
    rm -rf /srv/shiny-server/index.html /srv/shiny-server/sample-apps

# Set user and run shiny-server
USER shiny
CMD ["shiny-server"]

