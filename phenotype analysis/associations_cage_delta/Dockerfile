FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
apt-get upgrade -y && \
apt-get install libgl-dev curl -y && \
apt-get autoclean && \
apt-get autoremove && \
rm -rf /var/lib/apt/lists/*

RUN curl -L https://github.com/conda-forge/miniforge/releases/download/4.13.0-1/Miniforge3-4.13.0-1-Linux-x86_64.sh -o /tmp/miniforge.sh
RUN chmod 770 /tmp/miniforge.sh
RUN /tmp/miniforge.sh -b -p /opt/miniforge
ENV PATH=/opt/miniforge/bin:${PATH}
RUN conda init

RUN conda install jupyter scipy statsmodels optuna ipython catboost scikit-learn ipdb \
r-tidyverse r-lme4 r-emmeans r-iswr r-ggridges r-dplyr r-jtools r-broom r-irkernel \
r-dplyr r-ggplot2 r-ggthemes r-tidyr r-broom r-purrr r-plyr r-tibble r-systemfit r-ggpmisc \
r-rlang r-car r-magrittr r-minpack.lm r-scales \
r-ggdendro r-gridextra r-shiny r-miniui r-matching r-mass \
r-bitops r-rcurl r-rcppprogress -y

RUN Rscript -e 'install.packages(c("forestmangr", "MatchIt"), repos="http://cran.us.r-project.org")'

CMD jupyter notebook --allow-root

# How to use this image:
# docker build . -t cardiac_ageing
# docker run -v `pwd`:/opt/workdir -w /opt/workdir --rm -it --network=host cardiac_ageing
