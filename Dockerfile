FROM continuumio/miniconda3
COPY . /home/
WORKDIR /home
RUN apt-get update
RUN apt install -y build-essential
RUN conda env create --file environment.linux.yml
RUN conda init bash

ENV PATH /opt/conda/envs/bambu-qsar/bin:$PATH
ENV CONDA_DEFAULT_ENV bambu-qsar

RUN /bin/bash -c "source activate bambu-qsar"
RUN echo "source ~/.bashrc"
RUN echo "conda activate bambu-qsar" > ~/.bashrc
RUN echo "source activate bambu-qsar" > ~/.bashrc

SHELL ["conda", "run", "-n", "bambu-qsar", "/bin/bash"]
CMD ["/bin/bash"]
