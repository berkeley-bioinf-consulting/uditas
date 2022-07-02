FROM continuumio/miniconda3

WORKDIR /inst

COPY ./uditas_env_step1.yml .

RUN conda env create -f uditas_env_step1.yml

COPY ./uditas_env_step2.yml .

RUN conda install -n uditas_env -c anaconda -c conda-forge -c bioconda --file uditas_env_step2.yml

RUN apt-get update && apt-get install -y vim

SHELL ["conda", "run", "-n", "uditas_env", "/bin/bash", "-c"]

COPY . .

RUN python setup.py install

SHELL ["/bin/bash","-c"]

ENV NB_USER uditas

ENV NB_UID 1000

RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER

WORKDIR /home/${NB_USER}

USER $NB_USER

RUN conda init

RUN echo 'conda activate uditas_env' >> ~/.bashrc
