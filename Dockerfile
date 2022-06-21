FROM continuumio/miniconda3

WORKDIR /inst

# RUN git clone --depth=1 https://github.com/editasmedicine/uditas 
COPY ./uditas_env_step1.yml .

RUN conda env create -f uditas_env_step1.yml

COPY ./uditas_env_step2.yml .

RUN conda install -n uditas_env -c anaconda -c conda-forge -c bioconda --file uditas_env_step2.yml

SHELL ["conda", "run", "-n", "uditas_env", "/bin/bash", "-c"]

COPY . .

RUN python setup.py install

CMD ["conda", "run", "-n", "uditas_env", "/bin/bash"]
