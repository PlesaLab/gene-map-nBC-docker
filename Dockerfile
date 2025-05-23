# 1) Base image with Micromamba
FROM mambaorg/micromamba:1.5.6

# 2) Switch to root to install system compilers & libs needed by edlib/pyabpoa
USER root
#  - build-essential brings in gcc/g++; cmake is needed by pyabpoa;
#  - libssl-dev, zlib1g-dev cover common header deps
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        wget unzip build-essential cmake libssl-dev zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

 # 3) Use bash for RUN steps (so micromamba activation hooks work)
SHELL ["/bin/bash", "-lc"]

# 4) Copy in your environment spec
COPY newenv.yml /tmp/env.yml

# 5) Create the newenv conda env (this will now have gcc available for the pip section)
RUN micromamba create -y -n newenv -f /tmp/env.yml \
 && micromamba clean --all --yes \
 && rm /tmp/env.yml

 # 6) Auto-activate the newenv env
ENV CONDA_DEFAULT_ENV=newenv \
    MAMBA_DOCKERFILE_ACTIVATE=1 \
    PATH=/opt/conda/envs/newenv/bin:$PATH

    # 7) (Optional) add a non-root user matching host UID/GID for volume permissions
ARG USER_ID=1000
ARG GROUP_ID=1000
RUN if ! getent group "$GROUP_ID" >/dev/null; then addgroup --gid "$GROUP_ID" mambauser; fi \
 && usermod --uid $USER_ID --gid $GROUP_ID mambauser

 # 8) Switch to the non-root user
USER mambauser
WORKDIR /home/mambauser

# 9) Default to an interactive shell with newenv auto-activated
CMD ["bash"]