# FROM ghcr.io/feelpp/feelpp:v0.110.2-focal
FROM ghcr.io/feelpp/feelpp:v0.111.0-alpha.5-bullseye

USER root

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    LC_ALL=C.UTF-8 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

# This Dockerfile adds a non-root user with sudo access. Use the "remoteUser"
# property in devcontainer.json to use it. On Linux, the container user's GID/UIDs
# will be updated to match your local UID/GID (when using the dockerFile property).
# See https://aka.ms/vscode-remote/containers/non-root-user for details.
ARG USERNAME=feelpp
ARG USER_UID=1001
ARG USER_GID=$USER_UID

# for singularity
ARG SVERSION=3.7.4
ARG GVERSION=1.14.12
ARG OS=linux
ARG ARCH=amd64

# install Feelpp from BinTray Debian repository
RUN apt update \
    && apt install -y lsb-release \
    && apt install -y debian-keyring \
    && cp /usr/share/keyrings/debian-maintainers.gpg /etc/apt/trusted.gpg.d \
    # Configure apt and install packages \
    && apt-get update \
    && apt-get -y install --no-install-recommends apt-utils dialog 2>&1 \
    # \
    # Verify git, process tools, lsb-release (useful for CLI installs) installed\
    && apt-get -y install git iproute2 procps lsb-release iputils-ping nmap \
    #\
    # Install C++ tools\
    && apt-get -y install build-essential python3-minimal libpython3-dev \
    && apt-get -y install python-is-python3 python3-pip python3-venv python3-jinja2 \
    && apt-get -y install python3-autopep8 black yapf3 python3-bandit flake8 pydocstyle pylint python3-pytest mypy \
    && apt-get -y install python3-statsmodels \
    && apt-get -y install python3-matplotlib texlive-latex-base texlive-latex-extra dvipng \
    && apt-get -y install python3-pandas python3-seaborn \
    && apt-get -y install python3-numpy python3-scipy \
    && apt-get -y install python3-nlopt python3-tabulate \
    && apt-get -y install python3-pip python-is-python3 \
    && apt-get -y install python3-pint \
    && apt-get -y install curl \
    && apt-get -y install python3-decouple python3-requests python3-natsort \
    # Add specific lncmi stuff \
    && apt-get install python3-iapws \
    # [Optional] Add sudo support for the non-root user
    && apt-get install -y sudo \
    #&& if [ "$(lsb_release -cs)" = "focal" ]; then echo "Set disable_coredump false" > /etc/sudo.conf; fi \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME\
    && chmod 0440 /etc/sudoers.d/$USERNAME \
    # add github ssh key
    && mkdir /home/$USERNAME/.ssh/ \
    && ssh-keyscan github.com >> /home/$USERNAME/.ssh/known_hosts \
    && chown -R $USERNAME.$USER_GID /home/$USERNAME/.ssh \
    #
    # Clean up
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog

# Setup demo environment variables
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    LC_ALL=C.UTF-8 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

# install visidata in virtual env for $USERNAME
