# Dev Container for the flexaidds Python package (analysis CLI).
# The native C++ engine is optional here; Codespaces users get a working
# `flexaidds --help` without a C++26 toolchain (FLEXAIDDS_SKIP_CORE=1).
FROM mcr.microsoft.com/devcontainers/python:3.12

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      git cmake ninja-build libeigen3-dev libomp-dev g++ wget curl \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /workspaces/FlexAIDdS
