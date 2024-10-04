


# Use bioconductor/bioconductor_docker for linux/amd64 as base image
FROM --platform=linux/amd64 bioconductor/bioconductor_docker:RELEASE_3_17


MAINTAINER Gabriel Villamil <gabriel.villamil@mdc-berlin.de>


# Execute by default (needed when running with Singularity)
CMD /bin/bash "$@"


# Use bash as shell
SHELL ["/bin/bash"]