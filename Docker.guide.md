# Docker Guide

## News! EQdyna is now available through Docker.

To use the EQdyna docker image, you need three steps.  <br/>

## 1. Install Docker Desktop. <br/>
a. For Windows users, you can download Docker Desktop from this Link [Download Docker Desktop](https://www.docker.com/products/docker-desktop/). <br/>
       After the installation, you may need to run it as the administrator.  
b. For MacOS users, you can also download Docker Desktop from the same link but to choose the Mac option - either Apple Chip or Intel Chip.<br/>
c. For Linux users, to be explored. <br/>

## 2. Pull the docker image from GitHub Container Registry. <br/>

Open a terminal and run the following command: <br/>
```
docker run -it --name ContainerName ghcr.io/eqdyna/eqdyna:latest
```
The image is built directly from the tagged commit by CI
(`.github/workflows/publish.yml`, triggered on `push: tags: ['v*']`) -- there
is no hand-modified or hand-committed step between the repository and the
image, so `:latest` always tracks the most recent tag CI has built. To pin a
specific released version instead (for reproducibility), replace `:latest`
with `:vX.Y.Z` for any tag on the
[Releases page](https://github.com/EQDYNA/EQdyna/releases) -- do not hardcode
a version number in this guide itself, since a pinned example here is exactly
what let a past version drift seven minor releases out of date before anyone
noticed (`testsys/regression/test_docker_guide_no_pinned_version.py` guards
against that recurring).

Note: for Windows users, you may need admin access for the Powershell. For MacOS users, it seems you don't need admin access (to be explored further). <br/>

## 3. Using the newly created EQdyna docker container! <br/>

a. After pulling the image `ghcr.io/eqdyna/eqdyna:latest` (or a pinned `:vX.Y.Z`), you will find it in the 'Images' tab on the left control panel of Docker Desktop. <br/>
b. Also you will find a container with the name you give, i.e., ContainerName, running in the 'Containers' tab in Docker Desktop. <br/>
c. In the 'Actions' panel of the container, if you click the three dots, you will find an drop-down list. <br/>
       Click 'Open in terminal' and you will be navigated to the terminal of the Ubuntu system. <br/>
d. Type 'bash', you will land in `/opt/eqdyna` with `EQDYNAROOT`/`PATH` already set (see the root `Dockerfile`). <br/>
e. Now hope you can enjoy running the EQdyna! <br/>

## Two caveats

* **MPI inside the container is single-node only.** The image is fine for
  running the gated 4-rank cases on the cores available to the container, but
  it does not give you a cluster -- there is no multi-node MPI fabric across
  containers here. For a cluster (e.g. TACC Lonestar6), use
  `./install-eqdyna.sh -m ls6` directly on the cluster, not this image.
* **This image is not a citable artifact.** It answers "I cannot get netCDF to
  build", not "how do I cite this run of EQdyna". Citability is a separate
  Zenodo DOI + `CITATION.cff` effort that does not exist yet; do not treat
  pulling this image as equivalent to archiving a citable version of the code.

## Mounting additional storage

To mount the newly created docker container to other hard drives, additional to the above docker run command, you can specify a hard drive that is mounted to your system for additional space. It will be achieved by adding the -v option.

For example, in the following command, I mount the container to an external Windows OS drive G:/scratch, and rename it to /mount in my container environment.

```
docker run -it --name ContainerName -v G:/scratch:/mount ghcr.io/eqdyna/eqdyna:latest
```
