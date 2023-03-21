# Docker Guide

## News! EQdyna is now available through Docker.

To use the EQdyna docker image, you need three steps.  <br/>

## 1. Install Docker Desktop. <br/>
a. For Windows users, you can download Docker Desktop from this Link [Download Docker Desktop](https://www.docker.com/products/docker-desktop/). <br/>
       After the installation, you may need to run it as the administrator.  
b. For MacOS users, you can also download Docker Desktop from the same link but to choose the Mac option - either Apple Chip or Intel Chip.<br/>
c. For Linux users, to be explored. <br/>

## 2. Pull the docker image from the Docker Hub. <br/>

Open a terminal and run the following command: <br/>
```
docker run -it --name eqdyna dunyuliu/eqdyna.docker
```
Note: for Windows users, you may need admin access for the Powershell. For MacOS users, it seems you don't need admin access (to be explored further). <br/>
  
## 3. Using the newly created EQdyna docker container! <br/>

a. After some downloading (1.14 GB), you will find a new image called dunyuliu/eqdyna.docker in the 'Images' tab in the Docker Desktop. <br/>
b. Also you will find a container with the name you give - eqdyna in the command above - running in the 'Containers' tab in the Docker Desktop. <br/>
c. In the 'Actions' panel of the container, if you click the three dots, you will find an drop-down list. <br/>
       Click 'Open in terminal' and you will be navigated to the terminal of the Ubuntu system. <br/>
d. Type 'bash', you will enter the $HOME directory where you can find useful information in the 0README.md. <br/>
e. Now you can enjoy running the EQdyna! <br/>
    
## Appendix:
1. To commit a modified/developed container to a Docker image and push it to the Docker Hub.
```
docker commit container_name your_docker_hub_name/image_name:tag
```

2. To mount the newly created docker container to other hard drives. Additional to the above docker run command, you can specify a hard drive that is mounted to your system for additional space. It will be achieved by adding the -v option. 
	
For example, in the following command, I mount the container to an external drive G:/docker.data and rename it to /mount shown in the container environment. 
	
```
docker run -it --name eqdyna.docker -v G:/docker.data:/mount dunyuliu/eqdyna.docker
	
```
