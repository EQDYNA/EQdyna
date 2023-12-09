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
docker run -it --name ContainerName dunyuliu/eqdyna:v5.3.1
```
Note: for Windows users, you may need admin access for the Powershell. For MacOS users, it seems you don't need admin access (to be explored further). <br/>
  
## 3. Using the newly created EQdyna docker container! <br/>

a. After downloading the image dunyuliu/eqdyna:v5.3.1 (<1.5GB), you will find it in the 'Images' tab on the left control panel of Docker Desktop. <br/>
b. Also you will find a container with the name you give, i.e., ContainerName, running in the 'Containers' tab in Docker Desktop. <br/>
c. In the 'Actions' panel of the container, if you click the three dots, you will find an drop-down list. <br/>
       Click 'Open in terminal' and you will be navigated to the terminal of the Ubuntu system. <br/>
d. Type 'bash', you will enter the $HOME directory where you can find useful information in the 0README.md. <br/>
e. Now hope you can enjoy running the EQdyna! <br/>
    
## Appendix:
1. To commit a modified/developed container to a Docker image and push it to the Docker Hub.
```
docker commit ContainerName your_docker_hub_name/image_name:tag
```

2. To mount the newly created docker container to other hard drives. Additional to the above docker run command, you can specify a hard drive that is mounted to your system for additional space. It will be achieved by adding the -v option. 
	
For example, in the following command, I mount the container to an external Windows OS drive G:/scratch, and rename it to /mount in my container environment. 
	
```
docker run -it --name ContainerName -v G:/scratch:/mount dunyuliu/eqdyna:v5.3.1
	
```
