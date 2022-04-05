# GET PDAL 2.3.0
FROM daunnc/pdal:2.3.0 as pdal

# PYTHON BASE IMAGE
FROM python:3.6 as py

RUN apt-get update && apt-get install -y git && \
    apt install software-properties-common ffmpeg libsm6 libxext6 -y 

# DEBUG INFO
RUN set -ex && \
    apt-get install libgdal-dev -y \
    && /usr/local/bin/python -m pip install --upgrade pip

# INSTALL GDAL
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal
RUN pip3 install GDAL==3.1.4

# INSTALL PDAL - WRONG VERSION ! (2.3.0 REQUIRED)
# RUN apt-get update && apt-get install pdal -y

WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt

# COPY PDAL 2.3.0 TO CURRENT LAYER AND ADD TO PATH
COPY --from=pdal /opt/conda/envs/pdal /opt/conda/envs/pdal
ENV PATH=/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/conda/envs/pdal/bin

ENTRYPOINT ["/bin/bash"]
