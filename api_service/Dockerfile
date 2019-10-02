FROM python:3.7-stretch
RUN apt update && apt install nginx \
        python-dev \
        apt-transport-https \
        ca-certificates \
        curl \
        gnupg2 \
        software-properties-common \
        -y
    # Make project directory
RUN mkdir -p /app/seattle_flu_incidence_mapper \
        /model_store

ENV PYTHONPATH=/app:${PYTHONPATH}
ENV PATH=/app:${PATH}
ENV FLASK_APP=seattle_flu_incidence_mapper.app:app
WORKDIR /app

ENV MODEL_STORE=/model_store

ADD setup.py README.rst /app/seattle_flu_incidence_mapper/
RUN cd seattle_flu_incidence_mapper && \
    pip3 install -e .[production]

# Define we will have a model data folder
VOLUME /model_data
# Finish the setup
ADD docker_scripts/nginx.conf /etc/nginx/nginx.conf
ADD docker_scripts/entrypoint.sh docker_scripts/wsgi.ini ./
# Add migrations
ADD migrations ./migrations
# Add the rest of our library here
ADD seattle_flu_incidence_mapper ./seattle_flu_incidence_mapper

EXPOSE 80
#VOLUME [/model_store]

CMD /app/entrypoint.sh