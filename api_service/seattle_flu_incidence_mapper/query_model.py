# API Methods for the /query
import hashlib
import os
import tarfile
import time
from io import BytesIO
import docker
from flask import current_app, send_file, request
from sqlalchemy.orm.exc import NoResultFound
import json
from seattle_flu_incidence_mapper.models.generic_model import GenericModel
from seattle_flu_incidence_mapper.utils import get_model_id, ModelExecutionException

loaded_models = []
client = docker.DockerClient()
api_client = docker.APIClient()


def query():
    query_json = request.json
    file_format = 'csv' if 'csv' in request.headers.get('accept', 'json').lower() else 'json'

    model_id = get_model_id(query_json)
    model_path = os.path.join(current_app.config['WORKER_IMAGE'], model_id + '.' + file_format)

    return send_file(model_path,
                     as_attachment=False,
                     mimetype='application/json' if file_format == 'json' else 'text/csv'
                     )

