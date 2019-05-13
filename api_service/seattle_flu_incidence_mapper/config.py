import json
import os
import connexion
from flask import Response
from sqlalchemy.orm.exc import NoResultFound

from seattle_flu_incidence_mapper.orm_config import setup_db
from seattle_flu_incidence_mapper.utils import set_marshmallow, ModelExecutionException

basedir = os.path.abspath(os.path.dirname(__file__))
true_vals = ['1', 'y', 'yes', 't', 'true']
# Create the Connexion application instance
connex_app = connexion.App("seattle_flu_incidence_mapper.config", specification_dir=os.path.join(basedir, 'swagger'))

# Get the underlying Flask app instance
app = connex_app.app
app.config['WORKER_IMAGE'] = os.environ.get('WORKER_IMAGE', 'idm-docker-production.packages.idmod.org/sfim-worker:latest')
app.config['MODEL_STORE'] = os.environ.get('MODEL_STORE', os.path.abspath(os.path.join(os.getcwd(), "../../test_model_store")))
app.config['MODEL_HOST_PATH'] = os.environ.get('MODEL_HOST_PATH',  os.path.abspath(os.path.join(os.getcwd(), "../../test_model_store")))
app.config['WORKER_JOB_HOST_PATH'] = os.environ.get('WORKER_JOB_DIR_HOST_PATH',  os.path.abspath(os.path.join(os.getcwd(), "../../test_jobs")))
app.config['MODEL_JOB_PATH'] = os.environ.get('MODEL_JOB_PATH',  os.path.abspath(os.path.join(os.getcwd(), "../../test_jobs")))

db = setup_db(basedir, app)
set_marshmallow(app)

# DO NOT MOVE this line. The order matters here
# we need to init our db before loading our models
from seattle_flu_incidence_mapper.models import *

if os.environ.get('DEBUG', '0').lower() in true_vals or os.environ.get('CREATE_DB', '0').lower() in true_vals:
    db.create_all()


def sqlalchemy_error_handler(exception):
    return Response(response=json.dumps({'error': ' '.join(exception.args)}),
                    status=404, mimetype="application/json")


def model_exec_error_handler(exception):
    return Response(response=json.dumps({'error': ' '.join(exception.title)}),
                    status=500, mimetype="application/json")


connex_app.add_error_handler(ModelExecutionException, model_exec_error_handler)
connex_app.add_error_handler(NoResultFound, sqlalchemy_error_handler)