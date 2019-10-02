import os
from flask import current_app
from seattle_flu_incidence_mapper.models.generic_model import GenericModel


def create_id_from_query_str(query_str):
    return query_str


def get_model_id_from_query_str(query_str):
    model = GenericModel.query.filter(GenericModel.query_str == query_str).first()
    return model


def get_model_file(id):
    basedir = current_app.config.get('MODEL_STORE', '/model_store')
    return os.path.join(basedir, f"{id}.csv")


def save_model_file(file, id):
    basedir = current_app.config.get('MODEL_STORE', '/model_store')
    file_path = os.path.join(basedir, id)
    csv_path = file_path + '.csv'
    try:
        file.stream.seek(0)
    except:
        pass
    file.save(csv_path)
    # convert model to json as well
    import pandas as pd
    df = pd.read_csv(csv_path)
    df.to_json(file_path + '.json', orient='records')
