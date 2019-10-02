import argparse
import json
import os

import requests


def download_csv(base_url, model):
    headers = dict(accepts='text/csv')
    query = json.loads(model['query_str'])
    response = requests.post(f'{base_url}/v1/query', json=query, headers=headers)
    if response.status_code != 200 or len(response.text) < 20:
        print(f'Model {model["id"]} failed')


def get_models(url, model_store, model_type=None, only_latest=True):
    response = requests.get(url)
    if response.status_code == 200:
        models = response.json()
        if only_latest:
            n_models = {}
            for model in models:
                if model['query_str'] not in n_models:
                    n_models[model['query_str']] = model
            models = n_models.values()
        for model in models:
            query_obj = json.loads(model['query_str'])
            if model_type is None or query_obj['model_type'] == model_type:
                if not os.path.exists(os.path.join(model_store, f"{model['id']}.json")):
                    metafn = os.path.join(model_store, f"{model['id']}.json")
                    json.dump(model, open(metafn, 'w'))

            #response = model['query_str']


if __name__ == "__main__":

    default_model_store_path = os.path.abspath(os.path.join(os.path.os.getcwd(), '../../test_model_store'))
    parser = argparse.ArgumentParser(description='Uploads trained SF Models to production')
    parser.add_argument("--db-file", default=os.path.join(default_model_store_path, "modelDB.tsv"),
                        help="Where the modelDB.tsv produced during training is stored")
    parser.add_argument("--model-store", default=default_model_store_path)
    parser.add_argument("--api-url", default="http://40.112.165.255/v1/generic_models",
                        help="URL for Seattle FLU API Incidence Mapper Model Server API")
    parser.add_argument("--model-type", default=None, help="Filter models by model_type")
    parser.add_argument("--no-only-latest", dest='latest', action='store_false')
    parser.set_defaults(latest=True)
    args = parser.parse_args()
    get_models(args.api_url, args.model_store, args.model_type, args.latest)


