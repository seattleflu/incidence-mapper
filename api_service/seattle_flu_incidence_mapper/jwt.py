import time
import click
import six
import os
from werkzeug.exceptions import Unauthorized
from flask import current_app
from jose import JWTError, jwt
from seattle_flu_incidence_mapper.config import app

DEFAULT_TOKEN_DB_PATH = os.environ.get("API_TOKEN_DB", "/model_tokens")
true_values = ['1', 'y', 'yes', 't', 'true']


def generate_token(user_id):
    timestamp = _current_timestamp()
    payload = {
        "iss": current_app.config['JWT_ISSUER'],
        "iat": int(timestamp),
        "exp": int(timestamp + current_app.config['JWT_LIFETIME_SECONDS']),
        "sub": str(user_id),
    }

    return jwt.encode(payload, current_app.config['JWT_SECRET'], algorithm=current_app.config['JWT_ALGORITHM'])


@app.cli.command('generate-token')
@click.argument('user_id')
@click.argument('lifetime', default=60*60*24*364)
def generate_token_cli(user_id, lifetime):
    current_app.config['JWT_LIFETIME_SECONDS'] = lifetime
    token = generate_token(user_id)
    print(f"{token}")


def decode_token(token):
    try:
        return jwt.decode(token, current_app.config['JWT_SECRET'], algorithms=current_app.config['JWT_ALGORITHM'])
    except JWTError as e:
        six.raise_from(Unauthorized, e)


def _current_timestamp() -> int:
    return int(time.time())
