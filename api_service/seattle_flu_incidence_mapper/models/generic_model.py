import json
from datetime import datetime

from marshmallow import pre_load, post_dump
from sqlalchemy import String, Column, DateTime
from seattle_flu_incidence_mapper.orm_config import get_session, get_declarative_base, ma

base = get_declarative_base()


class GenericModel(base):
    __tablename__ = 'generic_model'
    id = Column(String,  primary_key=True)
    name = Column(String)
    query_str = Column(String)
    model_type = Column(String)
    model_key = Column(String, primary_key=True)
    created = Column(DateTime, default=datetime.utcnow)


class GenericModelSchema(ma.ModelSchema):

    @pre_load
    def correct_query_str(self, *args):
        pass

    @post_dump
    def correct_query_str(self, data):
        data['query_str'] = json.loads(data['query_str'])
        return data

    class Meta:
        model = GenericModel
        sqla_session = get_session()
