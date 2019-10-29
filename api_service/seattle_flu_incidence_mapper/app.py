

# Get the application instance
from seattle_flu_incidence_mapper.config import app
from seattle_flu_incidence_mapper.models.generic_model import GenericModel
items = GenericModel.query.all()
model_ids = set([o.id for o in items])


if __name__ == "__main__":
    app.cli()
