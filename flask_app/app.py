import flask
from flask import request
from model_api import run_simulation

app = flask.Flask(__name__)


@app.route("/", methods=["GET", "POST"])
def predict():
    print(request.args)
    if(request.args):
        results = run_simulation(request.args)
        return flask.render_template('index.html', datahere=results)
    else:
        return flask.render_template('index.html')


if __name__ == "__main__":
    # For local development, set to True:
    app.run(debug=True)
    app.run()
