import flask
from flask import request
from model_api import run_simulation
from SEIR.seir import SEIR

app = flask.Flask(__name__)

@app.route("/", methods=["GET","POST"])
def predict():
    # request.args contains all the arguments passed by our form
    # comes built in with flask. It is a dictionary of the form
    # "form name (as set in template)" (key): "string in the    
    # textbox" (value)
    print(request.args)
    if(request.args):
        results = run_simulation(request.args)
        return flask.render_template('index.html',
        datahere=results)
    else: 
         return flask.render_template('index.html')


if __name__=="__main__":
    # For local development, set to True:
    app.run(debug=True)
    # For public web serving:
    #app.run(host='0.0.0.0')
    app.run()