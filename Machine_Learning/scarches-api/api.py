from flask import Flask, request
import os
import init as scarches
app = Flask(__name__)


@app.route('/query', methods=['POST'])
def query():
    config = request.get_json(force=True)
    scarches.query(config)
    return "computed", 200


@app.route("/liveness")
def hello():
    return "up"


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))