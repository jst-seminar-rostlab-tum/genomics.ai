from flask import Flask, request
import os

app = Flask(__name__)


@app.route('/query', methods=['POST'])
def query():
    config = request.get_json(force=True)


@app.route("/hello")
def hello():
    return "Hello World!"


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))