from flask import Flask, request
import os
import init as scarches
from threading import Thread
from utils import utils, parameters
import traceback

app = Flask(__name__)


def get_from_config(configuration, key):
    if key in configuration:
        return configuration[key]
    return None


@app.route('/query', methods=['POST'])
def query():
    try:
        config = request.get_json(force=True)
        run_async = get_from_config(config, parameters.RUN_ASYNCHRONOUSLY)
        if run_async is not None and run_async:
            actual_config = scarches.merge_configs(config)
            thread = Thread(target=scarches.query, args=(config,))
            thread.start()
            return actual_config, 200
        else:
            actual_configuration = scarches.query(config)
            return actual_configuration, 200
    except Exception as e:
        traceback.print_exc()
        return {'error': str(e)}, 500


@app.route("/liveness")
def liveness():
    return "up"


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))
