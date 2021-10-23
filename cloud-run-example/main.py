import os
import time
import random
from threading import Thread

from pymongo import MongoClient
from flask import Flask, request

app = Flask(__name__)


@app.route("/")
def hello_world():

    try:    
        bucket_id = int(request.args.get('bucket_id'))

        thread = Thread(target= do_work, args=(bucket_id, ))
        thread.start()

        return f"Started working on {bucket_id}..."
    except:
        return f"Please provide a valid bucket ID!"

def do_work(bucket_id: int) -> None:
    time.sleep(random.random() * 5)
    db = MongoClient(os.environ["DATABASE_URI"]).get_default_database()

    db.jobs.update_one({'bucketId': bucket_id}, {"$set": {"isFinished": True }})

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))