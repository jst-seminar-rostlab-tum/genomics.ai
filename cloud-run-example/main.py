from os import environ
from time import sleep
from threading import Thread
from sys import exit
from pymongo import MongoClient
from flask import Flask, request

app = Flask(__name__)
db = MongoClient(environ["DATABASE_URI"]).get_default_database()

@app.route("/")
def hello_world():
    try:    
        upload_id = request.args.get('upload_id')

        project = db.projects.find_one({"uploadId": upload_id})

        if project is None:
            return f"There exists no project with upload ID {upload_id}"

        thread = Thread(target= do_work, args=(upload_id, ))
        thread.start()

        return f"Started working on {upload_id}..."
    except:
        return f"Please provide a valid bucket ID!"

def do_work(upload_id: str) -> None:
    for _ in range(5):
        project = db.projects.find_one({"uploadId": upload_id})

        if project is None:
            print("Project does not exist anymore. Terminating")
            exit()

        if project.status == 3:
            print("Project has been aborted. Terminating.")
            exit()

        sleep(3)

    db.projects.update_one({'uploadId': upload_id}, {"$set": {"status": 4 }})

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(environ.get("PORT", 8080)))